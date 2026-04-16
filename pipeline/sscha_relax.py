"""
SSCHA anharmonic phonon renormalization for the QE + EPW pipeline.

Usage:
    python sscha_relax.py compounds/CaH6_200GPa.yaml [--supercell 2 2 2]
                                                      [--nconfigs 50]
                                                      [--niter 20]
                                                      [--temp 300]
                                                      [--start-from-iter N]

What it does:
    1. Loads harmonic dynamical matrices from runs/<prefix>/2_phonons/<prefix>.dyn*
    2. Makes them positive-definite (replaces imaginary modes with small positive freqs)
       so SSCHA can use them as a starting trial force constants
    3. Runs SSCHA minimization: generates random structures, computes QE forces,
       updates trial force constants until convergence (dF < threshold)
    4. Writes renormalized dyn files to runs/<prefix>/sscha/<prefix>_sscha.dyn*
    5. These can be used in place of harmonic dyns for EPW (re-run from phonons step
       using --start-from phonons after symlinking or copying the sscha dyn files)

Notes:
    - Each SSCHA iteration runs nconfigs DFT single-point force calculations.
      With nconfigs=50 and ~2 min/config, one iteration takes ~100 min on Mac.
      Use nconfigs=10 for PoC, 50-100 for production.
    - SSCHA minimizes the free energy at a given temperature. Run at several
      temperatures to check T-dependence of phonon renormalization.
    - For CaH6 at 200 GPa: expect imaginary modes to become real after ~5-10 iter.
    - The output dyn files can feed directly back into the EPW pipeline.
"""

import sys
import os
import argparse
import yaml
import subprocess
import shutil
import numpy as np
from pathlib import Path

import cellconstructor as CC
import cellconstructor.Phonons as Phonons
import cellconstructor.Structure as Struct
import cellconstructor.Units as Units
import sscha
import sscha.Ensemble as Ensemble
import sscha.Relax as Relax


# ─────────────────────────────────────────────────────────────────────────────
# QE force calculator
# ─────────────────────────────────────────────────────────────────────────────

def _write_scf_input(structure, cfg, scf_path, pseudo_dir, pressure_kbar):
    """Write a QE pw.x input for a single-point force+stress calculation."""
    p = cfg["compound"]
    d = cfg["dft"]
    s_cfg = cfg["structure"]

    # Build ATOMIC_SPECIES
    species_lines = []
    for sp in s_cfg["species"]:
        species_lines.append(f"  {sp['symbol']:<4s}  {sp['mass']:.5f}  {sp['pseudo']}")

    # Build ATOMIC_POSITIONS from CellConstructor structure
    pos_lines = []
    for i, (sym, pos) in enumerate(zip(structure.atoms, structure.coords)):
        pos_lines.append(f"  {sym:<4s}  {pos[0]:.9f}  {pos[1]:.9f}  {pos[2]:.9f}")

    # Build CELL_PARAMETERS (angstrom)
    cell = structure.unit_cell  # in angstrom
    cell_lines = [f"  {row[0]:14.9f}  {row[1]:14.9f}  {row[2]:14.9f}" for row in cell]

    content = f"""\
&CONTROL
  calculation   = 'scf'
  restart_mode  = 'from_scratch'
  prefix        = '{p['prefix']}_sscha'
  outdir        = './tmp/'
  pseudo_dir    = '{pseudo_dir}'
  tprnfor       = .true.
  tstress       = .true.
/

&SYSTEM
  ibrav    = 0
  nat      = {structure.N_atoms}
  ntyp     = {len(s_cfg['species'])}
  input_dft = 'PBE'
  ecutwfc  = {d['ecutwfc']}
  ecutrho  = {d['ecutrho']}
  occupations = 'smearing'
  smearing    = 'mp'
  degauss     = {d['degauss']}
  press       = {pressure_kbar:.2f}
/

&ELECTRONS
  conv_thr    = 1.0d-8
  mixing_beta = 0.4
  diagonalization = 'david'
/

&IONS
/

&CELL
  cell_dynamics = 'none'
/

ATOMIC_SPECIES
{chr(10).join(species_lines)}

CELL_PARAMETERS {{angstrom}}
{chr(10).join(cell_lines)}

ATOMIC_POSITIONS {{angstrom}}
{chr(10).join(pos_lines)}

K_POINTS {{automatic}}
{d['scf_kgrid'][0]} {d['scf_kgrid'][1]} {d['scf_kgrid'][2]} 0 0 0
"""
    scf_path.write_text(content)


def _parse_qe_output(out_path, structure):
    """Parse forces and energy from pw.x output. Returns (energy_Ry, forces_Ry/bohr)."""
    text = out_path.read_text()

    # Energy
    import re
    m = re.search(r'!\s+total energy\s+=\s+([\d.+-]+)\s+Ry', text)
    if not m:
        raise RuntimeError(f"Could not find energy in {out_path}")
    energy_Ry = float(m.group(1))

    # Forces (in Ry/bohr)
    forces = np.zeros((structure.N_atoms, 3))
    force_block = re.search(
        r'Forces acting on atoms.*?atom.*?force\s*=\s*(.*?)(?=Total|$)',
        text, re.DOTALL
    )
    if force_block:
        for line in force_block.group(0).split('\n'):
            m2 = re.match(r'\s+atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)', line)
            if m2:
                idx = int(m2.group(1)) - 1
                forces[idx] = [float(m2.group(2)), float(m2.group(3)), float(m2.group(4))]
    else:
        raise RuntimeError(f"Could not find forces in {out_path}")

    return energy_Ry, forces


class QEForceCalculator:
    """Computes forces on an ensemble of structures using pw.x."""

    def __init__(self, cfg, run_dir, qe_bin, pseudo_dir):
        self.cfg = cfg
        self.run_dir = Path(run_dir)
        self.qe_bin = qe_bin
        self.pseudo_dir = pseudo_dir
        self.pressure_kbar = cfg["compound"]["pressure_GPa"] * 10.0

    def compute_ensemble(self, ensemble, iter_idx):
        """Run pw.x on all structures in ensemble, fill energies and forces."""
        work_dir = self.run_dir / f"iter_{iter_idx:03d}"
        work_dir.mkdir(exist_ok=True)

        energies = np.zeros(ensemble.N)
        forces = np.zeros((ensemble.N, ensemble.structures[0].N_atoms, 3))

        for i, structure in enumerate(ensemble.structures):
            conf_dir = work_dir / f"conf_{i:04d}"
            conf_dir.mkdir(exist_ok=True)
            (conf_dir / "tmp").mkdir(exist_ok=True)

            scf_in = conf_dir / "scf.in"
            scf_out = conf_dir / "scf.out"

            _write_scf_input(structure, self.cfg, scf_in, self.pseudo_dir, self.pressure_kbar)

            if scf_out.exists():
                try:
                    e, f = _parse_qe_output(scf_out, structure)
                    energies[i] = e
                    forces[i] = f
                    print(f"  conf {i:4d}/{ensemble.N}: cached", flush=True)
                    continue
                except Exception:
                    pass

            print(f"  conf {i:4d}/{ensemble.N}: running pw.x...", end=" ", flush=True)
            result = subprocess.run(
                [str(Path(self.qe_bin) / "pw.x"), "-in", str(scf_in)],
                stdout=open(scf_out, "w"),
                stderr=subprocess.STDOUT,
                cwd=str(conf_dir)
            )
            if result.returncode != 0:
                print("FAILED")
                raise RuntimeError(f"pw.x failed for conf {i}, check {scf_out}")

            e, f = _parse_qe_output(scf_out, structure)
            energies[i] = e
            forces[i] = f
            print(f"E={e:.6f} Ry", flush=True)

        # Fill ensemble
        ensemble.energies = energies * Units.RY_TO_EV
        ensemble.forces = forces * (Units.RY_TO_EV / Units.BOHR_TO_ANGSTROM)
        ensemble.has_stress = False


# ─────────────────────────────────────────────────────────────────────────────
# Main SSCHA loop
# ─────────────────────────────────────────────────────────────────────────────

def run_sscha(yaml_path, supercell, nconfigs, niter, temperature, start_from_iter):
    cfg = yaml.safe_load(open(yaml_path))
    prefix = cfg["compound"]["prefix"]
    qe_bin = cfg["paths"]["qe_bin"]
    pseudo_dir = cfg["paths"]["pseudo_dir"]

    pipeline_dir = Path(yaml_path).parent.parent
    phonon_dir = pipeline_dir / "runs" / prefix / "2_phonons"
    sscha_dir = pipeline_dir / "runs" / prefix / "sscha"
    sscha_dir.mkdir(exist_ok=True)

    # ── Load harmonic dynamical matrices ──────────────────────────────────────
    dyn_prefix = str(phonon_dir / f"{prefix}.dyn")
    # Count available dyn files (dyn0 is metadata, dyn1..N are q-points)
    nqirr = sum(1 for p in phonon_dir.glob(f"{prefix}.dyn[0-9]*") if p.stat().st_size > 100) - 1

    # Read q-grid from dyn0 (first line: nq1 nq2 nq3)
    dyn0_tokens = (phonon_dir / f"{prefix}.dyn0").read_text().split()
    qgrid = [int(dyn0_tokens[0]), int(dyn0_tokens[1]), int(dyn0_tokens[2])]
    print(f"\nLoading harmonic dyn matrices: {dyn_prefix}1 ... {dyn_prefix}{nqirr}")
    dyn = Phonons.Phonons(dyn_prefix, nqirr=nqirr)
    print(f"  {dyn.structure.N_atoms} atoms, {nqirr} irreducible q-points, q-grid {qgrid}")

    # Check for imaginary modes
    freqs_all = []
    for iq in range(dyn.nqirr):
        w2, _ = np.linalg.eigh(dyn.dynmats[iq])
        freqs_all.extend(np.sign(w2) * np.sqrt(np.abs(w2)))
    n_imag = sum(1 for f in freqs_all if f < -1e-3)
    print(f"  Imaginary modes in harmonic: {n_imag}")

    # ── Make positive-definite as SSCHA starting point ───────────────────────
    dyn_sscha = dyn.Copy()
    dyn_sscha.ForcePositiveDefinite()
    print("  Applied ForcePositiveDefinite (replaces imaginary modes with small positive)")

    # ── Build supercell dynamical matrix ────────────────────────────────────
    # Supercell must match the phonon q-grid — that's what defines the real-space FC range.
    # User --supercell arg is ignored; the q-grid from dyn0 is authoritative.
    sc = qgrid
    if supercell != [2, 2, 2]:
        print(f"  Note: --supercell {supercell} ignored; using q-grid {sc} from dyn0")
    print(f"\nBuilding {sc[0]}x{sc[1]}x{sc[2]} supercell dynamical matrix...")
    dyn_sc = dyn_sscha.GenerateSupercellDyn(sc)

    # ── SSCHA minimization ───────────────────────────────────────────────────
    calculator = QEForceCalculator(cfg, sscha_dir / "force_calcs", qe_bin, pseudo_dir)

    print(f"\nStarting SSCHA minimization:")
    print(f"  Temperature: {temperature} K")
    print(f"  Configs per iteration: {nconfigs}")
    print(f"  Max iterations: {niter}")
    print(f"  Supercell: {sc}")

    for iteration in range(start_from_iter, niter):
        iter_label = f"iter_{iteration:03d}"
        checkpoint = sscha_dir / f"{iter_label}_dyn"

        print(f"\n{'='*60}")
        print(f"  SSCHA iteration {iteration+1}/{niter}")
        print(f"{'='*60}")

        # Check for checkpoint
        if (sscha_dir / f"{iter_label}_dyn1").exists():
            print(f"  Loading checkpoint: {iter_label}")
            dyn_sc = Phonons.Phonons(str(sscha_dir / f"{iter_label}_dyn"), nqirr=1)

        # Generate random displacements from current trial force constants
        print(f"  Generating {nconfigs} random configurations at T={temperature} K...")
        ensemble = Ensemble.Ensemble(dyn_sc, T0=temperature, supercell=sc)
        ensemble.generate(nconfigs, sobol=True)

        # Compute forces with QE
        print(f"  Computing forces ({nconfigs} pw.x runs)...")
        calculator.compute_ensemble(ensemble, iteration)

        # SSCHA gradient step
        print(f"  Updating trial force constants...")
        import sscha.SchaMinimizer as SchaMin
        minim = SchaMin.SSCHA_Minimizer(
            ensemble,
            kong_liu_ratio=0.5,
            meaningful_factor=0.2,
            minimization_algorithm="sdes",
        )
        minim.init()
        minim.run(verbose=False)
        minim.finalize()

        # Save checkpoint
        dyn_sc = minim.dyn
        dyn_sc.save_qe(str(sscha_dir / f"{iter_label}_dyn"))
        print(f"  Saved checkpoint: {iter_label}_dyn")

        # Check convergence
        gradient_norm = np.max(np.abs(minim.gradient))
        print(f"  Gradient norm: {gradient_norm:.6e}")

        # Check imaginary modes in current dyn
        freqs_iter = []
        for iq in range(dyn_sc.nqirr):
            w2, _ = np.linalg.eigh(dyn_sc.dynmats[iq])
            freqs_iter.extend(np.sign(w2) * np.sqrt(np.abs(w2)))
        n_imag_iter = sum(1 for f in freqs_iter if f < -1e-3)
        print(f"  Imaginary modes remaining: {n_imag_iter}")

        if gradient_norm < 1e-4:
            print(f"\n  ✓ Converged at iteration {iteration+1}")
            break

    # ── Save SSCHA supercell dyn and interpolate back to original q-grid ─────
    print(f"\nSaving SSCHA force constants and interpolating to original q-grid...")
    # Save the supercell (full q-grid) dyn — this IS the renormalized result
    out_prefix = str(sscha_dir / f"{prefix}_sscha.dyn")
    dyn_sc.save_qe(out_prefix)

    # Also interpolate back to the original irreducible q-point representation
    # so q2r/matdyn/EPW can use them directly (same format as harmonic dyn files)
    dyn_renorm = dyn_sc.Interpolate(np.array(dyn.q_tot), False)
    dyn_renorm.save_qe(str(sscha_dir / f"{prefix}_sscha_irred.dyn"))
    print(f"  Saved supercell dyn:   {out_prefix}1...")
    print(f"  Saved irred q-pt dyn:  {sscha_dir}/{prefix}_sscha_irred.dyn1...")

    # Check final imaginary modes
    freqs_final = []
    for iq in range(dyn_renorm.nqirr):
        w2, _ = np.linalg.eigh(dyn_renorm.dynmats[iq])
        freqs_final.extend(np.sign(w2) * np.sqrt(np.abs(w2)))
    n_imag_final = sum(1 for f in freqs_final if f < -1e-3)

    print(f"\n{'='*60}")
    print(f"  SSCHA complete")
    print(f"  Harmonic imaginary modes: {n_imag}")
    print(f"  SSCHA imaginary modes:    {n_imag_final}")
    print(f"  Output: {out_prefix}*")
    print(f"{'='*60}")
    print(f"""
Next steps:
  1. Copy SSCHA dyn files over harmonic ones:
       cp {sscha_dir}/{prefix}_sscha.dyn* {phonon_dir}/
       (rename to {prefix}.dyn1 ... {prefix}.dyn{nqirr})
  2. Re-run q2r + matdyn + EPW with SSCHA force constants:
       python pipeline.py {yaml_path} --start-from phonons
     (phonon step will be cached since dyn files exist; q2r will rerun)
""")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SSCHA anharmonic phonon renormalization for a compound YAML."
    )
    parser.add_argument("yaml", help="Path to compound YAML (e.g. compounds/CaH6_200GPa.yaml)")
    parser.add_argument(
        "--supercell", nargs=3, type=int, default=[2, 2, 2], metavar=("N1", "N2", "N3"),
        help="Supercell dimensions (default: 2 2 2)"
    )
    parser.add_argument(
        "--nconfigs", type=int, default=50,
        help="Number of random configurations per SSCHA iteration (default: 50)"
    )
    parser.add_argument(
        "--niter", type=int, default=20,
        help="Maximum number of SSCHA iterations (default: 20)"
    )
    parser.add_argument(
        "--temp", type=float, default=300.0,
        help="Temperature in Kelvin for SSCHA (default: 300 K)"
    )
    parser.add_argument(
        "--start-from-iter", type=int, default=0, metavar="N",
        help="Resume from iteration N (default: 0)"
    )
    args = parser.parse_args()

    run_sscha(
        args.yaml,
        supercell=args.supercell,
        nconfigs=args.nconfigs,
        niter=args.niter,
        temperature=args.temp,
        start_from_iter=args.start_from_iter,
    )
