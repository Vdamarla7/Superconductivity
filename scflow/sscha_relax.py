"""
SSCHA anharmonic phonon renormalization for the QE + EPW pipeline.

Usage:
    python -m scflow.sscha_relax compounds/PdH_0GPa.yaml [--supercell 2 2 2]
                                                      [--nconfigs 50]
                                                      [--niter 20]
                                                      [--temp 300]
                                                      [--start-from-iter N]

What it does:
    1. Loads harmonic dynamical matrices from <run-dir>/2_phonons/<prefix>.dyn*
    2. Makes them positive-definite (replaces imaginary modes with small positive freqs)
       so SSCHA can use them as a starting trial force constants
    3. Runs SSCHA minimization: generates random structures, computes QE forces,
       updates trial force constants until convergence (dF < threshold)
    4. Writes renormalized dyn files to <run-dir>/sscha_<supercell>/
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

from scflow.config import default_outdir, load_config


# ─────────────────────────────────────────────────────────────────────────────
# QE force calculator
# ─────────────────────────────────────────────────────────────────────────────

def _write_scf_input(structure, cfg, scf_path, pseudo_dir, pressure_kbar, sc=None):
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
/

&ELECTRONS
  conv_thr    = 1.0d-7
  mixing_beta = 0.4
  diagonalization = 'david'
/

ATOMIC_SPECIES
{chr(10).join(species_lines)}

CELL_PARAMETERS {{angstrom}}
{chr(10).join(cell_lines)}

ATOMIC_POSITIONS {{angstrom}}
{chr(10).join(pos_lines)}

K_POINTS {{automatic}}
{max(1, -(-d['scf_kgrid'][0] // (sc[0] if sc else 1)))} {max(1, -(-d['scf_kgrid'][1] // (sc[1] if sc else 1)))} {max(1, -(-d['scf_kgrid'][2] // (sc[2] if sc else 1)))} 0 0 0
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

    def __init__(self, cfg, run_dir, qe_bin, pseudo_dir, nproc=1, sc=None, gpu_pw_bin=None, gpu_mpirun=None):
        self.cfg = cfg
        self.run_dir = Path(run_dir)
        self.qe_bin = qe_bin
        self.pseudo_dir = pseudo_dir
        self.pressure_kbar = cfg["compound"]["pressure_GPa"] * 10.0
        self.nproc = nproc
        self.sc = sc
        self.gpu_pw_bin = gpu_pw_bin    # if set, use GPU pw.x instead of qe_bin/pw.x
        self.gpu_mpirun = gpu_mpirun    # HPC SDK mpirun required for GPU build

    def _build_pw_cmd(self):
        """Return (pw_cmd_prefix, env) for a single pw.x call."""
        if self.gpu_pw_bin:
            mpirun = self.gpu_mpirun or "mpirun"
            nvhpc_lib  = str(Path(mpirun).parent.parent / "lib")
            nvhpc_comp = str(Path(mpirun).parents[3] / "compilers" / "lib")
            nvhpc_cuda = str(Path(mpirun).parents[3] / "cuda" / "lib64")
            env = {**os.environ,
                   "OMP_NUM_THREADS": "4",
                   "CUDA_VISIBLE_DEVICES": "0"}
            env["LD_LIBRARY_PATH"] = ":".join(
                [nvhpc_lib, nvhpc_comp, nvhpc_cuda, env.get("LD_LIBRARY_PATH", "")])
            return [mpirun, "-n", "1", self.gpu_pw_bin, "-in", "scf.in"], env
        else:
            pw_bin = str(Path(self.qe_bin) / "pw.x")
            env = {**os.environ, "OMP_NUM_THREADS": "1"}
            if self.nproc > 1:
                return ["mpirun", "-np", str(self.nproc), pw_bin, "-in", "scf.in"], env
            return [pw_bin, "-in", "scf.in"], env

    def _run_one(self, args):
        """Run pw.x for a single configuration. Returns (i, energy, forces) or raises."""
        i, structure, work_dir, ensemble_N = args
        conf_dir = work_dir / f"conf_{i:04d}"
        conf_dir.mkdir(exist_ok=True)
        (conf_dir / "tmp").mkdir(exist_ok=True)
        scf_in  = conf_dir / "scf.in"
        scf_out = conf_dir / "scf.out"
        _write_scf_input(structure, self.cfg, scf_in, self.pseudo_dir, self.pressure_kbar, sc=self.sc)
        if scf_out.exists():
            try:
                e, f = _parse_qe_output(scf_out, structure)
                return i, e, f, True   # cached
            except Exception:
                pass
        pw_cmd, env = self._build_pw_cmd()
        result = subprocess.run(pw_cmd, stdout=open(scf_out, "w"),
                                stderr=subprocess.STDOUT, cwd=str(conf_dir), env=env)
        if result.returncode != 0:
            return i, None, None, False   # signal failure without crashing
        e, f = _parse_qe_output(scf_out, structure)
        return i, e, f, False

    def compute_ensemble(self, ensemble, iter_idx):
        """Run pw.x on all structures in ensemble in parallel, fill energies and forces."""
        import concurrent.futures, math
        work_dir = self.run_dir / f"iter_{iter_idx:03d}"
        work_dir.mkdir(parents=True, exist_ok=True)

        energies = np.zeros(ensemble.N)
        forces = np.zeros((ensemble.N, ensemble.structures[0].N_atoms, 3))
        failed_indices = []

        # Parallel workers: GPU fits ~6 concurrent (3.4 GB each, 24 GB total);
        # CPU uses nproc MPI tasks each, so 1 concurrent job.
        if self.gpu_pw_bin:
            max_workers = 6
            label = "GPU×6"
        else:
            max_workers = 1
            label = f"np={self.nproc}"

        args = [(i, s, work_dir, ensemble.N) for i, s in enumerate(ensemble.structures)]
        completed = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {pool.submit(self._run_one, a): a[0] for a in args}
            for fut in concurrent.futures.as_completed(futures):
                i, e, f, cached = fut.result()
                if e is None:
                    # pw.x failed for this config — skip it, don't fill energies/forces
                    failed_indices.append(i)
                    completed += 1
                    print(f"  conf {i:4d}/{ensemble.N} ({label}): FAILED (skipped)  [{completed}/{ensemble.N} done]",
                          flush=True)
                    continue
                energies[i] = e
                forces[i]   = f
                completed += 1
                tag = "cached" if cached else f"E={e:.6f} Ry"
                print(f"  conf {i:4d}/{ensemble.N} ({label}): {tag}  [{completed}/{ensemble.N} done]",
                      flush=True)

        if failed_indices:
            n_ok = ensemble.N - len(failed_indices)
            if n_ok < max(ensemble.N * 0.75, 8):
                raise RuntimeError(
                    f"Too many failed configs: {len(failed_indices)}/{ensemble.N}. "
                    f"Check {work_dir} for pw.x errors."
                )
            print(f"  WARNING: {len(failed_indices)} configs failed (conf {failed_indices}), "
                  f"proceeding with {n_ok}/{ensemble.N} for gradient.")
            good = np.array([i for i in range(ensemble.N) if i not in set(failed_indices)])
            # Filter every per-config array so shapes stay consistent
            ensemble.structures  = [ensemble.structures[i]  for i in good]
            ensemble.energies    = energies[good] * Units.RY_TO_EV
            ensemble.forces      = forces[good]   * (Units.RY_TO_EV / Units.BOHR_TO_ANGSTROM)
            if hasattr(ensemble, 'xats') and isinstance(ensemble.xats, np.ndarray):
                ensemble.xats = ensemble.xats[good]
            if hasattr(ensemble, 'u_disps') and isinstance(ensemble.u_disps, np.ndarray):
                ensemble.u_disps = ensemble.u_disps[good]
            if hasattr(ensemble, 'u_disps_original') and isinstance(ensemble.u_disps_original, np.ndarray):
                ensemble.u_disps_original = ensemble.u_disps_original[good]
            if hasattr(ensemble, 'rho') and isinstance(ensemble.rho, np.ndarray):
                ensemble.rho = ensemble.rho[good]
            if hasattr(ensemble, 'sscha_energies') and isinstance(ensemble.sscha_energies, np.ndarray):
                ensemble.sscha_energies = ensemble.sscha_energies[good]
            if hasattr(ensemble, 'sscha_forces') and isinstance(ensemble.sscha_forces, np.ndarray):
                ensemble.sscha_forces = ensemble.sscha_forces[good]
            # force_computed must match ensemble.N; all kept configs are fully computed
            if hasattr(ensemble, 'force_computed') and isinstance(ensemble.force_computed, np.ndarray):
                ensemble.force_computed = np.ones(len(good), dtype=bool)
            ensemble.N           = len(good)
            ensemble.has_stress  = False
            return

        # Fill ensemble
        ensemble.energies = energies * Units.RY_TO_EV
        ensemble.forces = forces * (Units.RY_TO_EV / Units.BOHR_TO_ANGSTROM)
        ensemble.has_stress = False


# ─────────────────────────────────────────────────────────────────────────────
# MACE force calculator
# ─────────────────────────────────────────────────────────────────────────────

class MACEForceCalculator:
    """Computes forces on an ensemble of structures using a trained MACE model."""

    def __init__(self, model_path, device='cuda'):
        self.model_path = str(model_path)
        self.device = device
        self._calc = None

    def _calc_instance(self):
        if self._calc is None:
            from mace.calculators import MACECalculator
            self._calc = MACECalculator(
                model_paths=self.model_path,
                device=self.device,
                default_dtype='float32',
            )
            print(f"  MACE model loaded: {self.model_path} ({self.device})")
        return self._calc

    @staticmethod
    def _to_ase(structure):
        import ase
        return ase.Atoms(
            symbols=list(structure.atoms),
            positions=structure.coords.copy(),   # Å
            cell=structure.unit_cell.copy(),      # Å
            pbc=True,
        )

    def compute_ensemble(self, ensemble, iter_idx):
        """Run MACE on all structures, fill ensemble energies and forces (eV, eV/Å)."""
        calc = self._calc_instance()
        N = ensemble.N
        nat = ensemble.structures[0].N_atoms
        energies = np.zeros(N)
        forces   = np.zeros((N, nat, 3))

        for i, structure in enumerate(ensemble.structures):
            atoms = self._to_ase(structure)
            atoms.calc = calc
            energies[i] = atoms.get_potential_energy()   # eV
            forces[i]   = atoms.get_forces()             # eV/Å
            if (i + 1) % 6 == 0 or i == N - 1:
                print(f"  conf {i+1:4d}/{N} (MACE): E={energies[i]:.4f} eV", flush=True)

        ensemble.energies   = energies
        ensemble.forces     = forces
        ensemble.has_stress = False


# ─────────────────────────────────────────────────────────────────────────────
# Main SSCHA loop
# ─────────────────────────────────────────────────────────────────────────────

def run_sscha(yaml_path, supercell, nconfigs, niter, temperature, start_from_iter, nproc=1,
              dyn_prefix_override=None, nqirr_override=None, run_dir=None):
    cfg = load_config(yaml_path)
    prefix = cfg["compound"]["prefix"]
    qe_bin = cfg["paths"]["qe_bin"]
    pseudo_dir = cfg["paths"]["pseudo_dir"]
    gpu_pw_bin = os.environ.get("SCFLOW_GPU_PW_BIN") or cfg["paths"].get("gpu_pw_bin")
    gpu_mpirun = os.environ.get("SCFLOW_GPU_MPIRUN") or cfg["paths"].get("gpu_mpirun")

    run_dir = Path(run_dir) if run_dir else default_outdir(cfg)
    phonon_dir = run_dir / "2_phonons"
    sc_tag = "x".join(map(str, supercell)) if supercell else "default"
    sscha_dir  = run_dir / f"sscha_{sc_tag}"
    sscha_dir.mkdir(parents=True, exist_ok=True)

    # ── Load harmonic dynamical matrices ──────────────────────────────────────
    if dyn_prefix_override:
        dyn_prefix = dyn_prefix_override
        dyn_glob_prefix = Path(dyn_prefix_override)
        nqirr = nqirr_override if nqirr_override else (
            sum(1 for p in dyn_glob_prefix.parent.glob(f"{dyn_glob_prefix.name}[0-9]*")
                if p.stat().st_size > 100) - 1)
    else:
        dyn_prefix = str(phonon_dir / f"{prefix}.dyn")
        nqirr = sum(1 for p in phonon_dir.glob(f"{prefix}.dyn[0-9]*") if p.stat().st_size > 100) - 1

    # Read q-grid from dyn0 (first line: nq1 nq2 nq3)
    dyn0_path = Path(dyn_prefix + "0")
    dyn0_tokens = dyn0_path.read_text().split()
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
    dyn_sscha.Symmetrize()  # re-enforce acoustic sum rule broken by ForcePositiveDefinite
    print("  Applied ForcePositiveDefinite + Symmetrize (ASR restored)")

    # Use user-specified supercell; if smaller than q-grid, interpolate dyn down
    sc = supercell if supercell else qgrid
    if sc != qgrid:
        print(f"  Interpolating dyn from {qgrid} → {sc} supercell...")
        interp_prefix = str(sscha_dir / f"{prefix}_sc{'x'.join(map(str,sc))}.dyn")
        dyn_sscha.Interpolate(qgrid, sc)   # (coarse_grid, fine_grid)
        dyn_sscha.save_qe(interp_prefix)
        n_interp = sum(1 for p in Path(sscha_dir).glob(
            f"{prefix}_sc{'x'.join(map(str,sc))}.dyn[0-9]*") if p.stat().st_size > 100)
        dyn_sscha = Phonons.Phonons(interp_prefix, nqirr=n_interp)
        print(f"  Interpolated: {dyn_sscha.nqirr} irred q-pts on {sc} grid")

    # ── SSCHA minimization ───────────────────────────────────────────────────
    # dyn_sscha is the primitive-cell dyn (nqirr q-points). Ensemble takes this
    # plus the supercell size — it builds the supercell internally. minim.dyn
    # returns the updated primitive-cell dyn, so we track dyn_sscha throughout.
    sscha_cfg = cfg.get("sscha", {})
    fc_type   = sscha_cfg.get("force_calculator", "qe")

    if fc_type == "mace":
        model_path = sscha_cfg.get("mace_model")
        if not model_path or not Path(model_path).exists():
            raise FileNotFoundError(
                f"MACE model not found: {model_path!r}\n"
                f"Train one first: bash pipeline/train_mace.sh finetune"
            )
        device = sscha_cfg.get("mace_device", "cuda")
        print(f"  MACE mode: {model_path} ({device})")
        calculator = MACEForceCalculator(model_path, device=device)
    else:
        if gpu_pw_bin:
            print(f"  GPU mode: using {gpu_pw_bin}")
        calculator = QEForceCalculator(cfg, sscha_dir / "force_calcs", qe_bin, pseudo_dir,
                                       nproc=nproc, sc=sc, gpu_pw_bin=gpu_pw_bin, gpu_mpirun=gpu_mpirun)

    print(f"\nStarting SSCHA minimization:")
    print(f"  Temperature: {temperature} K")
    print(f"  Configs per iteration: {nconfigs}")
    print(f"  Max iterations: {niter}")
    print(f"  Supercell: {sc}")

    for iteration in range(start_from_iter, niter):
        iter_label = f"iter_{iteration:03d}"

        print(f"\n{'='*60}")
        print(f"  SSCHA iteration {iteration+1}/{niter}")
        print(f"{'='*60}")

        # If this iteration already completed (checkpoint exists), load dyn and skip
        if (sscha_dir / f"{iter_label}_dyn1").exists():
            print(f"  Skipping completed iteration (checkpoint {iter_label}_dyn exists)")
            dyn_sscha = Phonons.Phonons(str(sscha_dir / f"{iter_label}_dyn"), nqirr=nqirr)
            dyn_sscha.Symmetrize()
            continue

        # Ensure dyn is positive definite before sampling.
        # Order matters: Symmetrize() (ASR) first, then ForcePositiveDefinite() last —
        # doing it in reverse lets Symmetrize() re-introduce negative eigenvalues.
        dyn_sscha.Symmetrize()
        dyn_sscha.ForcePositiveDefinite()

        # Generate random displacements — pass primitive-cell dyn + supercell size
        print(f"  Generating {nconfigs} random configurations at T={temperature} K...")
        ensemble = Ensemble.Ensemble(dyn_sscha, T0=temperature, supercell=sc)
        ensemble.generate(nconfigs, sobol=True)

        # Compute forces with QE
        print(f"  Computing forces ({nconfigs} pw.x runs)...")
        calculator.compute_ensemble(ensemble, iteration)

        # SSCHA gradient step
        print(f"  Updating trial force constants...")
        import sscha.SchaMinimizer as SchaMin
        minim = SchaMin.SSCHA_Minimizer(
            ensemble,
            kong_liu_ratio=0.0,
            meaningful_factor=0.2,
            minimization_algorithm="sdes",
        )
        minim.max_ka = 2              # at most 2 gradient steps per outer SSCHA iteration
        minim.max_imag_freq_counter = 200  # allow more line-search attempts (default 50 too low for early iters with many imaginary modes)
        minim.init()
        minim.run(verbose=False)
        minim.finalize()

        # minim.dyn is the updated primitive-cell dyn; re-enforce ASR before saving
        dyn_sscha = minim.dyn
        dyn_sscha.Symmetrize()
        dyn_sscha.save_qe(str(sscha_dir / f"{iter_label}_dyn"))
        print(f"  Saved checkpoint: {iter_label}_dyn")

        # Check convergence — __gc__ stores ||dyn_grad|| per minimization step
        gradient_norm = minim.__gc__[-1] if minim.__gc__ else float("inf")
        print(f"  Gradient norm: {gradient_norm:.6e}")

        # Count imaginary modes in updated dyn
        freqs_iter = []
        for iq in range(dyn_sscha.nqirr):
            w2, _ = np.linalg.eigh(dyn_sscha.dynmats[iq])
            freqs_iter.extend(np.sign(w2) * np.sqrt(np.abs(w2)))
        n_imag_iter = sum(1 for f in freqs_iter if f < -1e-3)
        print(f"  Imaginary modes remaining: {n_imag_iter}")

        if gradient_norm < 1e-4:
            print(f"\n  ✓ Converged at iteration {iteration+1}")
            break

    # ── Save renormalized dyn files ───────────────────────────────────────────
    print(f"\nSaving SSCHA renormalized dynamical matrices...")
    # dyn_sscha is now the primitive-cell dyn on the original irreducible q-grid —
    # same format as the harmonic dyn files, ready for q2r → matdyn → EPW.
    out_prefix = str(sscha_dir / f"{prefix}_sscha_irred.dyn")
    dyn_sscha.save_qe(out_prefix)
    print(f"  Saved irred q-pt dyn:  {out_prefix}1...")

    # Check final imaginary modes
    freqs_final = []
    for iq in range(dyn_sscha.nqirr):
        w2, _ = np.linalg.eigh(dyn_sscha.dynmats[iq])
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
       cp {sscha_dir}/{prefix}_sscha_irred.dyn* {phonon_dir}/
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
    parser.add_argument(
        "--nproc", type=int, default=1,
        help="MPI ranks per pw.x call — keep at 1 on Mac (default: 1)"
    )
    parser.add_argument(
        "--dyn-prefix", type=str, default=None,
        help="Override path to dyn files (e.g. path/to/CaH6_200_333.dyn)"
    )
    parser.add_argument(
        "--nqirr", type=int, default=None,
        help="Override number of irreducible q-points (auto-detected if not set)"
    )
    parser.add_argument("--run-dir", default=None, help="scflow run directory (default: runs/<prefix>)")
    args = parser.parse_args()

    run_sscha(
        args.yaml,
        supercell=args.supercell,
        nconfigs=args.nconfigs,
        niter=args.niter,
        temperature=args.temp,
        start_from_iter=args.start_from_iter,
        nproc=args.nproc,
        dyn_prefix_override=args.dyn_prefix,
        nqirr_override=args.nqirr,
        run_dir=args.run_dir,
    )
