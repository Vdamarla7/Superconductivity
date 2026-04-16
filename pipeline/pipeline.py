#!/usr/bin/env python3
"""
EPW Superconductor Pipeline
===========================
Full pipeline from crystal structure to Eliashberg Tc.

Usage:
    python pipeline.py compounds/H3S_200GPa.yaml
    python pipeline.py compounds/H3S_200GPa.yaml --outdir /path/to/output
    python pipeline.py compounds/H3S_200GPa.yaml --start-from phonons
    python pipeline.py compounds/H3S_200GPa.yaml --dry-run

Steps (can resume from any step with --start-from):
  vcrelax   → vc-relax at target pressure (skip if A already in config)
  scf       → ground-state SCF
  phonons   → ph.x harmonic phonons + q2r + matdyn
  nscf      → NSCF on coarse k-grid for EPW
  epw1      → EPW: Wannierize, write .epmatwp
  epw2      → EPW: fine-grid interpolation + Eliashberg Tc
"""

import sys, os, yaml, subprocess, shutil, argparse
from pathlib import Path

# Add pipeline dir to path so imports work regardless of cwd
sys.path.insert(0, str(Path(__file__).parent))
from generators import gen_vcrelax, gen_scf, gen_phonon, gen_q2r, gen_matdyn, gen_nscf, gen_epw_pass1, gen_epw_pass2
from postprocess import job_done, parse_vcrelax_a, parse_fermi_energy, parse_sym_ops, parse_imaginary_modes, parse_allen_dynes_tc, parse_nqirr

STEPS = ["vcrelax", "scf", "phonons", "nscf", "epw1", "epw2"]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def banner(msg):
    print(f"\n{'─'*60}\n  {msg}\n{'─'*60}")

def status(step, n, total, msg):
    print(f"  [{n}/{total}] {step:<12} {msg}", flush=True)

def run_qe(binary, infile, outfile, qe_bin, background=False, cwd=None):
    """Run a QE binary. cwd sets working dir so QE resolves outdir/fildyn correctly."""
    cmd = [str(Path(qe_bin) / binary), "-in", str(Path(infile).resolve())]
    if background:
        proc = subprocess.Popen(cmd, stdout=open(outfile, "w"), stderr=subprocess.STDOUT, cwd=cwd)
        return True, proc.pid
    with open(outfile, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=cwd)
    if result.returncode != 0 or not job_done(outfile):
        print(f"\n  ✗ {binary} FAILED — check {outfile}")
        sys.exit(1)
    return True, None

def write_input(path, content):
    Path(path).write_text(content)


# ─────────────────────────────────────────────────────────────────────────────
# Individual step runners
# ─────────────────────────────────────────────────────────────────────────────

def run_vcrelax(cfg, dirs, qe_bin, dry_run):
    outfile = dirs["vcrelax"] / "vcrelax.out"
    if job_done(outfile):
        A = parse_vcrelax_a(outfile, cfg["structure"]["ibrav"])
        print(f"  [1/{len(STEPS)}] vcrelax      ✓ (cached) A = {A} Å")
        _update_A(cfg, A)
        return

    write_input(dirs["vcrelax"] / "vcrelax.in", gen_vcrelax(cfg))
    if dry_run:
        print(f"  [1/{len(STEPS)}] vcrelax      (dry-run) input written")
        return

    print(f"  [1/{len(STEPS)}] vcrelax      running...", end=" ", flush=True)
    run_qe("pw.x", dirs["vcrelax"] / "vcrelax.in", outfile, qe_bin, cwd=str(dirs["vcrelax"]))
    A = parse_vcrelax_a(outfile, cfg["structure"]["ibrav"])
    _update_A(cfg, A)
    print(f"done  A = {A} Å")


def _update_A(cfg, A):
    """Update cfg with relaxed geometry from vc-relax."""
    if cfg["structure"]["ibrav"] == 3 and A is not None:
        cfg["structure"]["A"] = A
    # For ibrav=0 the cell params are updated differently — future work


def run_scf(cfg, dirs, qe_bin, dry_run, step_n):
    outfile = dirs["phonons"] / "scf.out"
    if job_done(outfile):
        Ef = parse_fermi_energy(outfile)
        nsym = parse_sym_ops(outfile)
        print(f"  [{step_n}/{len(STEPS)}] scf          ✓ (cached) Ef={Ef:.2f} eV, {nsym} sym ops")
        return

    write_input(dirs["phonons"] / "scf.in", gen_scf(cfg))
    if dry_run:
        print(f"  [{step_n}/{len(STEPS)}] scf          (dry-run) input written")
        return

    print(f"  [{step_n}/{len(STEPS)}] scf          running...", end=" ", flush=True)
    (dirs["phonons"] / "tmp").mkdir(exist_ok=True)
    run_qe("pw.x", dirs["phonons"] / "scf.in", outfile, qe_bin, cwd=str(dirs["phonons"]))
    Ef = parse_fermi_energy(outfile)
    nsym = parse_sym_ops(outfile)
    print(f"done  Ef={Ef:.2f} eV, {nsym} sym ops")


def run_phonons(cfg, dirs, qe_bin, dry_run, step_n):
    ph_out   = dirs["phonons"] / "ph.out"
    q2r_out  = dirs["phonons"] / "q2r.out"
    mat_out  = dirs["phonons"] / "matdyn.out"

    prefix = cfg["compound"]["prefix"]

    if not job_done(ph_out):
        write_input(dirs["phonons"] / "ph.in", gen_phonon(cfg))
        write_input(dirs["phonons"] / "q2r.in", gen_q2r(cfg))
        write_input(dirs["phonons"] / "matdyn.in", gen_matdyn(cfg))
        if dry_run:
            print(f"  [{step_n}/{len(STEPS)}] phonons      (dry-run) inputs written")
            return
        print(f"  [{step_n}/{len(STEPS)}] phonons      running ph.x (this takes a while)...", end=" ", flush=True)
        run_qe("ph.x", dirs["phonons"] / "ph.in", ph_out, qe_bin, cwd=str(dirs["phonons"]))
        nqirr = parse_nqirr(ph_out)
        imag = parse_imaginary_modes(ph_out)
        flag = "✓ stable" if not imag else f"⚠ {len(imag)} imaginary modes"
        print(f"done  {nqirr} irred q-pts, {flag}")
    else:
        nqirr = parse_nqirr(ph_out)
        imag  = parse_imaginary_modes(ph_out)
        flag  = "✓ stable" if not imag else f"⚠ {len(imag)} imaginary modes"
        print(f"  [{step_n}/{len(STEPS)}] phonons      ✓ (cached) {nqirr} irred q-pts, {flag}")

    # pp.py: must run from tmp/ with dyn symlinks
    save_dir = dirs["phonons"] / "tmp" / "save"
    if not save_dir.exists():
        _run_pp_py(cfg, dirs, prefix, qe_bin)
        # Symlink for EPW
        save_link = dirs["phonons"] / "save"
        if not save_link.exists():
            save_link.symlink_to(save_dir.resolve())

    # q2r
    if not job_done(q2r_out):
        print(f"  [{step_n}/{len(STEPS)}] q2r          running...", end=" ", flush=True)
        run_qe("q2r.x", dirs["phonons"] / "q2r.in", q2r_out, qe_bin, cwd=str(dirs["phonons"]))
        print("done")

    # matdyn
    if not job_done(mat_out):
        print(f"  [{step_n}/{len(STEPS)}] matdyn       running...", end=" ", flush=True)
        run_qe("matdyn.x", dirs["phonons"] / "matdyn.in", mat_out, qe_bin, cwd=str(dirs["phonons"]))
        print("done")


def _run_pp_py(cfg, dirs, prefix, qe_bin):
    """
    Run pp.py from the phonons/tmp directory to create the EPW save/.

    pp.py expects:
      _ph0/{prefix}.phsave/control_ph.xml
      _ph0/{prefix}.q_N/{prefix}.dvscf1  (for N>=2)
      _ph0/{prefix}.dvscf1               (for q=1)
    These are created by ph.x when cwd=phonons/ (outdir=./tmp/).
    We build _ph0/ symlinks pointing at the actual ph.x output.
    """
    tmp_dir = dirs["phonons"] / "tmp"
    n_dyn = len(list(dirs["phonons"].glob(f"{prefix}.dyn[0-9]*")))

    # Symlink dyn files into tmp/ (pp.py reads them from cwd)
    for i in range(n_dyn + 1):
        src = (dirs["phonons"] / f"{prefix}.dyn{i}").resolve()
        dst = tmp_dir / f"{prefix}.dyn{i}"
        if src.exists() and not dst.exists():
            dst.symlink_to(src)

    # Build _ph0/ symlinks so pp.py can find phsave, dvscf, q_N dirs
    ph0_dir = tmp_dir / "_ph0"
    ph0_dir.mkdir(exist_ok=True)
    for src in tmp_dir.iterdir():
        if src.name.startswith(prefix) and src.name != f"{prefix}.save":
            dst = ph0_dir / src.name
            if not dst.exists():
                dst.symlink_to(src.resolve())

    pp_py = Path(qe_bin).parent / "EPW" / "bin" / "pp.py"
    proc = subprocess.run(
        ["python3", str(pp_py)],
        input=f"{prefix}\n", text=True,
        capture_output=True, cwd=str(tmp_dir)
    )
    if not (tmp_dir / "save").exists():
        print(f"\n  ✗ pp.py failed:\n{proc.stderr}")
        sys.exit(1)
    print(f"  pp.py done → save/ created ({n_dyn} dvscf files)")


def run_nscf(cfg, dirs, qe_bin, dry_run, step_n):
    nscf_out = dirs["epw"] / "nscf.out"

    # Copy SCF save data to EPW tmp if not already there
    scf_save = dirs["phonons"] / "tmp" / f"{cfg['compound']['prefix']}.save"
    epw_save = dirs["epw"] / "tmp" / f"{cfg['compound']['prefix']}.save"
    if not epw_save.exists() and scf_save.exists():
        (dirs["epw"] / "tmp").mkdir(exist_ok=True)
        shutil.copytree(str(scf_save), str(epw_save))

    if job_done(nscf_out):
        Ef = parse_fermi_energy(nscf_out)
        print(f"  [{step_n}/{len(STEPS)}] nscf         ✓ (cached) Ef={Ef:.2f} eV")
        return

    write_input(dirs["epw"] / "nscf.in", gen_nscf(cfg))
    if dry_run:
        print(f"  [{step_n}/{len(STEPS)}] nscf         (dry-run) input written")
        return

    print(f"  [{step_n}/{len(STEPS)}] nscf         running...", end=" ", flush=True)
    run_qe("pw.x", dirs["epw"] / "nscf.in", nscf_out, qe_bin, cwd=str(dirs["epw"]))
    Ef = parse_fermi_energy(nscf_out)
    print(f"done  Ef={Ef:.2f} eV")


def run_epw1(cfg, dirs, qe_bin, dry_run, step_n):
    epmatwp = dirs["epw"] / "tmp" / f"{cfg['compound']['prefix']}.epmatwp"
    out1    = dirs["epw"] / "epw_pass1.out"

    if epmatwp.exists():
        size_mb = epmatwp.stat().st_size / 1e6
        print(f"  [{step_n}/{len(STEPS)}] epw (pass1)  ✓ (cached) .epmatwp = {size_mb:.0f} MB")
        return

    write_input(dirs["epw"] / "epw_pass1.in", gen_epw_pass1(cfg))
    if dry_run:
        print(f"  [{step_n}/{len(STEPS)}] epw (pass1)  (dry-run) input written")
        return

    print(f"  [{step_n}/{len(STEPS)}] epw (pass1)  Wannierizing + writing .epmatwp...", end=" ", flush=True)
    run_qe("epw.x", dirs["epw"] / "epw_pass1.in", out1, qe_bin, cwd=str(dirs["epw"]))
    size_mb = epmatwp.stat().st_size / 1e6 if epmatwp.exists() else 0
    print(f"done  .epmatwp = {size_mb:.0f} MB")


def run_epw2(cfg, dirs, qe_bin, dry_run, step_n):
    out2 = dirs["epw"] / "epw_pass2.out"

    # EPW's write_freq() writes into {prefix}.ephmat/ before kqmap_fine creates it.
    # Pre-create so the file open doesn't fail.
    ephmat_dir = dirs["epw"] / "tmp" / f"{cfg['compound']['prefix']}.ephmat"
    ephmat_dir.mkdir(exist_ok=True)

    write_input(dirs["epw"] / "epw_pass2.in", gen_epw_pass2(cfg))

    if job_done(out2):
        tc, wlog = parse_allen_dynes_tc(out2)
        print(f"  [{step_n}/{len(STEPS)}] epw (pass2)  ✓ (cached)")
        _print_results(tc, wlog, cfg)
        return

    if dry_run:
        print(f"  [{step_n}/{len(STEPS)}] epw (pass2)  (dry-run) input written")
        return

    nkf = cfg["epw"]["fine_kgrid"]
    print(f"  [{step_n}/{len(STEPS)}] epw (pass2)  fine grid {nkf[0]}³ + Eliashberg...", end=" ", flush=True)
    run_qe("epw.x", dirs["epw"] / "epw_pass2.in", out2, qe_bin, cwd=str(dirs["epw"]))
    tc, wlog = parse_allen_dynes_tc(out2)
    print("done")
    _print_results(tc, wlog, cfg)


def _print_results(tc, wlog, cfg):
    muc = cfg["eliashberg"]["muc"]
    print(f"\n  {'═'*40}")
    print(f"  RESULTS")
    print(f"  {'═'*40}")
    if wlog:
        print(f"  ω_log         = {wlog:.1f} meV ({wlog*11.6045:.0f} K)")
    if tc:
        print(f"  Allen-Dynes Tc = {tc:.1f} K  (μ* = {muc})")
    print(f"  {'═'*40}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="EPW superconductor pipeline: structure → Eliashberg Tc",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  python pipeline.py compounds/H3S_200GPa.yaml\n"
               "  python pipeline.py compounds/MyCompound.yaml --start-from epw1"
    )
    parser.add_argument("config", help="Path to compound YAML (see compounds/template.yaml)")
    parser.add_argument("--outdir", default=None,
                        help="Output directory (default: runs/<prefix>)")
    parser.add_argument("--start-from", default="vcrelax", choices=STEPS,
                        help="Skip earlier steps (assumes they already ran)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate inputs only, don't run QE")
    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    prefix   = cfg["compound"]["prefix"]
    qe_bin   = cfg["paths"]["qe_bin"]
    outdir   = Path(args.outdir) if args.outdir else Path(__file__).parent / "runs" / prefix

    # Create directory structure
    dirs = {
        "base":    outdir,
        "vcrelax": outdir / "1_vcrelax",
        "phonons": outdir / "2_phonons",
        "epw":     outdir / "3_epw",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    (dirs["vcrelax"] / "tmp").mkdir(exist_ok=True)
    (dirs["phonons"] / "tmp").mkdir(exist_ok=True)
    (dirs["epw"]     / "tmp").mkdir(exist_ok=True)

    banner(f"EPW Pipeline: {prefix} @ {cfg['compound']['pressure_GPa']} GPa\n  Output: {outdir}")

    skip_until = args.start_from
    active = (skip_until == "vcrelax")
    step_n = 0

    def should_run(step_name):
        nonlocal active, step_n
        step_n += 1
        if step_name == skip_until:
            active = True
        return active

    # Step 1: vc-relax
    # Skip if: ibrav != 0 and A is given, OR skip_vcrelax: true in config
    if should_run("vcrelax"):
        s = cfg["structure"]
        skip_relax = (
            cfg.get("skip_vcrelax", False)
            or (s["ibrav"] != 0 and s.get("A") is not None)
        )
        if skip_relax:
            label = f"A = {s['A']} Å" if s.get("A") else "skip_vcrelax=true"
            print(f"  [1/{len(STEPS)}] vcrelax      ─ skipped ({label})")
        else:
            run_vcrelax(cfg, dirs, qe_bin, args.dry_run)

    # Step 2: SCF
    if should_run("scf"):
        run_scf(cfg, dirs, qe_bin, args.dry_run, step_n)

    # Step 3: phonons + pp.py + q2r + matdyn
    if should_run("phonons"):
        run_phonons(cfg, dirs, qe_bin, args.dry_run, step_n)

    # Step 4: NSCF
    if should_run("nscf"):
        run_nscf(cfg, dirs, qe_bin, args.dry_run, step_n)

    # Step 5: EPW pass 1 (Wannierize)
    if should_run("epw1"):
        run_epw1(cfg, dirs, qe_bin, args.dry_run, step_n)

    # Step 6: EPW pass 2 (Eliashberg)
    if should_run("epw2"):
        run_epw2(cfg, dirs, qe_bin, args.dry_run, step_n)

    print(f"\n  All steps complete. Output in: {outdir}\n")


if __name__ == "__main__":
    main()
