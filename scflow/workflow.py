"""Structure → Eliashberg Tc workflow: vc-relax → SCF → DFPT phonons → NSCF → EPW (Wannier) → EPW (Tc).

Every step is checkpointed on its own output file, so re-running skips finished steps.
"""

import json
import platform
import shutil
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

from . import __version__
from .executors import StepFailed
from .generators import (gen_epw_pass1, gen_epw_pass2, gen_matdyn, gen_nscf, gen_phonon,
                         gen_q2r, gen_scf, gen_vcrelax)
from .postprocess import (excluded_band_count, frozen_window_clearance, job_done, nscf_eigenvalues,
                          parse_allen_dynes_tc, parse_fermi_energy, parse_imaginary_modes, parse_nqirr,
                          parse_relaxed_structure, parse_sym_ops, scan_frozen_window)

STEPS = ["vcrelax", "scf", "phonons", "nscf", "epw1", "epw2"]


def largest_divisor(n, max_val):
    for d in range(min(n, max_val), 0, -1):
        if n % d == 0:
            return d
    return 1


class Workflow:
    def __init__(self, cfg, outdir, executor, dry_run=False, log=print):
        self.cfg = cfg
        self.prefix = cfg["compound"]["prefix"]
        self.ex = executor
        self.dry_run = dry_run
        self.log = log
        self.nproc = cfg["parallel"]["nproc"]
        base = Path(outdir).resolve()
        self.dirs = {"base": base, "vcrelax": base / "1_vcrelax", "phonons": base / "2_phonons", "epw": base / "3_epw"}
        for key in ("vcrelax", "phonons", "epw"):
            (self.dirs[key] / "tmp").mkdir(parents=True, exist_ok=True)
        self.manifest_path = base / "manifest.json"
        self.manifest = self._load_manifest()

    # ── bookkeeping ──────────────────────────────────────────────────────────

    def _load_manifest(self):
        if self.manifest_path.exists():
            return json.loads(self.manifest_path.read_text())
        return {
            "scflow_version": __version__,
            "created": datetime.now(timezone.utc).isoformat(),
            "host": socket.gethostname(),
            "platform": platform.platform(),
            "cpu": platform.processor() or platform.machine(),
            "steps": {},
        }

    def _save_manifest(self):
        self.manifest["nproc"] = self.nproc
        self.manifest["config"] = {k: v for k, v in self.cfg.items() if k != "paths"}
        self.manifest_path.write_text(json.dumps(self.manifest, indent=2, default=str))

    def _status(self, step, msg):
        n = STEPS.index(step) + 1 if step in STEPS else "-"
        self.log(f"  [{n}/{len(STEPS)}] {step:<12} {msg}")

    def _run(self, step, binary, infile, outfile, cwd, nproc=1, extra_args=None):
        t0 = time.time()
        self.ex.run_qe(binary, infile, outfile, cwd, nproc=nproc, extra_args=extra_args)
        rec = self.manifest["steps"].setdefault(step, {})
        rec[binary] = {"seconds": round(time.time() - t0, 1), "nproc": nproc,
                       "finished": datetime.now(timezone.utc).isoformat()}
        self._save_manifest()

    def _cached(self, step, infile, text, outfile, if_changed):
        """True if `outfile` is a finished run of exactly `text`; raise if it ran with different input."""
        infile = Path(infile)
        if not job_done(outfile):
            return False
        if infile.exists() and infile.read_text() == text:
            return True
        raise StepFailed(f"{step}: {infile.name} differs from the run in {outfile.parent}; {if_changed}")

    # ── steps ────────────────────────────────────────────────────────────────

    def vcrelax(self):
        if self.cfg.get("skip_vcrelax", False):
            self._status("vcrelax", "─ skipped (skip_vcrelax: true; using the geometry in the config)")
            return
        d = self.dirs["vcrelax"]
        out = d / "vcrelax.out"
        text = gen_vcrelax(self.cfg)
        if not self._cached("vcrelax", d / "vcrelax.in", text, out, "everything downstream depends on it — use a new --outdir"):
            (d / "vcrelax.in").write_text(text)
            if self.dry_run:
                self._status("vcrelax", "(dry-run) input written")
                return
            self._status("vcrelax", f"running at {self.cfg['compound']['pressure_GPa']} GPa...")
            self._run("vcrelax", "pw.x", d / "vcrelax.in", out, d, self.nproc)
        self._apply_relaxed_structure(out)

    def _apply_relaxed_structure(self, out):
        cell, positions = parse_relaxed_structure(out)
        if cell is None:
            raise StepFailed(f"vc-relax finished without final coordinates; see {out}")
        s = self.cfg["structure"]
        if s["ibrav"] in (2, 3):
            s["A"] = round(2.0 * max(abs(x) for x in cell[0]), 6)
            msg = f"A = {s['A']:.4f} Å"
        elif s["ibrav"] == 0:
            s["cell_parameters_angstrom"] = cell
            msg = "cell + positions updated"
        else:
            raise StepFailed(f"propagating relaxed geometry for ibrav={s['ibrav']} is not implemented")
        if positions:
            s["positions"] = positions
        (self.dirs["base"] / "relaxed_structure.yaml").write_text(yaml.safe_dump({"structure": s}, sort_keys=False))
        self._status("vcrelax", f"✓ {msg}")

    def scf(self):
        d = self.dirs["phonons"]
        out = d / "scf.out"
        text = gen_scf(self.cfg)
        if self._cached("scf", d / "scf.in", text, out, "everything downstream depends on it — use a new --outdir"):
            self._status("scf", f"✓ (cached) Ef={parse_fermi_energy(out):.2f} eV, {parse_sym_ops(out)} sym ops")
            return
        (d / "scf.in").write_text(text)
        if self.dry_run:
            self._status("scf", "(dry-run) input written")
            return
        self._status("scf", "running...")
        self._run("scf", "pw.x", d / "scf.in", out, d, self.nproc)
        self._status("scf", f"✓ Ef={parse_fermi_energy(out):.2f} eV, {parse_sym_ops(out)} sym ops")

    def phonons(self):
        d = self.dirs["phonons"]
        ph_out, q2r_out, mat_out = d / "ph.out", d / "q2r.out", d / "matdyn.out"
        text = gen_phonon(self.cfg)
        if not self._cached("phonons", d / "ph.in", text, ph_out, "everything downstream depends on it — use a new --outdir"):
            (d / "ph.in").write_text(text)
            (d / "q2r.in").write_text(gen_q2r(self.cfg))
            (d / "matdyn.in").write_text(gen_matdyn(self.cfg))
            if self.dry_run:
                self._status("phonons", "(dry-run) inputs written")
                return
            self._status("phonons", "running ph.x (the slow step)...")
            self._run("phonons", "ph.x", d / "ph.in", ph_out, d, self.nproc)
        imag = parse_imaginary_modes(ph_out)
        flag = "stable" if not imag else f"⚠ {len(imag)} imaginary modes (min {min(imag):.1f} cm⁻¹)"
        self._status("phonons", f"✓ {parse_nqirr(ph_out)} irreducible q-points, {flag}")

        if not (d / "tmp" / "save").exists():
            self._run_pp_py()
            link = d / "save"
            if not link.exists():
                link.symlink_to((d / "tmp" / "save").resolve())
        if not job_done(q2r_out):
            self._run("phonons", "q2r.x", d / "q2r.in", q2r_out, d)
        if not job_done(mat_out):
            self._run("phonons", "matdyn.x", d / "matdyn.in", mat_out, d)

    def _run_pp_py(self):
        """pp.py collects dvscf + dyn files into save/ for EPW; it expects them under tmp/_ph0/."""
        d, tmp = self.dirs["phonons"], self.dirs["phonons"] / "tmp"
        n_dyn = len(list(d.glob(f"{self.prefix}.dyn[0-9]*")))
        for i in range(n_dyn + 1):
            src, dst = (d / f"{self.prefix}.dyn{i}").resolve(), tmp / f"{self.prefix}.dyn{i}"
            if src.exists() and not dst.exists():
                dst.symlink_to(src)
        ph0 = tmp / "_ph0"
        ph0.mkdir(exist_ok=True)
        for src in tmp.iterdir():
            if src.name.startswith(self.prefix) and src.name != f"{self.prefix}.save":
                dst = ph0 / src.name
                if not dst.exists():
                    dst.symlink_to(src.resolve())
        proc = self.ex.run_pp_py(self.prefix, tmp)
        if not (tmp / "save").exists():
            raise StepFailed(f"pp.py failed:\n{proc.stderr}")
        self._status("phonons", f"✓ pp.py collected {n_dyn} dvscf files into save/")

    def nscf(self):
        d = self.dirs["epw"]
        out = d / "nscf.out"
        scf_save = self.dirs["phonons"] / "tmp" / f"{self.prefix}.save"
        epw_save = d / "tmp" / f"{self.prefix}.save"
        if not epw_save.exists() and scf_save.exists():
            shutil.copytree(scf_save, epw_save)
        text = gen_nscf(self.cfg)
        if not self._cached("nscf", d / "nscf.in", text, out,
                            "delete nscf.out and rerun with --start-from nscf --reset-epw"):
            (d / "nscf.in").write_text(text)
            if self.dry_run:
                self._status("nscf", "(dry-run) input written")
                return
            self._status("nscf", "running...")
            self._run("nscf", "pw.x", d / "nscf.in", out, d, self.nproc)
        self._status("nscf", f"✓ Ef={parse_fermi_energy(out):.2f} eV")
        w = self.cfg["wannier"]
        scan_frozen_window(out, w["nbndsub"], w.get("dis_win_min"), excluded_band_count(w.get("bands_skipped")))

    def _check_frozen_window(self):
        nscf_out = self.dirs["epw"] / "nscf.out"
        if not job_done(nscf_out):
            return
        _, eigs = nscf_eigenvalues(nscf_out)
        if eigs is None:
            return
        w = self.cfg["wannier"]
        clr = frozen_window_clearance(eigs[:, excluded_band_count(w.get("bands_skipped")):], w["dis_froz_max"])
        if clr < 0.010:
            raise StepFailed(f"epw1: dis_froz_max = {w['dis_froz_max']} eV lies {clr * 1000:.1f} meV from a band energy; "
                             f"disentanglement becomes ill-conditioned (NaN couplings). Use the value suggested by "
                             f"the NSCF window scan.")

    def reset_epw(self):
        d = self.dirs["epw"]
        doomed = [d / "epw_pass1.out", d / "epw_pass2.out"]
        doomed += list((d / "tmp").glob(f"{self.prefix}.epmatwp")) + list((d / "tmp").glob(f"{self.prefix}.epb*"))
        doomed += list((d / "tmp").glob(f"{self.prefix}.ephmat"))
        for p in doomed:
            if p.is_dir():
                shutil.rmtree(p)
            elif p.exists():
                p.unlink()
        self.log(f"  reset-epw: removed {len(doomed)} stale EPW files")

    def epw1(self, nproc):
        d = self.dirs["epw"]
        epmatwp = d / "tmp" / f"{self.prefix}.epmatwp"
        text = gen_epw_pass1(self.cfg)
        if epmatwp.exists():
            if (d / "epw_pass1.in").exists() and (d / "epw_pass1.in").read_text() != text:
                raise StepFailed("epw1: Wannier/EPW pass-1 settings changed since .epmatwp was written; "
                                 "rerun with --reset-epw")
            self._status("epw1", f"✓ (cached) .epmatwp = {epmatwp.stat().st_size / 1e6:.0f} MB")
            return
        (d / "epw_pass1.in").write_text(text)
        if self.dry_run:
            self._status("epw1", "(dry-run) input written")
            return
        self._check_frozen_window()
        self._status("epw1", "Wannierizing + writing e-ph matrices...")
        # EPW 5.9 requires one MPI task per k-point pool.
        self._run("epw1", "epw.x", d / "epw_pass1.in", d / "epw_pass1.out", d, nproc,
                  ["-npool", str(nproc)] if nproc > 1 else None)
        if not epmatwp.exists():
            raise StepFailed(f"EPW pass 1 finished but wrote no {epmatwp.name}")
        self._status("epw1", f"✓ .epmatwp = {epmatwp.stat().st_size / 1e6:.0f} MB")

    def epw2(self, nproc):
        d = self.dirs["epw"]
        out = d / "epw_pass2.out"
        # EPW writes into <prefix>.ephmat/ before creating it.
        (d / "tmp" / f"{self.prefix}.ephmat").mkdir(exist_ok=True)
        text = gen_epw_pass2(self.cfg)
        # Pass 2 only reads pass-1 data, so a changed pass-2 input (fine grids, μ*, ...) just reruns it.
        fresh = job_done(out) and (d / "epw_pass2.in").exists() and (d / "epw_pass2.in").read_text() == text
        (d / "epw_pass2.in").write_text(text)
        if not fresh:
            if self.dry_run:
                self._status("epw2", "(dry-run) input written")
                return
            nkf = self.cfg["epw"]["fine_kgrid"]
            nqf = self.cfg["epw"].get("fine_qgrid", nkf)
            self._status("epw2", f"fine grids k {nkf[0]}×{nkf[1]}×{nkf[2]}, q {nqf[0]}×{nqf[1]}×{nqf[2]} + Tc...")
            self._run("epw2", "epw.x", d / "epw_pass2.in", out, d, nproc,
                      ["-npool", str(nproc)] if nproc > 1 else None)
        tc, wlog = parse_allen_dynes_tc(out)
        if tc is not None:
            self._status("epw2", f"✓ EPW Allen-Dynes Tc = {tc:.1f} K (μ*={self.cfg['eliashberg']['muc']})"
                         + (f", ω_log = {wlog:.1f} meV" if wlog else ""))

    # ── driver ───────────────────────────────────────────────────────────────

    def run(self, start_from="vcrelax", stop_after=None, reset_epw=False):
        coarse = self.cfg["dft"].get("epw_coarse_kgrid", self.cfg["dft"]["coarse_kgrid"])
        fine = self.cfg["epw"]["fine_kgrid"]
        epw1_n = largest_divisor(coarse[0] * coarse[1] * coarse[2], self.nproc)
        epw2_n = largest_divisor(fine[0] * fine[1] * fine[2], self.nproc)
        self.log(f"\n  scflow {__version__}: {self.prefix} @ {self.cfg['compound']['pressure_GPa']} GPa"
                 f"\n  output: {self.dirs['base']}\n  MPI tasks: {self.nproc} (EPW pass1={epw1_n}, pass2={epw2_n})\n")
        if reset_epw:
            self.reset_epw()
        relaxed = self.dirs["base"] / "relaxed_structure.yaml"
        if start_from != "vcrelax" and relaxed.exists():
            self.cfg["structure"] = yaml.safe_load(relaxed.read_text())["structure"]
            self.log("  using relaxed geometry from relaxed_structure.yaml")
        actions = {
            "vcrelax": self.vcrelax, "scf": self.scf, "phonons": self.phonons, "nscf": self.nscf,
            "epw1": lambda: self.epw1(epw1_n), "epw2": lambda: self.epw2(epw2_n),
        }
        first = STEPS.index(start_from)
        last = STEPS.index(stop_after) if stop_after else len(STEPS) - 1
        for step in STEPS[first:last + 1]:
            actions[step]()
        self._save_manifest()
