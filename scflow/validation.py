"""Physics validation gate: automated checks run on a finished (or partial) run directory.

Each check returns pass / warn / fail / skip. A run with any `fail` is recorded as
`invalid`; otherwise it is `preliminary` until convergence across runs is shown
(see db.confidence_tier).
"""

import itertools
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from . import postprocess as pp
from .eliashberg import A2F, summarize

CM1_PER_MEV = 8.0655439


@dataclass
class Check:
    name: str
    status: str
    message: str
    value: object = None

    def to_dict(self):
        return asdict(self)


@dataclass
class Report:
    checks: list = field(default_factory=list)
    metrics: dict = field(default_factory=dict)

    def add(self, name, status, message, value=None):
        self.checks.append(Check(name, status, message, value))

    @property
    def worst(self):
        order = {"fail": 3, "warn": 2, "pass": 1, "skip": 0}
        return max((c.status for c in self.checks), key=order.get, default="skip")

    def to_dict(self):
        return {"checks": [c.to_dict() for c in self.checks], "metrics": self.metrics, "worst": self.worst}


# ─────────────────────────────────────────────────────────────────────────────
# Structure (no QE output needed — runs before anything is submitted)
# ─────────────────────────────────────────────────────────────────────────────

def lattice_vectors(cfg):
    s = cfg["structure"]
    if s["ibrav"] == 0:
        return np.array(s["cell_parameters_angstrom"], dtype=float)
    a = float(s["A"])
    if s["ibrav"] == 2:
        return a / 2 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]], dtype=float)
    if s["ibrav"] == 3:
        return a / 2 * np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]], dtype=float)
    if s["ibrav"] == 1:
        return a * np.eye(3)
    raise ValueError(f"lattice_vectors: ibrav={s['ibrav']} not supported")


def _frac_positions(cfg):
    pos = cfg["structure"]["positions"]
    return [p[0] for p in pos], np.array([[float(x) for x in p[1:4]] for p in pos])


def min_distances(cfg):
    """{(el_a, el_b): shortest periodic distance in Å}."""
    A = lattice_vectors(cfg)
    syms, frac = _frac_positions(cfg)
    cart = frac @ A
    shifts = np.array(list(itertools.product([-1, 0, 1], repeat=3))) @ A
    out = {}
    for i, j in itertools.combinations_with_replacement(range(len(cart)), 2):
        d = np.linalg.norm(cart[j] + shifts - cart[i], axis=1)
        if i == j:
            d = d[d > 1e-8]
        key = tuple(sorted((syms[i], syms[j])))
        out[key] = min(out.get(key, np.inf), float(d.min()))
    return out


def check_structure(cfg, report):
    from ase.data import atomic_masses, atomic_numbers, covalent_radii

    syms, frac = _frac_positions(cfg)
    for (a, b), d in sorted(min_distances(cfg).items()):
        rsum = covalent_radii[atomic_numbers[a]] + covalent_radii[atomic_numbers[b]]
        ratio = d / rsum
        name = f"structure.distance.{a}-{b}"
        msg = f"shortest {a}-{b} = {d:.3f} Å ({ratio:.2f} × covalent-radius sum {rsum:.2f} Å)"
        if ratio < 0.5:
            report.add(name, "fail", msg + " — atoms overlap; positions or Wyckoff site are almost certainly wrong", d)
        elif ratio < 0.65:
            report.add(name, "warn", msg + " — unusually short even for high pressure", d)
        else:
            report.add(name, "pass", msg, d)

    for sp in cfg["structure"]["species"]:
        ref = atomic_masses[atomic_numbers[sp["symbol"]]]
        dev = abs(sp["mass"] - ref) / ref
        status = "pass" if dev < 0.01 else "warn"
        report.add(f"structure.mass.{sp['symbol']}", status,
                   f"mass {sp['mass']} vs standard {ref:.4f}" + ("" if status == "pass" else " — isotope substitution intended?"),
                   sp["mass"])

    try:
        import spglib
    except ImportError:
        report.add("structure.spacegroup", "skip", "spglib not installed")
        return
    spglib.error.OLD_ERROR_HANDLING = False
    numbers = [atomic_numbers[s] for s in syms]
    ds = spglib.get_symmetry_dataset((lattice_vectors(cfg), frac, numbers), symprec=1e-3)
    n_rot = len({r.tobytes() for r in ds.rotations})
    report.metrics["spacegroup"] = ds.international
    report.metrics["spglib_point_ops"] = n_rot
    expected = cfg["structure"].get("expected_spacegroup")
    if expected and expected.replace(" ", "") != ds.international:
        report.add("structure.spacegroup", "fail",
                   f"spglib finds {ds.international}, config expects {expected}", ds.international)
    else:
        report.add("structure.spacegroup", "pass",
                   f"{ds.international} ({n_rot} point operations)", ds.international)


# ─────────────────────────────────────────────────────────────────────────────
# DFT ground state
# ─────────────────────────────────────────────────────────────────────────────

def check_scf(cfg, scf_out, report):
    if not Path(scf_out).exists():
        report.add("scf.converged", "skip", f"{Path(scf_out).name} not found")
        return
    ok = pp.job_done(scf_out) and pp.scf_converged(scf_out)
    report.add("scf.converged", "pass" if ok else "fail",
               "SCF converged" if ok else "SCF did not converge or did not finish")
    report.metrics["qe_version"] = pp.parse_program_version(scf_out)
    report.metrics["fermi_energy_eV"] = pp.parse_fermi_energy(scf_out)

    nsym = pp.parse_sym_ops(scf_out)
    expected = report.metrics.get("spglib_point_ops")
    if nsym is not None and expected:
        ignored = pp.parse_ignored_sym_ops(scf_out)
        if nsym >= expected:
            report.add("scf.symmetry", "pass", f"QE uses {nsym} symmetry operations (spglib: {expected})", nsym)
        else:
            why = f"; {ignored} dropped for FFT-incommensurate fractional translations" if ignored else ""
            report.add("scf.symmetry", "warn",
                       f"QE uses {nsym} of {expected} operations{why} — phonon irreps and q-point reduction will differ",
                       nsym)

    p_kbar = pp.parse_pressure_kbar(scf_out)
    if p_kbar is None:
        report.add("scf.pressure", "skip", "no stress in SCF output (tstress=.false.)")
        return
    p_gpa = p_kbar / 10.0
    target = float(cfg["compound"]["pressure_GPa"])
    report.metrics["scf_pressure_GPa"] = p_gpa
    dev = abs(p_gpa - target)
    warn_tol = max(5.0, 0.05 * abs(target))
    fail_tol = max(20.0, 0.20 * abs(target))
    msg = f"DFT pressure {p_gpa:.1f} GPa at the simulated geometry; target {target:.1f} GPa"
    if dev > fail_tol:
        report.add("scf.pressure", "fail",
                   msg + " — the structure is not at the stated pressure (wrong lattice or atomic positions)", p_gpa)
    elif dev > warn_tol:
        report.add("scf.pressure", "warn", msg + " — consider vc-relax at the target pressure", p_gpa)
    else:
        report.add("scf.pressure", "pass", msg, p_gpa)


# ─────────────────────────────────────────────────────────────────────────────
# Phonons (coarse DFPT grid)
# ─────────────────────────────────────────────────────────────────────────────

def check_phonons(cfg, ph_out, report, imag_fail_cm1=20.0, asr_warn_cm1=30.0):
    if not Path(ph_out).exists():
        report.add("phonons.stability", "skip", f"{Path(ph_out).name} not found")
        return
    qpts = pp.parse_phonon_qpoints(ph_out)
    if not qpts:
        report.add("phonons.stability", "fail", "no phonon frequencies found in ph.x output")
        return
    worst, worst_q, acoustic = 0.0, None, None
    for q, freqs in qpts:
        f = sorted(freqs)
        if all(abs(x) < 1e-6 for x in q):
            by_abs = sorted(freqs, key=abs)
            acoustic = abs(by_abs[2])
            f = sorted(by_abs[3:])
        if f and f[0] < worst:
            worst, worst_q = f[0], q
    report.metrics["phonon_max_cm1"] = max(max(fr) for _, fr in qpts)
    report.metrics["phonon_min_cm1"] = worst
    if worst < -imag_fail_cm1:
        report.add("phonons.stability", "fail",
                   f"imaginary mode {worst:.1f} cm⁻¹ at q={worst_q}: structure is dynamically unstable in the "
                   f"harmonic approximation, so a harmonic Tc is meaningless. Check the structure/pressure checks "
                   f"first; if those pass, an anharmonic (SSCHA) treatment is required.", worst)
    elif worst < -5.0:
        report.add("phonons.stability", "warn", f"small imaginary mode {worst:.1f} cm⁻¹ at q={worst_q}", worst)
    else:
        report.add("phonons.stability", "pass", f"all {len(qpts)} irreducible q-points stable", worst)
    if acoustic is not None:
        status = "pass" if acoustic < asr_warn_cm1 else "warn"
        report.add("phonons.acoustic_sum_rule", status,
                   f"Γ acoustic modes |ω| ≤ {acoustic:.1f} cm⁻¹ before ASR correction"
                   + ("" if status == "pass" else " — tighten tr2_ph / ecutrho"), acoustic)


# ─────────────────────────────────────────────────────────────────────────────
# Wannier + EPW
# ─────────────────────────────────────────────────────────────────────────────

def check_wannier(cfg, wout, report, spread_warn=10.0):
    if not Path(wout).exists():
        report.add("wannier.spreads", "skip", f"{Path(wout).name} not found")
        return
    spreads, dis_ok = pp.parse_wannier_spreads(wout)
    report.add("wannier.disentanglement", "pass" if dis_ok else "warn",
               "disentanglement converged" if dis_ok else
               "disentanglement hit the iteration limit without converging — band interpolation may be off")
    if spreads:
        big = max(spreads)
        report.metrics["wannier_max_spread_A2"] = big
        report.add("wannier.spreads", "pass" if big < spread_warn else "warn",
                   f"{len(spreads)} WFs, max spread {big:.2f} Å²", big)
        if len(spreads) != cfg["wannier"]["nbndsub"]:
            report.add("wannier.count", "fail", f"{len(spreads)} WFs in .wout but nbndsub={cfg['wannier']['nbndsub']}")


def check_frozen_window(cfg, nscf_out, report):
    if not Path(nscf_out).exists():
        return
    _, eigs = pp.nscf_eigenvalues(nscf_out)
    if eigs is None:
        report.add("wannier.frozen_window_clearance", "skip", "NSCF has no band energies (verbosity='high' needed)")
        return
    w = cfg["wannier"]
    clr = pp.frozen_window_clearance(eigs[:, pp.excluded_band_count(w.get("bands_skipped")):], w["dis_froz_max"]) * 1000
    status = "fail" if clr < 10 else ("warn" if clr < 50 else "pass")
    report.add("wannier.frozen_window_clearance", status,
               f"dis_froz_max is {clr:.0f} meV from the nearest band energy"
               + ("" if status == "pass" else " — an edge on a band makes disentanglement ill-conditioned"), clr)


def check_epw(cfg, epw_out, a2f_file, report):
    if not Path(epw_out).exists() or not pp.job_done(epw_out):
        report.add("epw.completed", "skip", f"{Path(epw_out).name} missing or unfinished")
        return
    report.add("epw.completed", "pass", "EPW pass 2 finished")
    if "NaN" in Path(epw_out).read_text():
        report.add("epw.finite", "fail", "EPW output contains NaN electron-phonon couplings — check the frozen-window "
                   "clearance, then raise wannier.nbnd (the band set must not end inside a degenerate multiplet)")
        return
    report.metrics["epw_version"] = pp.parse_program_version(epw_out)
    muc = float(cfg["eliashberg"]["muc"])

    modes = pp.parse_epw_mode_couplings(epw_out)
    if modes:
        neg = [(w, l) for w, l in modes if w < -0.5]
        frac = len(neg) / len(modes)
        report.metrics["fine_grid_negative_mode_fraction"] = frac
        if frac > 0.05:
            report.add("epw.fine_grid_phonons", "fail",
                       f"{len(neg)}/{len(modes)} interpolated modes imaginary — phonon interpolation is unreliable", frac)
        elif neg:
            report.add("epw.fine_grid_phonons", "warn",
                       f"{len(neg)}/{len(modes)} interpolated modes imaginary (min {min(w for w, _ in neg):.2f} meV); "
                       f"refine the coarse q-grid", frac)
        else:
            report.add("epw.fine_grid_phonons", "pass", f"all {len(modes)} interpolated modes real", 0.0)

    if not a2f_file or not Path(a2f_file).exists():
        report.add("eliashberg.a2f", "skip", "no .a2f file — λ cannot be independently recomputed")
        return

    lam_sum, _ = pp.parse_epw_a2f_summary(epw_out)
    lams = pp.parse_a2f_smearing_lambdas(a2f_file)
    smear = pp.parse_a2f_smearings(a2f_file)
    # α²F at too small a phonon smearing is spiky on a coarse fine-q grid and overestimates λ.
    # Use the narrowest smearing whose λ agrees with EPW's mode-summed λ (the plateau).
    column, plateau = 1, False
    if lam_sum and lams:
        for i, lam in enumerate(lams):
            if abs(lam - lam_sum) / lam_sum < 0.05:
                column, plateau = i + 1, True
                break
        else:
            column = len(lams)
    spec = A2F.from_epw(a2f_file, column=column)
    ours = summarize(spec, muc, wc_mev=float(cfg["epw"]["wscut"]) * 1000.0)
    report.metrics.update({k: float(v) for k, v in ours.items()})
    report.metrics["lambda_epw_mode_sum"] = lam_sum
    report.metrics["a2f_smearing_meV"] = smear[column - 1] if smear else None
    tc_table = pp.parse_epw_ad_tc_table(epw_out)
    if tc_table:
        report.metrics["tc_allen_dynes_epw_K"] = tc_table.get(round(muc, 2))

    if lams:
        dl = abs(ours["lambda"] - lams[column - 1]) / lams[column - 1]
        report.add("eliashberg.cross_check", "pass" if dl < 0.02 else "fail",
                   f"recomputed λ={ours['lambda']:.3f} from α²F matches EPW's integral {lams[column - 1]:.3f}"
                   if dl < 0.02 else f"recomputed λ={ours['lambda']:.3f} but EPW reports {lams[column - 1]:.3f} "
                   f"for the same α²F column — parsing or units problem", dl)
    if lam_sum and lams:
        if plateau:
            report.add("eliashberg.smearing_plateau", "pass",
                       f"α²F at {report.metrics['a2f_smearing_meV']} meV phonon smearing reproduces the mode-summed "
                       f"λ={lam_sum:.3f}", column)
        else:
            report.add("eliashberg.smearing_plateau", "fail",
                       f"no phonon smearing reproduces the mode-summed λ={lam_sum:.3f} within 5% "
                       f"(α²F λ ranges {min(lams):.3f}–{max(lams):.3f}); the fine q-grid is too coarse", column)

    wmax_ev = ours["omega_max_meV"] / 1000.0
    fsthick = float(cfg["epw"]["fsthick"])
    report.add("epw.fsthick", "pass" if fsthick >= 3 * wmax_ev else "warn",
               f"fsthick {fsthick} eV vs ω_max {wmax_ev * 1000:.0f} meV"
               + ("" if fsthick >= 3 * wmax_ev else " — window narrower than 3×ω_max may truncate scattering"),
               fsthick / wmax_ev)
    wscut = float(cfg["epw"]["wscut"])
    report.add("epw.wscut", "pass" if wscut >= 5 * wmax_ev else "warn",
               f"Matsubara cutoff {wscut} eV = {wscut / wmax_ev:.1f}×ω_max" + ("" if wscut >= 5 * wmax_ev else " (want ≥5×)"),
               wscut / wmax_ev)

    if not 0.08 <= muc <= 0.16:
        report.add("eliashberg.muc", "warn",
                   f"μ*={muc} outside the conventional 0.08–0.16 range — a μ* chosen to match experiment makes the "
                   f"Tc a fit, not a prediction", muc)
    else:
        report.add("eliashberg.muc", "pass", f"μ*={muc}", muc)

    if ours["lambda"] > 1.5:
        report.add("eliashberg.strong_coupling", "pass",
                   f"λ={ours['lambda']:.2f} > 1.5: Allen-Dynes without f1·f2 underestimates Tc here; "
                   f"the isotropic Eliashberg value ({ours['tc_eliashberg_iso_K']:.1f} K) is the better estimate", ours["lambda"])


def check_reference(cfg, report):
    ref = cfg.get("reference") or {}
    tc_ref = ref.get("tc_exp_K") or ref.get("tc_calc_K")
    tc = report.metrics.get("tc_eliashberg_iso_K")
    if not tc_ref or tc is None:
        return
    rel = (tc - tc_ref) / tc_ref
    report.metrics["tc_reference_K"] = tc_ref
    report.metrics["tc_reference_rel_dev"] = rel
    kind = "experiment" if ref.get("tc_exp_K") else "published calculation"
    p = float(cfg["compound"]["pressure_GPa"])
    p_ref = float(ref.get("pressure_GPa", p))
    at = f" at {p_ref:g} GPa" if p_ref != p else ""
    same_pressure = abs(p_ref - p) <= max(10.0, 0.10 * p)
    msg = f"Eliashberg Tc {tc:.1f} K vs {kind} {tc_ref} K{at} ({rel * 100:+.0f}%)"
    if not same_pressure:
        report.add("reference.tc", "warn", msg + f" — reference is at a different pressure than {p:g} GPa, "
                   "so agreement is not evidence", rel)
    else:
        report.add("reference.tc", "pass" if abs(rel) <= 0.25 else "warn", msg, rel)


def validate(cfg, run_dir):
    run_dir = Path(run_dir)
    prefix = cfg["compound"]["prefix"]
    report = Report()
    check_structure(cfg, report)
    check_scf(cfg, run_dir / "2_phonons" / "scf.out", report)
    check_phonons(cfg, run_dir / "2_phonons" / "ph.out", report)
    check_frozen_window(cfg, run_dir / "3_epw" / "nscf.out", report)
    check_wannier(cfg, run_dir / "3_epw" / f"{prefix}.wout", report)
    a2f = sorted((run_dir / "3_epw").glob(f"{prefix}.a2f.0*")) or sorted((run_dir / "3_epw").glob(f"{prefix}.a2f"))
    check_epw(cfg, run_dir / "3_epw" / "epw_pass2.out", a2f[0] if a2f else None, report)
    check_reference(cfg, report)
    return report


def format_report(report, show_skipped=False):
    icon = {"pass": "✓", "warn": "⚠", "fail": "✗", "skip": "·"}
    shown = [c for c in report.checks if show_skipped or c.status != "skip"]
    lines = [f"  {icon[c.status]} {c.name:<34} {c.message}" for c in shown]
    skipped = len(report.checks) - len(shown)
    if skipped:
        lines.append(f"  · {skipped} check(s) skipped — outputs for those steps not present yet")
    m = report.metrics
    if "lambda" in m:
        lines.append("")
        lines.append(f"  λ = {m['lambda']:.3f}   ω_log = {m['omega_log_meV']:.2f} meV   ω₂ = {m['omega2_meV']:.2f} meV")
        lines.append(f"  Tc  Allen-Dynes {m['tc_allen_dynes_K']:.1f} K   AD·f1f2 {m['tc_allen_dynes_f1f2_K']:.1f} K   "
                     f"Eliashberg (isotropic) {m['tc_eliashberg_iso_K']:.1f} K")
    return "\n".join(lines)
