"""
Parse QE/EPW output files to extract key results and enable checkpointing.
"""

import re
from pathlib import Path

import numpy as np


def job_done(outfile):
    """Return True if a QE/EPW output file completed successfully.
    pw.x/ph.x/q2r.x/matdyn.x print 'JOB DONE'; epw.x prints 'Total program execution'."""
    p = Path(outfile)
    if not p.exists():
        return False
    text = p.read_text()
    return "JOB DONE" in text or "Total program execution" in text


BOHR_ANG = 0.529177210903


def parse_relaxed_structure(outfile):
    """Final cell (Å, rows = lattice vectors) and crystal positions from a vc-relax output.

    Returns (cell, [[symbol, x, y, z], ...]) or (None, None) if the relaxation did not finish.
    """
    content = Path(outfile).read_text()
    m = re.search(r'Begin final coordinates(.*?)End final coordinates', content, re.DOTALL)
    if not m:
        return None, None
    block = m.group(1)
    cp = re.search(r'CELL_PARAMETERS \((alat=\s*([\d.]+)|angstrom|bohr)\)\s*\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})', block)
    if not cp:
        return None, None
    rows = [[float(x) for x in line.split()] for line in cp.group(3).strip().splitlines()]
    if cp.group(2):
        scale = float(cp.group(2)) * BOHR_ANG
    elif cp.group(1) == "bohr":
        scale = BOHR_ANG
    else:
        scale = 1.0
    cell = [[x * scale for x in r] for r in rows]
    ap = re.search(r'ATOMIC_POSITIONS \(crystal\)\s*\n(.*?)(?:\n\s*\n|End final|$)', block + "\n\n", re.DOTALL)
    positions = None
    if ap:
        positions = [[t[0]] + [float(x) for x in t[1:4]] for t in (l.split() for l in ap.group(1).strip().splitlines()) if len(t) >= 4]
    return cell, positions


def parse_fermi_energy(outfile):
    """Extract Fermi energy in eV from SCF/NSCF output."""
    content = Path(outfile).read_text()
    m = re.search(r'the Fermi energy is\s+([\d.+-]+)\s+ev', content)
    return float(m.group(1)) if m else None


def parse_sym_ops(outfile):
    """Extract number of symmetry operations from SCF output."""
    content = Path(outfile).read_text()
    m = re.search(r'(\d+)\s+Sym\. Ops\.', content)
    return int(m.group(1)) if m else None


def parse_imaginary_modes(ph_outfile):
    """
    Return list of imaginary mode frequencies (cm-1) from ph.x output.
    Empty list means fully stable.
    """
    content = Path(ph_outfile).read_text()
    imag = re.findall(r'freq\s*\(\s*\d+\)\s*=\s*([-\d.]+)\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]', content)
    return [float(f) for _, f in imag if float(f) < -1.0]


def parse_allen_dynes_tc(epw_outfile):
    """Extract Allen-Dynes Tc and omega_log from EPW pass 2 output.

    Handles two EPW output formats:
      Eliashberg path: 'Estimated Allen-Dynes Tc = X K'
                       'Estimated w_log in Allen-Dynes Tc = X meV'
      a2f path:        'Estimated Allen-Dynes Tc\\n  logavg = X  ...\\n  mu = 0.10 Tc = Y K'
                       logavg is in Ry; convert to meV (1 Ry = 13605.7 meV)
    """
    content = Path(epw_outfile).read_text()
    # Eliashberg path: single-line "Allen-Dynes Tc = X K"
    tc_m = re.search(r'(?:Estimated )?Allen-Dynes Tc\s*=\s*([\d.]+)\s*K', content)
    # a2f path: "mu =  0.10 Tc =  X K" (first mu* entry)
    if not tc_m:
        tc_m = re.search(r'mu\s*=\s*[\d.]+\s+Tc\s*=\s*([\d.]+)\s*K', content)
    # omega_log: explicit meV line (Eliashberg path)
    wlog_m = re.search(
        r'(?:Estimated w_log in Allen-Dynes Tc|omega_log \(a2f\))\s*=\s*([\d.]+)\s*meV',
        content)
    # a2f path: logavg in Ry — convert to meV
    if not wlog_m:
        logavg_m = re.search(r'logavg\s*=\s*([\d.eE+-]+)', content)
        if logavg_m:
            wlog_ry = float(logavg_m.group(1))
            wlog = wlog_ry * 13605.7  # Ry → meV
            return float(tc_m.group(1)) if tc_m else None, wlog
    tc   = float(tc_m.group(1))   if tc_m   else None
    wlog = float(wlog_m.group(1)) if wlog_m else None
    return tc, wlog


def parse_lambda(epw_outfile):
    """Extract e-ph coupling λ from EPW pass 2 output.

    a2f path:        'l_a2f =  X'  (from logavg line)
    Eliashberg path: 'Electron-phonon coupling strength =  X'
    """
    content = Path(epw_outfile).read_text()
    m = re.search(r'l_a2f\s*=\s*([\d.eE+-]+)', content)
    if m:
        return float(m.group(1))
    m = re.search(r'Electron-phonon coupling strength\s*=\s*([\d.eE+-]+)', content)
    return float(m.group(1)) if m else None


def excluded_band_count(bands_skipped):
    """Number of bands removed by an EPW `bands_skipped = 'exclude_bands = 1-4,7'` string."""
    if not bands_skipped:
        return 0
    spec = bands_skipped.split("=", 1)[1] if "=" in bands_skipped else bands_skipped
    n = 0
    for part in spec.replace("'", "").replace('"', "").split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-")
            n += int(hi) - int(lo) + 1
        elif part:
            n += 1
    return n


def nscf_eigenvalues(nscf_outfile):
    """(Ef, array[nk, nbnd] in eV) from an NSCF run with verbosity='high', or (Ef, None)."""
    content = Path(nscf_outfile).read_text()
    ef = re.search(r'the Fermi energy is\s+([\d.+-]+)\s+ev', content)
    blocks = re.findall(r'bands \(ev\):\s*([\d\s.+-]+?)(?=\n\s*\n|\n\s*k =)', content, re.DOTALL)
    if not blocks:
        return (float(ef.group(1)) if ef else None), None
    rows = [list(map(float, b.split())) for b in blocks]
    width = min(len(r) for r in rows)
    return (float(ef.group(1)) if ef else None), np.array([r[:width] for r in rows])


def frozen_window_clearance(eigs, froz_max):
    """Smallest |E_nk − dis_froz_max| in eV. A frozen-window edge sitting on a band (often a
    degenerate pair) makes disentanglement ill-conditioned and can turn every e-ph element into NaN."""
    return float(np.abs(eigs - froz_max).min())


def scan_frozen_window(nscf_outfile, nbndsub, dis_win_min=None, n_excluded=0):
    """Print how many bands each candidate dis_froz_max would freeze, and suggest the value with the
    largest clearance from every band energy among those that freeze ≤ nbndsub bands.

    Returns (Ef, suggested_froz_max, clearance_eV) or (Ef, None, None).
    """
    Ef, eigs = nscf_eigenvalues(nscf_outfile)
    if eigs is None:
        print("\n  ⚠ Window scan: no band energies in the NSCF output (needs verbosity='high' for >100 k-points)")
        return Ef, None, None
    eigs = eigs[:, n_excluded:]
    if dis_win_min is not None:
        eigs = np.where(eigs >= dis_win_min, eigs, np.inf)

    def frozen(e_test):
        return int((eigs < e_test).sum(axis=1).max())

    print(f"\n  Window scan (NSCF coarse grid, nbndsub={nbndsub}, {n_excluded} bands excluded):")
    print(f"  {'dis_froz_max':>14}   max bands   clearance   status")
    print(f"  {'─' * 56}")
    for e_test in [Ef + x for x in [-1, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10]]:
        n = frozen(e_test)
        clr = frozen_window_clearance(eigs[np.isfinite(eigs)], e_test) * 1000
        print(f"  {e_test:>14.3f} eV   {n:>9}   {clr:>6.0f} meV   {'✓' if n <= nbndsub else '✗ EXCEEDS'}")

    finite = eigs[np.isfinite(eigs)]
    best = None
    for e_test in np.arange(Ef + 0.5, Ef + 10.0, 0.01):
        if frozen(e_test) > nbndsub:
            break
        clr = frozen_window_clearance(finite, e_test)
        if best is None or clr > best[1]:
            best = (round(float(e_test), 2), clr)
    print(f"\n  Coarse Ef  = {Ef:.3f} eV")
    if best is None:
        return Ef, None, None
    print(f"  Suggested dis_froz_max = {best[0]:.2f} eV (largest clearance from any band: {best[1] * 1000:.0f} meV)")
    return Ef, best[0], best[1]


def parse_nqirr(ph_outfile):
    """Count number of irreducible q-points from ph.x output."""
    content = Path(ph_outfile).read_text()
    matches = re.findall(r'Calculation of q =', content)
    return len(matches)


# ─────────────────────────────────────────────────────────────────────────────
# Parsers used by the validation gate
# ─────────────────────────────────────────────────────────────────────────────

def parse_pressure_kbar(outfile):
    """Last 'P=' value (kbar) from a pw.x run with tstress=.true."""
    hits = re.findall(r'\(kbar\)\s+P=\s*([-\d.]+)', Path(outfile).read_text())
    return float(hits[-1]) if hits else None


def scf_converged(outfile):
    text = Path(outfile).read_text()
    return "convergence has been achieved" in text and "convergence NOT achieved" not in text


def parse_program_version(outfile):
    m = re.search(r'Program (\S+) v\.(\S+) starts', Path(outfile).read_text())
    return f"{m.group(1)} {m.group(2)}" if m else None


def parse_ignored_sym_ops(outfile):
    m = re.search(r'(\d+) additional sym\.ops\. were found but ignored', Path(outfile).read_text())
    return int(m.group(1)) if m else 0


def parse_phonon_qpoints(ph_outfile):
    """[(q_cart_2pi/a, [freq_cm-1, ...]), ...] for every q-point diagonalized by ph.x."""
    blocks = Path(ph_outfile).read_text().split("Diagonalizing the dynamical matrix")[1:]
    out = []
    for b in blocks:
        q = re.search(r'q = \(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\)', b)
        freqs = re.findall(r'freq\s*\(\s*\d+\)\s*=\s*[-\d.]+\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]', b)
        if q and freqs:
            out.append((tuple(float(x) for x in q.groups()), [float(f) for f in freqs]))
    return out


def parse_wannier_spreads(wout):
    """(final spreads in Å², disentanglement_converged) from a Wannier90 .wout file."""
    text = Path(wout).read_text()
    final = text.rsplit("Final State", 1)
    spreads = []
    if len(final) == 2:
        spreads = [float(x) for x in re.findall(
            r'WF centre and spread\s+\d+\s+\([^)]*\)\s+([\d.]+)', final[1].split("Sum of centres")[0])]
    # .wout accumulates across reruns; judge only the last disentanglement.
    last_dis = text.rsplit("Extraction of optimally-connected subspace", 1)[-1]
    dis_ok = "Disentanglement convergence criteria not satisfied" not in last_dis
    return spreads, dis_ok


def parse_epw_mode_couplings(epw_outfile):
    """[(omega_meV, lambda_qν), ...] for every fine-grid phonon mode printed by phonselfen."""
    pairs = re.findall(
        r'lambda___\(\s*\d+\s*\)=\s*([-\d.Ee+]+)\s+gamma___=\s*[-\d.Ee+]+\s*meV\s+omega=\s*([-\d.Ee+]+)\s*meV',
        Path(epw_outfile).read_text())
    return [(float(w), float(l)) for l, w in pairs]


def parse_epw_a2f_summary(epw_outfile):
    """λ and ω_log (meV) as printed by EPW's a2f / Eliashberg summary, or (None, None)."""
    text = Path(epw_outfile).read_text()
    lam = re.search(r'lambda\s*:\s*([\d.Ee+-]+)', text) or re.search(r'l_a2f\s*=\s*([\d.Ee+-]+)', text)
    logavg = re.search(r'logavg\s*=\s*([\d.Ee+-]+)', text)
    wlog = float(logavg.group(1)) * 13605.693 if logavg else None
    if wlog is None:
        m = re.search(r'Estimated w_log in Allen-Dynes Tc\s*=\s*([\d.]+)\s*meV', text)
        wlog = float(m.group(1)) if m else None
    return (float(lam.group(1)) if lam else None), wlog


def parse_epw_ad_tc_table(epw_outfile):
    """{mu*: Tc_K} from EPW's a2f-path Allen-Dynes table."""
    return {float(m): float(t) for m, t in
            re.findall(r'mu\s*=\s*([\d.]+)\s+Tc\s*=\s*([\d.]+)\s*K', Path(epw_outfile).read_text())}


def parse_a2f_smearings(a2f_file):
    """Phonon smearing (meV) of each column, from the footer of a .a2f file."""
    lines = Path(a2f_file).read_text().splitlines()
    for i, line in enumerate(lines):
        if "Phonon smearing" in line and i + 1 < len(lines):
            return [float(x) for x in lines[i + 1].replace("#", "").split()]
    return []


def parse_a2f_smearing_lambdas(a2f_file):
    """λ for each phonon smearing column from the 'Integrated el-ph coupling' footer of a .a2f file."""
    lines = Path(a2f_file).read_text().splitlines()
    for i, line in enumerate(lines):
        if "Integrated el-ph coupling" in line and i + 1 < len(lines):
            return [float(x) for x in lines[i + 1].replace("#", "").split()]
    return []
