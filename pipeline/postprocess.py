"""
Parse QE/EPW output files to extract key results and enable checkpointing.
"""

import re
from pathlib import Path


def job_done(outfile):
    """Return True if a QE/EPW output file completed successfully.
    pw.x/ph.x/q2r.x/matdyn.x print 'JOB DONE'; epw.x prints 'Total program execution'."""
    p = Path(outfile)
    if not p.exists():
        return False
    text = p.read_text()
    return "JOB DONE" in text or "Total program execution" in text


def parse_vcrelax_a(outfile, ibrav):
    """
    Extract the relaxed lattice parameter from vc-relax output.
    Returns A in Angstrom for ibrav=3 (BCC), or (cell_params, positions) for ibrav=0.
    """
    content = Path(outfile).read_text()

    if ibrav == 3:
        # For ibrav=3: CELL_PARAMETERS block has vectors as fractions of initial alat.
        # QE reports "lattice parameter (alat) = X a.u." which is the initial value.
        # The new primitive vector is component * alat_initial.
        # For BCC: a1 = A/2*(-1,1,1) in Cartesian. In units of alat, a1=(−0.5,0.5,0.5).
        # After relax: a1 = (s,s,s) in units of initial alat.
        # New A = 2*s*alat_initial (Bohr) * 0.529177 (Å/Bohr)
        alat_match = re.search(r'lattice parameter \(alat\)\s*=\s*([\d.]+)\s*a\.u\.', content)
        cp_match = re.search(
            r'CELL_PARAMETERS \(alat=\s*[\d.]+\)\s*\n'
            r'\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)',
            content
        )
        if alat_match and cp_match:
            alat = float(alat_match.group(1))       # Bohr (initial)
            scale = float(cp_match.group(1))        # e.g. 0.499346 for the x-component
            # For ibrav=3, a1=(s,s,s)*alat. Conventional cube edge A = 2*s*alat (Bohr)
            A_bohr = 2.0 * scale * alat
            A_ang = A_bohr * 0.529177
            return round(A_ang, 6)
        return None

    elif ibrav == 0:
        # Extract final CELL_PARAMETERS and ATOMIC_POSITIONS from output
        final_match = re.search(
            r'Begin final coordinates(.*?)End final coordinates',
            content, re.DOTALL
        )
        if not final_match:
            return None
        block = final_match.group(1)

        # Cell parameters (in angstrom or alat — QE prints in alat here, need to convert)
        # Actually QE prints them in units of alat in vc-relax final output
        # Better: look for CELL_PARAMETERS block and the alat header
        cp_match = re.search(
            r'CELL_PARAMETERS \(alat=\s*([\d.]+)\)\s*\n'
            r'\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s*\n'
            r'\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s*\n'
            r'\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)',
            block
        )
        if cp_match:
            alat = float(cp_match.group(1))
            scale = alat * 0.529177  # Bohr → Ang
            cell = [
                [float(cp_match.group(2))*scale, float(cp_match.group(3))*scale, float(cp_match.group(4))*scale],
                [float(cp_match.group(5))*scale, float(cp_match.group(6))*scale, float(cp_match.group(7))*scale],
                [float(cp_match.group(8))*scale, float(cp_match.group(9))*scale, float(cp_match.group(10))*scale],
            ]
            return cell
        return None

    return None


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


def scan_frozen_window(nscf_outfile, nbndsub, dis_win_min=None):
    """Scan NSCF eigenvalues to find the highest dis_froz_max that keeps ≤ nbndsub
    bands in the frozen window [dis_win_min, e_test] at every k-point.

    dis_win_min: only count bands above this energy (matches the EPW window floor).
                 If None, counts all bands below e_test (old behaviour — can overcount
                 when deep semicore bands are in the NSCF but below the EPW window).

    Returns (Ef, suggested_froz_max, gap_to_next_band) or (Ef, None, None) on failure.
    Prints a summary table for the user.
    """
    content = Path(nscf_outfile).read_text()
    Ef_m = re.search(r'the Fermi energy is\s+([\d.+-]+)\s+ev', content)
    if not Ef_m:
        return None, None, None
    Ef = float(Ef_m.group(1))

    blocks = re.findall(
        r'bands \(ev\):\s*([\d\s.+-]+?)(?=\n\s*\n|\n\s*k =)',
        content, re.DOTALL
    )
    if not blocks:
        return Ef, None, None

    def count_bands(block, e_test):
        evals = map(float, block.split())
        if dis_win_min is not None:
            return sum(1 for e in evals if dis_win_min <= e < e_test)
        return sum(1 for e in evals if e < e_test)

    print(f"\n  Window scan (NSCF coarse grid, nbndsub={nbndsub}"
          + (f", win_min={dis_win_min:.1f} eV" if dis_win_min is not None else "") + "):")
    print(f"  {'dis_froz_max':>14}   max bands   status")
    print(f"  {'─'*42}")

    last_ok = None
    for e_test in [Ef + x for x in [-1, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10]]:
        max_n = max(count_bands(b, e_test) for b in blocks)
        ok = max_n <= nbndsub
        marker = "✓" if ok else "✗ EXCEEDS"
        print(f"  {e_test:>14.3f} eV   {max_n:>9}   {marker}")
        if ok:
            last_ok = e_test

    gap = None
    if last_ok is not None:
        for e_test in [last_ok + x for x in [0.5, 1, 2, 3, 5, 10]]:
            max_n = max(count_bands(b, e_test) for b in blocks)
            if max_n > nbndsub:
                gap = e_test - last_ok
                break

    suggested = last_ok - 0.5 if last_ok else None
    print(f"\n  Coarse Ef  = {Ef:.3f} eV")
    if suggested:
        gap_str = f"  (gap to band {nbndsub+1}: ~{gap:.1f} eV)" if gap else ""
        print(f"  Suggested dis_froz_max = {suggested:.2f} eV{gap_str}")
    return Ef, suggested, gap


def parse_nqirr(ph_outfile):
    """Count number of irreducible q-points from ph.x output."""
    content = Path(ph_outfile).read_text()
    matches = re.findall(r'Calculation of q =', content)
    return len(matches)
