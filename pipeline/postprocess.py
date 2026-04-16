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
    """Extract Allen-Dynes Tc and omega_log from EPW pass 2 output."""
    content = Path(epw_outfile).read_text()
    tc_m = re.search(r'Estimated Allen-Dynes Tc\s*=\s*([\d.]+)\s*K', content)
    wlog_m = re.search(r'Estimated w_log in Allen-Dynes Tc\s*=\s*([\d.]+)\s*meV', content)
    tc  = float(tc_m.group(1))  if tc_m  else None
    wlog = float(wlog_m.group(1)) if wlog_m else None
    return tc, wlog


def parse_nqirr(ph_outfile):
    """Count number of irreducible q-points from ph.x output."""
    content = Path(ph_outfile).read_text()
    matches = re.findall(r'Calculation of q =', content)
    return len(matches)
