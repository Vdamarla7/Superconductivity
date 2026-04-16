"""
K-point grid utilities for EPW pipeline.
"""


def gen_kpoints_crystal(n1, n2, n3):
    """
    Return K_POINTS {crystal} block with all n1*n2*n3 explicit k-points.
    EPW requires the full unshifted grid — symmetry-reduced K_POINTS {automatic}
    will cause 'inconsistent nscf and elph k-grids' in wannierize.
    """
    lines = [f"K_POINTS {{crystal}}", f"{n1 * n2 * n3}"]
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                lines.append(f"  {i/n1:.6f}  {j/n2:.6f}  {k/n3:.6f}  1.0")
    return "\n".join(lines)


def gen_kpoints_automatic(n1, n2, n3, s1=0, s2=0, s3=0):
    """Return K_POINTS {automatic} block (for SCF only, NOT for NSCF fed into EPW)."""
    return f"K_POINTS {{automatic}}\n{n1} {n2} {n3} {s1} {s2} {s3}"
