"""Independent Eliashberg-theory quantities computed from an α²F(ω) spectrum.

All energies in meV, temperatures in K. These are used to cross-check EPW's own
reported λ / ω_log / Tc, and to give a full isotropic Eliashberg Tc on the a2f path
(where EPW only reports Allen-Dynes).
"""

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

K_PER_MEV = 11.604518


@dataclass
class A2F:
    omega: np.ndarray   # meV, strictly positive
    a2f: np.ndarray     # dimensionless

    @classmethod
    def from_epw(cls, path, column=1):
        """Read an EPW `<prefix>.a2f*` file. Column 1 is the smallest phonon smearing."""
        rows = []
        for line in Path(path).read_text().splitlines():
            parts = line.split()
            if not parts or not _is_float(parts[0]):
                if rows:
                    break
                continue
            rows.append([float(x) for x in parts])
        data = np.array(rows)
        mask = data[:, 0] > 0
        return cls(omega=data[mask, 0], a2f=data[mask, column])

    def moments(self):
        w, f = self.omega, self.a2f
        lam = 2.0 * np.trapezoid(f / w, w)
        wlog = math.exp(2.0 / lam * np.trapezoid(f * np.log(w) / w, w))
        w2 = math.sqrt(2.0 / lam * np.trapezoid(f * w, w))
        return lam, wlog, w2

    def lambda_matsubara(self, nu):
        """λ(iν) = 2∫ α²F(ω) ω / (ω² + ν²) dω for an array of bosonic frequencies ν (meV)."""
        w, f = self.omega, self.a2f
        nu = np.atleast_1d(nu)[:, None]
        return 2.0 * np.trapezoid(f * w / (w ** 2 + nu ** 2), w, axis=1)

    def omega_max(self, frac=1e-3):
        nz = np.nonzero(self.a2f > frac * self.a2f.max())[0]
        return float(self.omega[nz[-1]])


def _is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def allen_dynes(lam, wlog_mev, muc, w2_mev=None):
    """Allen-Dynes Tc (K). With w2_mev, includes the f1·f2 strong-coupling corrections."""
    denom = lam - muc * (1.0 + 0.62 * lam)
    if denom <= 0:
        return 0.0
    tc = wlog_mev * K_PER_MEV / 1.2 * math.exp(-1.04 * (1.0 + lam) / denom)
    if w2_mev is None:
        return tc
    lam1 = 2.46 * (1.0 + 3.8 * muc)
    lam2 = 1.82 * (1.0 + 6.3 * muc) * (w2_mev / wlog_mev)
    f1 = (1.0 + (lam / lam1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1.0 + (w2_mev / wlog_mev - 1.0) * lam ** 2 / (lam ** 2 + lam2 ** 2)
    return f1 * f2 * tc


def _gap_kernel_eigenvalue(spec, T, muc, wc_mev):
    """Largest eigenvalue of the linearized isotropic Eliashberg kernel at temperature T.

    Uses gap evenness to fold onto positive Matsubara frequencies ω_n = πT(2n+1):
      Z_n Δ_n = πT Σ_{m≥0} [λ(n−m) + λ(n+m+1) − 2μ*] Δ_m / ω_m
      Z_n     = 1 + (πT/ω_n) Σ_{m≥0} [λ(n−m) − λ(n+m+1)]
    with the sum cut at ω_m ≤ ω_c. Tc is where the eigenvalue crosses 1.
    """
    piT = math.pi * T / K_PER_MEV
    n_max = max(8, int(math.ceil(wc_mev / (2.0 * piT))))
    n = np.arange(n_max)
    wn = piT * (2 * n + 1)
    lam_k = spec.lambda_matsubara(2.0 * piT * np.arange(2 * n_max + 1))
    diff = np.abs(n[:, None] - n[None, :])
    summ = n[:, None] + n[None, :] + 1
    lam_minus, lam_plus = lam_k[diff], lam_k[summ]
    Z = 1.0 + piT / wn * (lam_minus - lam_plus).sum(axis=1)
    K = piT * (lam_minus + lam_plus - 2.0 * muc) / wn[None, :] / Z[:, None]
    # D K D⁻¹ with D = sqrt(Z/ω) is symmetric, so eigvalsh applies.
    D = np.sqrt(Z / wn)
    S = K * D[:, None] / D[None, :]
    return float(np.linalg.eigvalsh(0.5 * (S + S.T))[-1])


def eliashberg_tc(spec, muc, wc_mev, t_lo=0.5, t_hi=None, tol=0.05):
    """Isotropic Eliashberg Tc (K) by bisection on the kernel eigenvalue. Returns 0 if Tc < t_lo."""
    lam, wlog, w2 = spec.moments()
    t_ad = allen_dynes(lam, wlog, muc, w2)
    if t_hi is None:
        t_hi = max(5.0, 2.0 * t_ad + 10.0)
    # Low T means many Matsubara points; start the bracket near the AD estimate when it is valid.
    guess = 0.5 * t_ad
    if guess > t_lo and _gap_kernel_eigenvalue(spec, guess, muc, wc_mev) > 1.0:
        t_lo = guess
    elif _gap_kernel_eigenvalue(spec, t_lo, muc, wc_mev) < 1.0:
        return 0.0
    while _gap_kernel_eigenvalue(spec, t_hi, muc, wc_mev) > 1.0:
        t_hi *= 2.0
    while t_hi - t_lo > tol:
        mid = 0.5 * (t_lo + t_hi)
        if _gap_kernel_eigenvalue(spec, mid, muc, wc_mev) > 1.0:
            t_lo = mid
        else:
            t_hi = mid
    return 0.5 * (t_lo + t_hi)


def summarize(spec, muc, wc_mev):
    lam, wlog, w2 = spec.moments()
    return {
        "lambda": lam,
        "omega_log_meV": wlog,
        "omega2_meV": w2,
        "omega_max_meV": spec.omega_max(),
        "tc_allen_dynes_K": allen_dynes(lam, wlog, muc),
        "tc_allen_dynes_f1f2_K": allen_dynes(lam, wlog, muc, w2),
        "tc_eliashberg_iso_K": eliashberg_tc(spec, muc, wc_mev),
    }
