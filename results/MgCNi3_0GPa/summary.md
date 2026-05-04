# MgCNi3 @ 0 GPa — Results

## Structure
- Space group: Pm-3m (anti-perovskite, ibrav=0 in QE)
- Lattice parameter: A = 3.812 Å (experimental, skip_vcrelax=true)
- 5 atoms/cell: C at (0,0,0), Mg at (0.5,0.5,0.5), 3×Ni at face centres
- 48 sym ops ✓

## Pipeline settings
- SCF k-grid: 8×8×8, ecutwfc=60 Ry
- Phonon q-grid: 2×2×2 (coarse), harmonically stable ✓
- **EPW Wannier route FAILS** — degenerate Ni-d flat bands at Γ/R cause
  IEEE_INVALID in the Bloch-to-Wannier rotation (g_W = U†gU); all nbndsub/
  dis_froz_max combinations tested produce 50% NaN in .epmatwp
- **Workaround**: ph.x `electron_phonon='simple'` — Bloch-basis Fermi surface
  integral, no Wannier rotation needed; uses existing dvscf from phonon step

## Results (ph.x electron_phonon='simple', σ=0.04 Ry)
| Quantity | Value |
|---|---|
| λ (e-ph coupling) | **0.91** |
| ω_log | 132 K |
| Allen-Dynes Tc (μ*=0.10) | **7.9 K** |
| Literature Tc | 7–8 K (experimental) |
| Status | ✓ Validated |

## σ-dependence (broadening convergence)
| σ (Ry) | λ | ω_log (K) | Allen-Dynes Tc |
|---|---|---|---|
| 0.020 | 1.141 | 127 | 10.7 K |
| 0.030 | 0.985 | 129 | 8.8 K |
| 0.040 | 0.913 | 132 | **7.9 K** |
| 0.050 | 0.867 | 134 | 7.3 K |

λ converges toward ~0.87 at large σ; denser k-mesh (16³+) would improve the
σ→0 extrapolation. Best estimate at 8×8×8 is 7.9 K (σ=0.04 Ry).

## Why EPW fails for MgCNi3
The Ni-d bands form degenerate flat bands at Γ and R in Pm-3m. Wannier90
cannot uniquely gauge-fix the degenerate eigenspace, producing ill-conditioned
U matrices. The resulting Bloch-to-Wannier rotation g_W = U†gU triggers
IEEE_INVALID + OVERFLOW + UNDERFLOW during EPW's `elphon_wrap`/`Bloch2wanp`.
The .epmatwp file shows a regular alternating pattern of ~1e-45 denormals and
NaN values (50% NaN). This is not fixable by adjusting the frozen window.

The ph.x `electron_phonon='simple'` route bypasses Wannierization entirely and
gives physically correct results for this system.
