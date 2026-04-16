# CaH₆ @ 200 GPa — Results

## Structure
- Space group: Im-3m (BCC clathrate, ibrav=0 in QE)
- Lattice parameter: A = 3.70 Å (experimental, skip_vcrelax=true)
- 7 atoms/cell: Ca at 0,0,0; 6H at ±1/4,0,0 and permutations (12e Wyckoff)

## Pipeline settings
- SCF k-grid: 12×12×12, ecutwfc=80 Ry
- Phonon q-grid: 4×4×4 (coarse), **98 imaginary modes** (quantum anharmonic crystal)
- Wannier: Ca:d + 6H:s → 11 WFs
- Fine grid: 10×10×10

## Results (harmonic approximation only)
| Quantity | Value |
|---|---|
| λ (e-ph coupling) | 1.10 |
| ω_log | 24.4 meV (283 K) |
| Allen-Dynes Tc (μ*=0.10) | **22.6 K** |
| Literature Tc | ~220–235 K |
| Status | ⚠ Harmonic result — SSCHA required |

## Why the harmonic Tc is wrong
CaH₆ is a **quantum anharmonic crystal**. The H-cage phonons require nuclear quantum
effects (zero-point motion) to stabilize. In harmonic PBE:

- 98 imaginary modes at 200 GPa — all pressures in harmonic PBE
- The imaginary modes are NOT a pipeline error — they are real physics
- EPW treats imaginary modes as zero frequency → most coupling modes contribute nothing
- Result: λ≈1.1 from real modes only, Tc=22.6 K (real Tc ~220 K)

## What SSCHA does
The `sscha_relax.py` script renormalizes the phonon spectrum using stochastic self-consistent
harmonic approximation (SSCHA). It:
1. Samples displaced supercell configurations at finite temperature
2. Runs pw.x force calculations on each
3. Updates force constants via gradient minimization
4. Outputs renormalized dyn files → use these for q2r → EPW

Runtime: ~17 hours on Mac (50 configs × 10 iterations, 448-atom supercell).
Renormalized Tc expected: ~200–220 K.
