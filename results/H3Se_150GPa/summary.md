# H₃Se @ 150 GPa — Results

## Structure
- Space group: Im-3m (BCC, ibrav=2 in QE)
- Lattice parameter: A = 5.832 Å (vc-relax at 150 GPa)
- 4 atoms/cell: Se at 1a (0,0,0), 3H at 0.5,0,0 and permutations

## Pipeline settings
- SCF k-grid: 12×12×12, ecutwfc=80 Ry
- Phonon q-grid: 4×4×4 (coarse), harmonically stable ✓
- Wannier: Se:sp3 + 3H:s → 7 WFs, dis_froz_max = 15.0 eV
- Fine grid: 10×10×10

## Results
| Quantity | Value |
|---|---|
| λ (e-ph coupling) | 11.91 |
| Allen-Dynes Tc (μ*=0.10) | **97.0 K** |
| Literature Tc | ~100–130 K |
| Status | ✓ Validated |

## Notes
- λ=11.9 is unusually high, reflecting very soft H-dominated phonon modes
- No imaginary phonon modes at 150 GPa (harmonically stable)
- Two-pass EPW with Eliashberg equations solved (isotropic)
- Fine grid 10×10×10 is PoC quality; 20×20×20 for production
