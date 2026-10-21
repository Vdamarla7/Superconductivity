# H₃S @ 200 GPa — Results

> **Status (2026-10-19): needs re-run under scflow.** The outputs of this run were not kept, so the validation gate cannot re-check it. The configured fsthick (0.6 eV) is also below 3 × ω_max for H₃S (~0.63 eV); the gate will warn about this.

## Structure
- Space group: Im-3m (BCC, ibrav=3 in QE)
- Lattice parameter: A = 5.739 Å (vc-relax at 200 GPa)
- 4 atoms/cell: S at 1a (0,0,0), 3H at 0.5,0,0 and permutations

## Pipeline settings
- SCF k-grid: 12×12×12, ecutwfc=80 Ry
- Phonon q-grid: 4×4×4 (coarse), harmonically stable ✓
- Wannier: S:sp3 + 3H:s → 7 WFs, dis_froz_max = Ef + 3 eV
- Fine grid: 20×20×20

## Results
| Quantity | Value |
|---|---|
| Allen-Dynes Tc (μ*=0.10) | **216 K** |
| Literature Tc | ~200–220 K (exp: 203 K) |
| Status | ✓ Validated |

## Notes
- H₃S is the benchmark compound for high-pressure hydride superconductivity
- Experimental Tc = 203 K (Drozdov et al., Nature 2015)
- Our result at 200 GPa recovers the literature value within expected DFT/k-grid error
- No imaginary phonon modes at 200 GPa (harmonically stable)
- This is the primary pipeline validation run
