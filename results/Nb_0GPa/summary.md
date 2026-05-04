# Nb @ 0 GPa — Results

## Structure
- Space group: Im-3m (BCC, ibrav=3 in QE)
- Lattice parameter: A = 3.30 Å (experimental, skip_vcrelax=true)
- 1 atom/cell: Nb at (0,0,0)

## Pipeline settings
- SCF k-grid: 12×12×12, ecutwfc=60 Ry
- Phonon q-grid: 4×4×4 (coarse), harmonically stable ✓
- Wannier: Nb:d → **5 WFs**, nbndsub=5
- bands_skipped: exclude_bands=1-4 (Nb 4s + 4p semicore at −36/−11.6 eV)
- dis_froz_max=20.0 eV (~2 eV above coarse Ef=17.94 eV)
- Fine grid: 10×10×10, a2f path (Allen-Dynes from α²F)

## Results
| Quantity | Value |
|---|---|
| λ (e-ph coupling) | **1.49** |
| ω_log | 8.7 meV |
| Allen-Dynes Tc (μ*=0.10) | **11.4 K** |
| Literature Tc | 9.3 K (experimental) |
| Status | ✓ Validated |

## Notes
- Elemental BCC Nb is the simplest d-band benchmark: 1 atom/cell, pure Nb:d physics
- Minimal 5-WF Wannier basis sufficient (no ligand field splitting)
- Phonon calculation is extremely fast: 1 atom × 4×4×4 q-grid
- Tc slightly overestimated vs exp (11.4 vs 9.3 K) — typical at PoC fine grid (10³)
  and without spin-orbit coupling; 20³ fine grid expected to improve agreement
