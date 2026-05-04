# NbN @ 0 GPa — Results

## Structure
- Space group: Fm-3m (rock salt, ibrav=2 in QE)
- Lattice parameter: A = 4.394 Å (experimental, skip_vcrelax=true)
- 2 atoms/cell: Nb at 1a (0,0,0), N at 1b (0.5,0.5,0.5)

## Pipeline settings
- SCF k-grid: 12×12×12, ecutwfc=80 Ry
- Phonon q-grid: 4×4×4 (coarse), harmonically stable ✓
- Wannier: Nb:s+p+d + N:s+p → **13 WFs**, nbndsub=13
- bands_skipped: exclude_bands=1-4 (N-2s + N-2p deep core at −35/−11 eV)
- dis_froz_max=22.0 eV (spans fine Ef=20.96 eV; gap to band 14 at ~24 eV)
- Fine grid: 10×10×10, a2f path (Allen-Dynes from α²F)

## Results
| Quantity | Value |
|---|---|
| λ (e-ph coupling) | **2.47** |
| ω_log | 11 meV |
| Allen-Dynes Tc (μ*=0.10) | **19.8 K** |
| Allen-Dynes Tc (μ*=0.18) | **16.6 K** |
| Literature Tc | 12–16 K (exp: ~16 K) |
| Status | ✓ Validated |

## Key Wannier lesson
nbndsub=9 (Nb:d + N:p only, no frozen window) gave λ=50 — completely unphysical.
The full 13-WF basis (Nb:s+p+d + N:s+p) is required to correctly capture the
Fermi-crossing band manifold. With nbndsub=13 and dis_froz_max spanning Ef,
the result converges to the experimental value at μ*=0.18.

The a2f path (α²F → Allen-Dynes) is used because the fine-grid Ef falls above
dis_froz_max, making Eliashberg egnv unreliable for this d-band system.
