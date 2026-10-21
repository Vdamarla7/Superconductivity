# Nb3Sn @ 0 GPa — Results

> **Status (2026-10-19): not completed.**

## Structure
- Space group: Pm-3n (A15, ibrav=0 in QE)
- Lattice parameter: A = 5.289 Å (experimental, skip_vcrelax=true)
- 8 atoms/cell: 2×Sn at 2a Wyckoff, 6×Nb at 6c Wyckoff
- 48 sym ops ✓

## Pipeline settings
- SCF k-grid: 12×12×12, ecutwfc=60 Ry
- Phonon q-grid: 2×2×2 (coarse), harmonically stable ✓
- Wannier: Nb:s+d (6/Nb × 6) + Sn:s+p (4/Sn × 2) → **44 WFs**, nbndsub=44
- bands_skipped: exclude_bands=1-34 (Nb-4s×6, Sn-4d×2, Nb-4p×6 deep core)
- EPW coarse k-grid: 6×6×6 (216 k-pts for 44 WFs; must be multiple of q-grid)
- Fine grid: 10×10×10, a2f path (Allen-Dynes from α²F)
- dis_froz_max: pending NSCF window scan (placeholder 15.0 eV in YAML)

## Status: IN PROGRESS — phonons running

Phonon calculation in progress (2×2×2 q-grid, 8 atoms, ~4 irred. q-points).
Parallelised with 6 MPI k-point pools (OMP_NUM_THREADS=1).

Progress as of last check: Γ done ✓, X-point in progress.

Results will be added here once EPW completes. Expected Tc: 16–20 K (exp: 18.3 K).
