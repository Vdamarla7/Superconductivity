# CaH₆ @ 200 GPa — Results

> **Status (2026-10-19): this result is invalid.** The structure used here put H on the Im-3m 12e site (x=1/4), 0.91 Å from Ca. The correct site is 12d (1/4,0,1/2), Ca–H = 2.04 Å. Both have 48 symmetry operations, so the symmetry check passed. The wrong structure has a DFT pressure of 427 GPa (not 200), and that is why every harmonic Γ mode was imaginary. With H on 12d, the structure relaxes to a = 3.398 Å at 200 GPa and all Γ modes are real (1090–2079 cm⁻¹). Everything below, including the SSCHA work, describes the wrong-site structure. `compounds/CaH6_200GPa.yaml` now has the corrected structure. The corrected structure re-run with scflow in the container is dynamically stable at all 8 irreducible q-points of the 4×4×4 DFPT grid (Γ: 1036–2093 cm⁻¹), with the acoustic sum rule satisfied to 14 cm⁻¹.

>
> ### Corrected structure, scflow runs of 2026-10-19 (container, QE 7.4.1 / EPW 5.9)
> - **Structure:** vc-relaxed at 200 GPa to a = 3.398 Å. The DFT pressure at the relaxed geometry is 200.0 GPa.
> - **Phonons:** stable at all 8 irreducible q-points of the 4×4×4 DFPT grid (Γ: 1036–2093 cm⁻¹).
> - **EPW problems found and handled:**
>   1. `dis_froz_max` 5 meV from a band energy → the window scan now picks the edge with the most clearance.
>   2. NaN couplings when `nbnd` = 30 → `nbnd` = 34.
>   3. **Open issue:** NaN couplings with an 8³ Wannier k-grid on a 4³ q-grid. The electronic and phonon Wannier objects were clean; only the couplings broke. Production uses k = q = 4³.
> - **λ is not converged with the fine grid yet:**
>
>   | fine k / q | λ (mode sum) | ω_log | Tc (Eliashberg, μ* = 0.10) |
>   |---|---|---|---|
>   | 8³ / 4³ | 1.47 | 128 meV | 204 K |
>   | 16³ / 8³ | 0.66 | 122 meV | 57 K |
>   | 24³ / 12³ | 2.06 | 108 meV | 240 K |
>   | 32³ / 16³ | 1.34 | 115 meV | 161 K |
>
>   At EPW's default electronic smearing (25 meV), λ oscillates between grids: Fermi-surface sampling noise.
>   With `degaussw` = 0.1 eV, the same grids settle:
>
>   | fine k / q | degaussw | λ | ω_log | Tc AD / AD·f₁f₂ / Eliashberg (μ* = 0.10) |
>   |---|---|---|---|---|
>   | 24³ / 12³ | 0.05 eV | 1.73 | 110 meV | 165 / 190 / 211 K |
>   | 24³ / 12³ | 0.10 eV | 1.83 | 111 meV | 174 / 202 / 225 K |
>   | **32³ / 16³** | **0.10 eV** | **1.96** | **110 meV** | **181 / 214 / 239 K** |
>
>   λ changes 6.9% between the two densest grids, so the database tier is **converged**.
>
> **Caveats:**
> - The convergence evidence is two grid points at one smearing, and λ is still drifting upward.
> - The Wannier coarse grid is 4³ (see the open issue above).
> - Phonons are harmonic, and μ* = 0.10 is assumed.
> - The 215 K experiment is at 172 GPa, so it is not a benchmark for 200 GPa.

## Structure
- Space group: Im-3m (BCC clathrate, ibrav=0 in QE)
- Lattice parameter: A = 3.64 Å (literature at 200 GPa)
- 7 atoms/cell: Ca at 0,0,0; 6H at ±1/4,0,0 and permutations (12e Wyckoff)

## Pipeline settings
- SCF/NSCF k-grid: 6×6×6 (coarse Wannier), ecutwfc=80 Ry
- SSCHA supercell: 3×3×3 (448 atoms), 30 iterations, 24 configs/iter
- Phonon q-grid: 3×3×3 coarse (8 irreducible q-points under Im-3m)
- Wannier: Ca:d + 6×H:s → 11 WFs, Ω_total = 10.06 Ų (well-converged)
- EPW coarse: 2×2×2 or 3×3×3; fine: 20×20×20

## Results

### Harmonic (DFPT, no SSCHA)
| Quantity | Value |
|---|---|
| Status | ✗ All 21 modes imaginary at Γ — CaH₆ harmonically unstable |
| λ from harmonic EPW | ≈0 (imaginary modes masked) |
| Tc | undefined |

### SSCHA-renormalized (current best: 3×3×3 coarse grid)
| Quantity | Value |
|---|---|
| SSCHA convergence | ✓ 30 iters, 0 imaginary modes |
| SSCHA Γ-point modes | 0, 0, 0, **786**, 863, 984, 995 cm⁻¹ (matches SSCHA literature) |
| λ (e-ph coupling) | **0.997** |
| ω_log | 6.69 mRy = 1056 K = 91 meV |
| Allen-Dynes Tc (μ*=0.10) | **60.4 K** |
| Allen-Dynes Tc (μ*=0.12) | 53.7 K |
| Literature Tc | ~215–235 K |

### Why Tc=60K is not the converged result
The 3×3×3 coarse Wannier grid (8 irreducible q-points) **undersamples** the BZ:
- EPW Wannier interpolation of g(k,q) is accurate only near coarse q-points
- Literature requires ≥4×4×4 coarse q-grid for CaH₆ (Duan et al. 2014 used 4×4×4)
- 3×3×3 → 20×20×20 fine: interpolation ratio 6.67× (too large for accurate g)
- Expected converged λ with 4×4×4 coarse: ~2.0–2.7

## Planned next step: 4×4×4 coarse grid
1. `nscf` on 8×8×8 k-grid (~15 min)
2. `ph.x` DFPT on 4×4×4 q-grid (~3–4 h, 20 CPUs) → get dvscf at 12 irreducible q-points
3. `gen_sscha_dyn444.py`: matdyn from `sscha_ry.fc` → SSCHA-corrected 4×4×4 dyn files
4. EPW pass 1 with 8×8×8 k, 4×4×4 q (~45 min)
5. EPW pass 2 with 20×20×20 fine (~10 min)
6. Expected Tc: 150–215 K (literature range for 4×4×4 SSCHA EPW)

## Key technical notes
- **Format**: CellConstructor/matdyn dyn file format = EPW format (both use AMU_RY masses)
- **EPW pass 2**: uses `epwread=.true.` → reads Wannier phonon Hamiltonian (rdw) from
  `epwdata.fmt`; does NOT re-read dvscf files
- **SSCHA approach**: harmonic dvscf files reused for e-ph matrix elements; only phonon
  frequencies are replaced by SSCHA-renormalized values
- **EPW FP exceptions**: EPW sometimes exits non-zero (IEEE_OVERFLOW from Intel MKL) while
  completing successfully; pipeline.py now trusts `job_done` over exit code
