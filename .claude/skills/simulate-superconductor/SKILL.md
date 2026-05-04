---
name: simulate-superconductor
description: Run the full EPW superconductivity pipeline for a compound — from crystal structure to Allen-Dynes or Eliashberg Tc. Handles structure setup, pseudopotential selection, Wannier window tuning, and pipeline execution.
argument-hint: "[chemical-formula] [pressure_GPa] OR [path/to/compound.yaml]"
---

You are helping the user simulate a superconductor using the QE + EPW pipeline in this repository.

The argument is either:
- A compound formula + pressure: e.g. `H3Se 150` or `NbN 0`
- A path to an existing YAML config: e.g. `compounds/NbN_0GPa.yaml`

---

## OVERVIEW

The pipeline lives in `pipeline/` and consists of:
- `pipeline.py` — orchestrator (6 steps: vcrelax → scf → phonons → nscf → epw1 → epw2)
- `generators.py` — QE/EPW input file generators
- `postprocess.py` — output parsers (Ef, Tc, sym ops, window scan)
- `kpoints.py` — explicit k-point grid generator for NSCF/EPW
- `sscha_relax.py` — anharmonic phonon renormalization (for clathrate hydrides)

Run from the repo root:
```bash
python pipeline/pipeline.py compounds/MyCompound.yaml [--start-from STEP] [--dry-run] [--reset-epw]
```

Output goes to `pipeline/runs/<prefix>/` with subdirectories `1_vcrelax/`, `2_phonons/`, `3_epw/`.

**IMPORTANT**: Each compound needs a unique `prefix` in its YAML — two compounds sharing a prefix will use each other's cached outputs.

### Key flags:
- `--reset-epw`: deletes stale `.epmatwp`, `.epb*`, `.ephmat`, `epw_pass1/2.out` — use whenever you change Wannier windows or nbndsub
- `--start-from epw1`: re-Wannierize without rerunning phonons (requires `--reset-epw` if windows changed)
- `--dry-run`: write input files only, don't execute

### Automatic window scan:
When `--start-from epw1` runs, the pipeline scans NSCF eigenvalues and prints a `dis_froz_max` recommendation table for your nbndsub. Always check this before finalising the YAML.

---

## STEP 1 — DETERMINE STRUCTURE AND ibrav

**The single most important decision: ibrav=3 vs ibrav=0**

### Use ibrav=3 ONLY when:
- Space group is Im-3m (BCC)
- AND all atomic positions use only 0 and 0.5 fractional coordinates

This covers: H3S, H3Se, H3Te, Nb (elemental).

### Use ibrav=0 for EVERYTHING ELSE:
- CaH6, LaH10 (H at 12e Wyckoff: ±0.25,0,0)
- Rock-salt FCC (NbN, PdH): use ibrav=2 (positions ARE 0 and 0.5)
- Pm-3m simple cubic (MgCNi3): use ibrav=0 with `cell_parameters_angstrom: [[a,0,0],[0,a,0],[0,0,a]]`
- Pm-3n A15 (Nb3Sn, V3Si): 6c Wyckoff has 0.25/0.75 coords → must use ibrav=0

With ibrav=0, QE does full symmetry detection from coordinates — always finds the correct sym ops.

### Verify with spglib before running:
```python
import spglib
ds = spglib.get_symmetry_dataset((cell, positions, numbers), symprec=1e-3)
print(ds['international'], len(ds['rotations']))  # e.g. Im-3m, 48
```

### Always verify sym ops after SCF:
```bash
grep "Sym\. Ops\." pipeline/runs/<prefix>/2_phonons/scf.out
```
Expected: 48 for Im-3m, Fm-3m, Pm-3m, Pm-3n. Wrong sym ops = wrong phonons = wrong Tc.

---

## STEP 2 — CHOOSE PSEUDOPOTENTIALS

Always use USPP (rrkjus) for production. Set paths in your compound YAML:
```yaml
paths:
  qe_bin:     "/path/to/qe-7.x/bin"
  pseudo_dir: "/path/to/pslibrary/pbe/PSEUDOPOTENTIALS/"
```

Pseudopotential naming convention:
- `{Element}.pbe-n-rrkjus_psl.1.0.0.UPF` — no semicore (H, C, N, O, light elements)
- `{Element}.pbe-dn-rrkjus_psl.1.0.0.UPF` — semicore d (Sn, Pb, Cd, ...)
- `{Element}.pbe-spn-rrkjus_psl.1.0.0.UPF` — semicore s+p (Nb, Ni, Pd, Mg, Ca, ...)

For heavy elements (Ir, Pt, Os, Rh, Ru, Re): use `rel-pbe/PSEUDOPOTENTIALS/` (relativistic).

---

## STEP 3 — SET WANNIER PARAMETERS

This is the most compound-specific and failure-prone step.

### nbndsub and projections:
| Compound type | Projections | nbndsub |
|---|---|---|
| H₃X family (H3S, H3Se, H3Te) | S/Se/Te:sp3 + nH:s | 4 + n |
| Ca/Mg clathrate hydrides | M:d + nH:s | 5 + n |
| Elemental d-metal (Nb, V) | M:d | 5 |
| Rock-salt d-metal nitride (NbN, VN) | M:s+p+d + N:s+p | 13 |
| Rock-salt d-metal hydride (PdH) | M:d + H:s | 6 |
| Ternary anti-perovskite (MgCNi3) | EPW fails — use ph.x workaround (see below) | — |
| A15 ternary (Nb3Sn, V3Si, Nb3Al) | (M:s+M:d)×6 + (X:s+X:p)×2 | 44 |

Count projections carefully: `"Ni:d"` at 3 Ni sites → 15 WFs (5×3). nbndsub must equal the total WF count.

Use explicit orbitals (`"Se:s"` + `"Se:p"`) instead of `"Se:sp3"` if sp3 gives NaN coupling.

### Deep core bands with spn/dn pseudopotentials:
Any compound using `pbe-spn-rrkjus` or `pbe-dn-rrkjus` has semicore bands in the valence — deep low-index bands well below −10 eV that must be excluded from the Wannier window:

| Pseudo type | Semicore | Bands per atom | Electrons per atom |
|---|---|---|---|
| spn — 4d metals (Nb, Mo, Pd) | 4s+4p | 4 | 8 |
| spn — 3d metals (Ni, Co, Fe) | 3s+3p | 4 | 8 |
| spn — alkaline earth (Mg, Ca) | 2s+2p or 3s+3p | 4 | 8 |
| dn — (Sn, Pb, Cd) | 4d | 5 | 10 |
| n — (H, C, N, O) | none | 0 | 0 |

**Recipe for bands_skipped:**
1. List all pseudos; identify semicore bands per atom type × atom count.
2. Order by binding energy (deepest first) — typically: metal-4s, then dn-4d, then metal-4p.
3. `bands_skipped: "exclude_bands = 1-K"` where K = total semicore band count.
4. Verify: (total valence e⁻) − (excluded e⁻) = 2 × (active filled bands), and active filled bands ≤ nbndsub.

**Example — NbN** (Nb.spn × 1, N.n × 1):
- Nb-4s (1 band) + Nb-4p (3 bands) = 4 bands → `exclude_bands = 1-4`
- Removes 8e. Active: 13e(Nb) + 5e(N) − 8e = 10e → 5 active filled bands ✓

**Example — Nb3Sn** (Nb.spn × 6, Sn.dn × 2):
- Nb-4s×6 (bands 1-6) → Sn-4d×2 (bands 7-16) → Nb-4p×6 (bands 17-34) = 34 bands
- `exclude_bands = 1-34`, removes 68e. Active: 38e → 19 filled active bands (nbndsub=44).

### Wannier windows — CRITICAL RULES:

**Rule 1**: Run SCF first; get coarse Ef from scf.out.  
**Rule 2**: After NSCF, get fine-grid Ef from nscf.out (often 0.3–0.8 eV above coarse Ef).  
**Rule 3**: dis_froz_max must contain ≤ nbndsub bands at EVERY k-point.

```
dis_win_min:  Ef - 15 eV   (well below all relevant bands)
dis_win_max:  Ef + 20 eV   (must contain ≥ nbndsub bands at every k-point)
dis_froz_min: same as win_min
dis_froz_max: Ef + 2–4 eV  (MUST be verified from NSCF scan — pipeline prints this automatically)
```

**The pipeline prints a window scan table automatically at epw1.** Or run manually:
```python
import re
content = open("pipeline/runs/<prefix>/3_epw/nscf.out").read()
Ef = float(re.search(r'the Fermi energy is\s+([\d.+-]+)', content).group(1))
blocks = re.findall(r'bands \(ev\):\s*([\d\s.+-]+?)(?=\n\s*\n|\n\s*k =)', content, re.DOTALL)
for e_test in [Ef+1, Ef+2, Ef+3, Ef+4, Ef+5]:
    max_n = max(sum(1 for e in map(float, b.split()) if e < e_test) for b in blocks)
    print(f"dis_froz_max={e_test:.2f}: max {max_n} bands (need ≤ {nbndsub})")
```

**Pass1 and pass2 MUST use identical Wannier windows.** If you change windows, use `--reset-epw`.

### d-band systems: frozen window and the a2f path

For d-band systems where too many bands enter the frozen window before it spans Ef, the fix is to **expand nbndsub** until all Fermi-crossing bands fit inside the frozen window. Example: NbN needs nbndsub=13 (not 9) because bands 12-13 both cross Ef and require the full Nb+N manifold.

Set `use_a2f: true` in the YAML for d-band compounds. This selects the α²F → Allen-Dynes route in pass2, which works even when the fine-grid Ef sits slightly above dis_froz_max (making full Eliashberg egnv unreliable).

**If nbndsub cannot be extended** (degenerate d-bands cause NaN regardless — see MgCNi3): use `ph.x electron_phonon='simple'` instead of EPW.

### Validation in epw_pass2.out:
- `lambda : X.XX` — real, not NaN, physically reasonable (< 5 for most materials)
- `N bands within the Fermi window` > 0
- Fine-grid Ef shift < 0.5 eV vs coarse (larger → poor Wannier interpolation)

---

## STEP 4 — EPW PARAMETERS

```yaml
epw:
  fine_kgrid:  [10, 10, 10]   # PoC; [20, 20, 20] for production
  fsthick:     1.0             # eV — Fermi surface thickness; start at 1.0
  wscut:       X.X             # Phonon frequency cutoff in eV (~8–12× max phonon freq):
                               #   H₃X hydrides: 2.0–3.0 eV
                               #   Clathrates (CaH6, LaH10): 3.0–4.0 eV
                               #   d-metals (NbN, Nb, Nb3Sn): 0.3–0.5 eV
  use_a2f:     true/false      # true for d-band (α²F path); false for H₃X/clathrates
```

Temperature scan: bracket the expected Tc. If unknown, use `[50, 100, 150, 200, 250, 300]` K.

---

## STEP 5 — CREATE YAML AND RUN

1. Copy `compounds/template.yaml` → `compounds/{Formula}_{P}GPa.yaml`
2. Fill in structure, pseudopotentials, Wannier parameters, and QE paths
3. Dry-run: `python pipeline/pipeline.py compounds/X.yaml --dry-run`
4. Verify 48 sym ops after SCF before committing hours to phonons
5. After NSCF: check printed window scan, update `dis_froz_max`, then `--reset-epw --start-from epw1`
6. Full run: `python pipeline/pipeline.py compounds/X.yaml`

---

## ANHARMONIC PHONONS (SSCHA)

CaH₆, LaH₁₀, PdH, and other quantum anharmonic crystals show imaginary harmonic phonon modes. This is not a pipeline error — it is real physics (nuclear quantum effects required). The harmonic EPW result gives λ from real modes, but Tc is suppressed. SSCHA renormalizes the phonon spectrum and recovers the correct Tc.

```bash
python pipeline/sscha_relax.py compounds/CaH6_200GPa.yaml --nconfigs 50 --niter 20 --temp 300
```

After SSCHA: copy renormalized dyn files over harmonic ones, then:
```bash
python pipeline/pipeline.py compounds/CaH6_200GPa.yaml --start-from nscf
```

**Key parameters**: `--nconfigs 10` for a quick PoC (~20 min/iter); `50` for production.  
**Runtime**: ~17 hours on a modern laptop for CaH₆ (448-atom supercell, 50 configs × 10 iters).

---

## VALIDATED RESULTS

| Compound | Pressure | Tc (this pipeline) | Literature | Method |
|----------|----------|--------------------|------------|--------|
| H₃S  | 200 GPa | **216 K** | ~200–220 K (exp: 203 K) ✓ | Full Eliashberg, 20³ fine grid |
| H₃Se | 150 GPa | **97 K**  | ~100–130 K ✓ | Full Eliashberg, 10³ fine grid |
| CaH₆ | 200 GPa | 22.6 K ⚠ | ~220–235 K | Harmonic only — SSCHA required |
| Nb   | 0 GPa   | **11.4 K** (μ\*=0.10) | 9.3 K exp ✓ | a2f, 10³ fine grid |
| NbN  | 0 GPa   | **16.6 K** (μ\*=0.18) | 12–16 K ✓ | a2f, 10³ fine grid, nbndsub=13 |
| MgCNi3 | 0 GPa | **7.9 K** | 7–8 K exp ✓ | ph.x workaround (EPW NaN) |
| Nb₃Sn | 0 GPa  | in progress | 18.3 K exp | a2f, 44 WFs |
| PdH  | 0 GPa   | blocked ⚠ | ~9 K exp | Imaginary phonons — SSCHA needed |

---

## COMPOUND-SPECIFIC NOTES

### NbN — rock-salt nitride (nbndsub=13 required)

nbndsub=9 (Nb:d + N:p, no frozen window) gives λ=50 — completely unphysical. The full 13-WF basis (Nb:s+p+d + N:s+p) is required so the frozen window can span Ef. With nbndsub=13 and dis_froz_max=22.0 eV (> fine Ef=20.96 eV, gap to band 14 at ~24 eV), EPW gives λ=2.47, Tc=16.6 K (μ\*=0.18) ✓.

**Band structure** (after exclude_bands=1-4):
- Band 5: Nb-5s/N-2s at ~3 eV
- Bands 6-11: Nb-t2g + Nb-eg at ~17–19 eV
- Bands 12-13: crossing Ef at ~20.6 eV
- Band 14+: above ~24 eV (clean gap)

### A15 compounds (Nb3Sn, V3Si, Nb3Al)

Pm-3n (#223), 8 atoms/cell, ibrav=0 required. Standard atomic positions:
```
X at 2a: (0,0,0), (0.5,0.5,0.5)
M at 6c: (0.25,0,0.5), (0.75,0,0.5), (0.5,0.25,0), (0.5,0.75,0), (0,0.5,0.25), (0,0.5,0.75)
```

Use `epw_coarse_kgrid: [6,6,6]` (must be a multiple of `coarse_kgrid: [2,2,2]`) to give ~5:1 ratio of k-points to WFs for 44 WFs. dis_froz_max must be scanned from NSCF — run through nscf, check printed table, update YAML, then `--reset-epw --start-from epw1`. Runtime: ~12–16 hr (2×2×2 phonons, 6³ Wannier, 10³ fine).

### MgCNi3 — EPW NaN workaround (ph.x route)

EPW Wannierization fails for MgCNi3 regardless of nbndsub or dis_froz_max. Root cause: degenerate Ni-d flat bands at Γ/R in Pm-3m produce ill-conditioned Wannier rotation matrices → IEEE_INVALID in g_W = U†gU → 50% NaN in .epmatwp. Not fixable via window tuning.

**Fix**: use `ph.x electron_phonon='simple'` (no Wannier rotation needed):
```fortran
&INPUTPH
  prefix='MgCNi3', outdir='./tmp/', fildyn='MgCNi3.dyn', fildvscf='dvscf'
  ldisp=.true., trans=.false., nq1=2, nq2=2, nq3=2
  electron_phonon='simple', el_ph_sigma=0.005, el_ph_nsigma=10
/
```
Then integrate λ from `.dyn*.elph.*` output files. Note: `lambda.x` fails for metallic systems (inconsistent DOS across q-points) — compute λ, ω_log, and Tc manually. See `results/MgCNi3_0GPa/summary.md` for the σ-convergence series.

---

## PLATFORM NOTES

### Apple Silicon (ARM64) — QE 7.4.1 hexagonal crash
ph.x segfaults for hexagonal space groups (P6/mmm, P6₃/mmc) on Apple Silicon. Occurs before the irreps table is printed (`SIGSEGV`, no useful backtrace).

**Not affected**: all cubic structures (Im-3m, Fm-3m, Pm-3m, Pm-3n).  
**Workaround**: run hexagonal compounds on Linux/cloud.

### MPI parallelism for long phonon runs
ph.x supports MPI k-point pools via `-nk N`. ph.x checkpoints to `.phsave/` after each completed representation — safe to kill mid-run and restart with `recover = .true.` in ph.in:
```bash
OMP_NUM_THREADS=1 mpirun -np N ph.x -nk N -in ph.in >> ph.out
```
Use N = number of physical cores you want to dedicate (leave ≥1 free for the OS).

---

## KNOWN BUGS (fixed in current code)

| Symptom | Root cause | Status |
|---------|-----------|--------|
| NSCF "inconsistent k-grids" in EPW | `K_POINTS {automatic}` gives symmetry-reduced set | Fixed: explicit `K_POINTS {crystal}` N³ list |
| EPW `write_freq` error 108 | `tmp/{prefix}.ephmat/` missing before pass2 | Fixed: pipeline pre-creates directory |
| `A: null` crash in generator | Missing lattice constant in YAML | Fixed: generator raises clear error |
| EPW pass1 falsely reported FAILED | `job_done()` only matched `JOB DONE` | Fixed: also accepts `Total program execution` |

---

## QUICK DIAGNOSIS

| Symptom | Cause | Fix |
|---------|-------|-----|
| 4 sym ops instead of 48 | ibrav=3 with 0.25/0.75 positions | Switch to ibrav=0 |
| `efermig: cannot bracket Ef` | dis_froz_max below Ef | Run window scan; raise dis_froz_max above fine Ef |
| EPW NaN coupling (λ=NaN) | Frozen window wrong | Recompute dis_froz_max from NSCF eigenvalue scan |
| λ=50 or other unphysical value | nbndsub too small; distorted Fermi-surface WFs | Expand nbndsub to span full Fermi manifold (see NbN) |
| 50% NaN pattern in .epmatwp | Degenerate d-bands → IEEE_INVALID in Bloch-to-Wannier | Use ph.x `electron_phonon='simple'` (see MgCNi3) |
| EOF reading .epmatwp | Wannier windows changed between pass1 and pass2 | `--reset-epw --start-from epw1` |
| `egnv` = 0 bytes | dis_froz_max doesn't span fine Ef | Set `use_a2f: true` for Allen-Dynes from α²F |
| Fine-grid Ef shift > 0.5 eV | Poor Wannier interpolation | Change projections (e.g. sp3 → explicit s+p) |
| Eliashberg Broyden factorisation crash | All temperatures above Allen-Dynes Tc | Lower `temperatures_K` below Tc |
| ph.x SIGSEGV before irreps table | Hexagonal space group on Apple Silicon | Run on Linux; avoid hexagonal on Mac |
| ph.x slow convergence (100+ iters/rep) | alpha_mix too aggressive for d-metals | Set `alpha_mix(1)=0.1`, `nmix_ph=8`, `tr2_ph=1.0d-12` |
| Imaginary phonon modes | Quantum anharmonic crystal (CaH₆, PdH, ...) | Not an error — SSCHA required |
| CaH₆/LaH₁₀ Tc << literature | Imaginary modes contribute zero to λ | Run sscha_relax.py to renormalize phonon spectrum |
