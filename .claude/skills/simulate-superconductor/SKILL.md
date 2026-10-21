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

The tool is the `scflow` package (`pip install -e .`, or the container in `docker/`):
- `workflow.py` — orchestrator (6 steps: vcrelax → scf → phonons → nscf → epw1 → epw2)
- `generators.py` — QE/EPW input file generators
- `postprocess.py` — output parsers (Ef, Tc, sym ops, window scan)
- `validation.py` — the physics gate (structure, pressure, phonons, Wannier, EPW checks)
- `eliashberg.py` — independent λ / ω_log / Allen-Dynes / isotropic Eliashberg Tc from α²F
- `db.py` — SQLite results DB + confidence tiers
- `sscha_relax.py` — anharmonic phonon renormalization (only after the structure checks pass)

Run from the repo root:
```bash
scflow check compounds/MyCompound.yaml          # pre-flight, seconds, no QE needed
scflow run   compounds/MyCompound.yaml [--start-from STEP] [--dry-run] [--reset-epw]
scflow validate compounds/MyCompound.yaml       # re-run the physics gate on a run dir
```

Output goes to `runs/<prefix>/` with subdirectories `1_vcrelax/`, `2_phonons/`, `3_epw/`.
Machine paths come from `SCFLOW_QE_BIN` and `SCFLOW_PSEUDO_DIR`, not from the YAML.

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
- CaH6: Ca at 2a (0,0,0), **H at 12d (1/4,0,1/2)** in the conventional cell. NOT 12e (1/4,0,0):
  that puts H 0.91 Å from Ca, still has Im-3m/48 sym ops, and gives a 427 GPa structure with
  every Γ mode imaginary. This exact mistake cost weeks on this project.
- LaH10: H on 8c and 32f (x,x,x) — general x, so ibrav=0
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

### Symmetry alone does not validate a structure
Different Wyckoff sites of the same space group give identical sym-op counts. Always run
`scflow check` — it fails on interatomic distances < 0.5 × covalent-radius sum, and after SCF
`scflow validate` fails if the DFT pressure is far from the target. Set
`structure.expected_spacegroup` so a typo in the positions cannot silently change the space group.

### Always verify sym ops after SCF:
```bash
grep "Sym\. Ops\." runs/<prefix>/2_phonons/scf.out
```
Expected: 48 for Im-3m, Fm-3m, Pm-3m, Pm-3n. Wrong sym ops = wrong phonons = wrong Tc.

---

## STEP 2 — CHOOSE PSEUDOPOTENTIALS

Always use USPP (rrkjus) for production. Set paths once per machine (the container sets them):
```bash
export SCFLOW_QE_BIN=/path/to/qe-7.x/bin
export SCFLOW_PSEUDO_DIR=/path/to/pslibrary/pbe/PSEUDOPOTENTIALS
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
**Rule 4**: dis_froz_max must not sit on a band energy at any k-point. CaH₆ at 22.5 eV was 5 meV
from a degenerate pair, which made the disentanglement ill-conditioned; every element of `.epmatwp`
came out NaN. The scan reports each candidate's clearance and suggests the value with the most
(CaH₆: 21.2 eV, 223 meV). `scflow` refuses to start EPW pass 1 under 10 meV.

```
dis_win_min:  Ef - 15 eV   (well below all relevant bands)
dis_win_max:  Ef + 20 eV   (must contain ≥ nbndsub bands at every k-point)
dis_froz_min: same as win_min
dis_froz_max: Ef + 2–4 eV  (MUST be verified from NSCF scan — pipeline prints this automatically)
```

**The pipeline prints a window scan table automatically after NSCF.** Or run manually:
```python
import re
content = open("runs/<prefix>/3_epw/nscf.out").read()
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
  use_a2f:     true/false      # true for d-band (α²F path); scflow also solves Eliashberg from α²F
  fine_qgrid:  [8, 8, 8]       # optional (default = fine_kgrid); cost ∝ Nk × Nq
  degaussw:    0.1             # optional electron smearing (eV); EPW default 0.025
```

**Convergence of λ is the hard part.** With EPW's default `degaussw = 0.025` eV, CaH₆ gave
λ = 0.69, 1.99, 1.32 on fine grids k16/q8, k24/q12, k32/q16: sampling noise, not physics. Scan
`degaussw` together with the grid, and record each combination in the database. A smearing forms
its own convergence series, and only a series that settles within 10% reaches `converged`.

The NSCF is written with `verbosity='high'`. Without it, QE omits band energies for more than 100
k-points, and the window scan would have nothing to read.

Temperature scan: bracket the expected Tc. If unknown, use `[50, 100, 150, 200, 250, 300]` K.

---

## STEP 5 — CREATE YAML AND RUN

1. Copy `compounds/template.yaml` → `compounds/{Formula}_{P}GPa.yaml`
2. Fill in structure (+ `expected_spacegroup`), pseudopotentials, Wannier parameters, and a `reference:` block
3. Pre-flight: `scflow check compounds/X.yaml` — must have no ✗
4. Through SCF only: `scflow run compounds/X.yaml --stop-after scf`, then `scflow validate` —
   check sym ops AND `scf.pressure` before committing hours to phonons
5. After NSCF: check printed window scan, update `dis_froz_max`, then `scflow run X.yaml --reset-epw --start-from epw1`
6. Full run: `scflow run compounds/X.yaml` — records into `runs/scflow.db` with a confidence tier
7. Convergence: rerun with a denser `epw.fine_kgrid`; the DB marks a method `converged` when λ agrees within 10%

---

## ANHARMONIC PHONONS (SSCHA)

**Imaginary harmonic modes are a structure problem until proven otherwise.** Before reaching for
SSCHA, confirm that `structure.distance.*` and `scf.pressure` pass. Hundreds-of-cm⁻¹ imaginary
modes at Γ in a high-symmetry hydride usually mean wrong positions (see CaH₆ above), not anharmonicity.

Genuine cases exist: PdH is the textbook example (strong H anharmonicity; harmonic X-point
instability). For those, SSCHA renormalizes the phonon spectrum:

```bash
pip install -e '.[sscha]'
python -m scflow.sscha_relax compounds/PdH_0GPa.yaml --nconfigs 24 --niter 30 --temp 300 --supercell 2 2 2
```

After SSCHA: copy renormalized dyn files over harmonic ones, then:
```bash
scflow run compounds/PdH_0GPa.yaml --start-from nscf
```

**Key parameters**:
- `--nconfigs 24`: faster (~100 min/iter), gradient ~40% noisier, 30 iters recommended
- `--nconfigs 50`: production quality (~210 min/iter), 20 iters sufficient
- `--supercell 2 2 2`: required for CaH₆ (2×2×2 = 56-atom supercell, 7 atoms/cell × 8)

**Actual runtimes on RTX 5090 GPU** (56-atom CaH₆, 6 parallel workers):
- Per SCF config: ~25 min (not 7 min — the fast path was cached hits)
- 24 configs/iter: ~100 min/iter
- 50 configs/iter: ~210 min/iter

**Critical settings** (found while running CaH₆ — on the wrong-site structure, but the SSCHA mechanics carry over):
- `conv_thr = 1.0d-7` in SCF input (1e-8 is too tight for 56-atom cell — takes 3× longer, not needed)
- `kong_liu_ratio=0.0` — disables KL resampling; use 0.0 for CaH₆ (KL ratio is borderline and triggers infinite resampling loop at 0.3–0.5)
- `minim.max_ka = 2` then `minim.run()` — limits outer gradient steps to 2 per SSCHA iteration (calling `minimization_step()` directly fails: it requires `self.minimizer` initialized inside `run()`)
- `ForcePositiveDefinite()` + `Symmetrize()` on initial dyn — removes imaginary modes and restores ASR before first iteration
- `dyn_sscha.Symmetrize()` after every gradient step — SSCHA gradient breaks the Acoustic Sum Rule

**MACE acceleration**: `sscha_relax.py` supports `sscha.force_calculator: mace` in the YAML. The MACE
training scripts are not part of the released tool yet. Target force RMSE < 50 meV/Å before
trusting an MLIP for SSCHA, and revalidate against DFT every few iterations.

---

## SHARING RESULTS

`scflow` sends finished runs to the public database `Vdamarla7/scflow-results` by default, as a
PR, when GitHub credentials are available. Otherwise they queue in `runs/outbox/` for
`scflow share`. Tell the user this before running anything for them, and use `--no-share` if they
don't want a particular compound published (for example, unpublished work).

## RESULTS STATUS

Under the automated validation gate (see README for details), CaH₆ is the only compound in the
`converged` tier. Earlier hand-validated numbers are kept in `results/*/summary.md`
with notes on what the gate found:

| Compound | Status | Why |
|---|---|---|
| H₃S 200 GPa | needs re-run | earlier run's outputs not retained; cannot be re-validated |
| H₃Se 150 GPa | needs re-run | vc-relax never actually ran (old orchestrator bug); pressure unknown |
| CaH₆ 200 GPa | **converged** | corrected 12d structure, stable at 4×4×4. λ = 1.96, Eliashberg Tc = 239 K at μ* = 0.10 (k32/q16, degaussw 0.1). The old 12e result is invalid |
| Nb 0 GPa | preliminary | λ=1.84 vs ~1.0–1.3 expected: unconverged disentanglement, 0.8% imaginary fine-grid modes, single 10³ grid |
| NbN 0 GPa | invalid (harmonic) | −320 cm⁻¹ mode at X in the DFPT run; quoted Tc used μ*=0.18 fitted to experiment |
| MgCNi3 0 GPa | needs re-run | ph.x workaround outputs not retained |
| Nb₃Sn 0 GPa | not run | — |
| PdH 0 GPa | needs SSCHA | genuine anharmonic instability at X |

---

## COMPOUND-SPECIFIC NOTES

### NbN — rock-salt nitride (nbndsub=13 required)

nbndsub=9 (Nb:d + N:p, no frozen window) gives λ=50 — completely unphysical. The full 13-WF basis (Nb:s+p+d + N:s+p) is required so the frozen window can span Ef. With nbndsub=13 and dis_froz_max=22.0 eV (> fine Ef=20.96 eV, gap to band 14 at ~24 eV), EPW gives λ=2.47, Tc=19.8 K (μ\*=0.10). The often-quoted 16.6 K used μ\*=0.18, i.e. μ\* fitted to experiment — and the harmonic DFPT run has a −320 cm⁻¹ mode at X, so this λ is not trustworthy yet.

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
| SSCHA `run()` loops for 7+ hours (173k "Redoing the step") | `kong_liu_ratio=0.3–0.5` triggers infinite resampling for CaH₆ (KL score borderline) | Fixed: `kong_liu_ratio=0.0` |
| SSCHA `AttributeError: 'NoneType' has no attribute 'update_dyn'` | `minimization_step()` requires `self.minimizer` initialized by `run()` | Fixed: `minim.max_ka=2; minim.run()` instead of calling `minimization_step()` directly |
| `ASR violation / translational modes = 0` after SSCHA gradient | SSCHA gradient step breaks Acoustic Sum Rule | Fixed: `Symmetrize()` after ForcePositiveDefinite, checkpoint load, and each gradient step |
| SCF takes 20+ min per config for 56-atom cell | `conv_thr=1e-8` too tight | Fixed: `conv_thr=1.0d-7` |

---

## QUICK DIAGNOSIS

| Symptom | Cause | Fix |
|---------|-------|-----|
| 4 sym ops instead of 48 | ibrav=3 with 0.25/0.75 positions | Switch to ibrav=0 |
| `efermig: cannot bracket Ef` | dis_froz_max below Ef | Run window scan; raise dis_froz_max above fine Ef |
| EPW NaN coupling (λ=NaN) | Frozen-window edge on a band energy, or window wrong | Use the scan's max-clearance dis_froz_max, then `--reset-epw --start-from epw1` |
| NaN couplings only at some coarse q (all modes) | `nbnd` ends inside a degenerate multiplet, which breaks EPW's Bloch-space rotation (CaH₆: nbnd 30 → NaN, 34 → clean) | Raise `wannier.nbnd` by a few bands, rerun from nscf with `--reset-epw` |
| α²F λ ≫ mode-summed λ at small phonon smearing | Fine q-grid too coarse; α²F is spiky | scflow uses the narrowest smearing on the λ plateau; densify `fine_qgrid` |
| λ=50 or other unphysical value | nbndsub too small; distorted Fermi-surface WFs | Expand nbndsub to span full Fermi manifold (see NbN) |
| 50% NaN pattern in .epmatwp | Degenerate d-bands → IEEE_INVALID in Bloch-to-Wannier | Use ph.x `electron_phonon='simple'` (see MgCNi3) |
| EOF reading .epmatwp | Wannier windows changed between pass1 and pass2 | `--reset-epw --start-from epw1` |
| `egnv` = 0 bytes | dis_froz_max doesn't span fine Ef | Set `use_a2f: true` for Allen-Dynes from α²F |
| Fine-grid Ef shift > 0.5 eV | Poor Wannier interpolation | Change projections (e.g. sp3 → explicit s+p) |
| Eliashberg Broyden factorisation crash | All temperatures above Allen-Dynes Tc | Lower `temperatures_K` below Tc |
| ph.x SIGSEGV before irreps table | Hexagonal space group on Apple Silicon | Run on Linux; avoid hexagonal on Mac |
| ph.x slow convergence (100+ iters/rep) | alpha_mix too aggressive for d-metals | Set `alpha_mix(1)=0.1`, `nmix_ph=8`, `tr2_ph=1.0d-12` |
| Imaginary phonon modes | FIRST suspect the structure: wrong Wyckoff site, wrong lattice/pressure | `scflow validate`: check `structure.distance.*` and `scf.pressure`. Only if both pass: genuine anharmonicity (PdH) → SSCHA |
| Hydride Tc << literature | Wrong structure, unconverged coarse grids, or Allen-Dynes at λ>1.5 | Check pressure first; compare the gate's Eliashberg Tc with Allen-Dynes |
| Tc only matches experiment at μ*≈0.18+ | λ overestimated (unconverged) | Do not tune μ*; converge λ (q-grid, fine grid, disentanglement) |
