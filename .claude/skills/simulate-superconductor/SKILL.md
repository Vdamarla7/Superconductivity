---
name: simulate-superconductor
description: Run the full EPW superconductivity pipeline for a hydride compound — from crystal structure to Eliashberg Tc. Handles structure setup, pseudopotential selection, Wannier window tuning, and pipeline execution.
argument-hint: "[chemical-formula] [pressure_GPa] OR [path/to/compound.yaml]"
---

You are helping the user simulate a new hydride superconductor using the QE + EPW pipeline at `/Users/vdark/work/quantum/Simulations/pipeline/`.

The argument is either:
- A compound formula + pressure: e.g. `H3Se 150` or `CaH6 150`
- A path to an existing YAML config: e.g. `compounds/CaH6_150GPa.yaml`

---

## OVERVIEW

The pipeline lives at `/Users/vdark/work/quantum/Simulations/pipeline/` and consists of:
- `pipeline.py` — orchestrator (6 steps: vcrelax → scf → phonons → nscf → epw1 → epw2)
- `generators.py` — QE/EPW input file generators
- `postprocess.py` — output parsers
- `compounds/` — YAML configs (one per compound)
- `runs/<prefix>/` — output for each run

Run with: `python pipeline.py compounds/MyCompound.yaml [--start-from STEP] [--dry-run]`

**IMPORTANT**: Each compound needs a unique `prefix` — two compounds with the same prefix share `runs/<prefix>/` and will use each other's cached outputs.

---

## STEP 1 — DETERMINE STRUCTURE AND ibrav

**The single most important decision: ibrav=3 vs ibrav=0**

### Use ibrav=3 ONLY when:
- Space group is Im-3m (BCC)
- AND all atomic positions use only 0 and 0.5 fractional coordinates

This covers: H3S, H3Se, H3Te — the simple H₃X family.

### Use ibrav=0 for EVERYTHING ELSE including:
- CaH6 (H at 12e Wyckoff: ±1/4,0,0 → fractional coords 0.25/0.75)
- LaH10, YH6, clathrate hydrides, any non-trivial Wyckoff site
- Any compound where you're not 100% sure about ibrav=3 compatibility

With ibrav=0, QE does full symmetry detection from coordinates → always finds correct sym ops.

### How to set up ibrav=0:
1. Get BCC primitive vectors from conventional cube edge A:
   ```
   a1 = A/2 * (-1, 1, 1)
   a2 = A/2 * (1, -1, 1)
   a3 = A/2 * (1, 1, -1)
   ```
2. Verify with spglib before running:
   ```python
   import spglib
   ds = spglib.get_symmetry_dataset((cell_matrix, positions, numbers), symprec=1e-3)
   print(ds['international'], len(ds['rotations']))  # must show Im-3m, 48
   ```
3. Add `skip_vcrelax: true` if using literature cell parameters

### Sym ops check — always verify after SCF:
```bash
grep "Sym\. Ops\." scf.out
```
Must show expected number (48 for Im-3m). Wrong sym ops = wrong phonons = wrong Tc.

---

## STEP 2 — CHOOSE PSEUDOPOTENTIALS

Always USPP (rrkjus) for production:
```bash
ls /Users/vdark/work/quantum/pslibrary/pbe/PSEUDOPOTENTIALS/ | grep -i "^ELEMENT\."
```
- Standard elements: `pbe/PSEUDOPOTENTIALS/`
- Heavy elements (Ir, Pt, Os, Rh, Ru, Re): `rel-pbe/PSEUDOPOTENTIALS/` (relativistic)
- Pattern: `{Element}.pbe-{n|dn|spn}-rrkjus_psl.1.0.0.UPF`

---

## STEP 3 — SET WANNIER PARAMETERS

This is the most compound-specific and failure-prone step.

### nbndsub and projections:
| Compound type | Formula | nbndsub |
|---|---|---|
| H₃X family | S/Se/Te:sp3 + nH:s | 4 + nH |
| Ca/Mg hydrides | M:d + nH:s | 5 + nH |
| La/Y clathrates | M:d + nH:s | 5 + nH |

Use explicit orbitals (`"Se:s"` + `"Se:p"`) instead of `"Se:sp3"` if sp3 gives NaN coupling.

### Deep core bands in d-metal nitrides (NbN, VN, TiN):
N-2s (~-35 eV) and N-2p (~-11 eV) appear as low-index bands below the valence manifold. Use `bands_skipped: "exclude_bands = 1-4"` and set `dis_win_min` just below band 5. Projections Nb:s+Nb:d+N:p = 9 WFs (note Nb:s required to span the full band manifold from 3 eV to Ef).

### Windows — CRITICAL RULES:

**Rule 1: Run SCF first, get Ef from scf.out.**
**Rule 2: Run NSCF (or let pipeline run it), get fine-grid Ef from nscf.out.**
**Rule 3: Fine-grid Ef often shifts 0.3–0.8 eV above coarse Ef due to Wannier interpolation.**

```
dis_win_min:  Ef - 15 eV   (well below all relevant bands)
dis_win_max:  Ef + 20 eV   (must contain ≥ nbndsub bands at EVERY k-point)
dis_froz_min: same as win_min
dis_froz_max: Ef + 3 eV    (based on SCF Ef; must NOT contain > nbndsub bands at any k)
```

**To find correct dis_froz_max — run this after NSCF:**
```python
import re
content = open("3_epw/nscf.out").read()
Ef = float(re.search(r'the Fermi energy is\s+([\d.+-]+)', content).group(1))
blocks = re.findall(r'bands \(ev\):\s*([\d\s.+-]+?)(?=\n\s*\n|\n\s*k =)', content, re.DOTALL)
for e_test in [Ef+2, Ef+3, Ef+4, Ef+5]:
    max_n = max(sum(1 for e in map(float,b.split()) if e < e_test) for b in blocks)
    print(f"dis_froz_max={e_test:.1f}: max {max_n} bands (need ≤ {nbndsub})")
```

**Pass1 and pass2 MUST use identical Wannier windows.** If you change windows, delete `.epmatwp` and rerun pass1.

### No-frozen-window (pure disentanglement) for entangled d-band systems:

**Problem**: For systems like NbN where the Fermi level sits inside a dense d-band complex, setting `dis_froz_max` above Ef causes too many bands (> nbndsub) to enter the frozen window at some k-points. Setting it below Ef means the Fermi-crossing bands are free to be distorted by Wannier90.

**Key insight about egnv**: egnv is ONLY written when the frozen window spans the Fermi level — i.e., `dis_froz_max > Ef`. Without the frozen window spanning Ef, the Wannier-interpolated bands on the fine grid won't have the right band structure for the egnv Eliashberg approach. This affects any system where the frozen window can't span Ef.

**Solution for d-band systems**: Use `a2f=.true.` + `eliashberg=.false.` to get Allen-Dynes Tc from the spectral function directly. This bypasses the egnv/ephmat requirement entirely:

```yaml
dis_win_min:  3.0
dis_win_max:  30.0
dis_froz_min: 3.0
dis_froz_max: 3.0   # = dis_froz_min → no frozen window (pure disentanglement)

phonselfen: true
a2f: true           # ← compute α²F(ω) → λ → Allen-Dynes Tc (no ephmat/egnv needed)
eliashberg: false   # ← don't solve full Eliashberg (avoids egnv EOF crash)
```

EPW will print "Electron-phonon coupling strength = X.XX" and "Estimated Allen-Dynes Tc = YY K".

This works when: (a) the system has strongly entangled bands crossing Ef where frozen window can't span Ef, and (b) Allen-Dynes Tc is sufficient (not full Eliashberg). For full Eliashberg in d-band systems, the frozen window MUST span Ef (dis_froz_max > Ef_fine), which requires finding a narrow enough band manifold at Ef.

### Validation in epw_pass2.out:
- `Electron-phonon coupling strength = X.XX` (real number, not NaN)
- `N bands within the Fermi window` > 0
- Fine-grid Ef shift < 0.5 eV vs coarse (larger shift → Wannier interpolation poor)
- `egnv` file in `tmp/{prefix}.ephmat/` must be non-zero bytes after pass2

---

## STEP 4 — EPW PARAMETERS

```yaml
epw:
  fine_kgrid:  [10, 10, 10]   # PoC; [20,20,20] for production
  fsthick:     1.0             # eV — start wide (1.0-1.5); narrow only if too slow
  wscut:       X.X             # 10× max phonon frequency in eV:
                               # H₃X: 2.0–3.0 eV; clathrates (CaH6, LaH10): 3.0–4.0 eV
```

Temperature scan: bracket the expected Tc. If unknown, scan 50–300 K.

---

## STEP 5 — CREATE YAML AND RUN

1. Copy `compounds/template.yaml` → `compounds/{Formula}_{P}GPa.yaml`
2. For ibrav=0: set `skip_vcrelax: true` and provide `cell_parameters_angstrom`
3. Dry-run first: `python pipeline.py compounds/X.yaml --dry-run`
4. Verify sym ops after SCF before committing to phonons (phonons take hours)
5. Check NSCF Ef and compute correct dis_froz_max before EPW
6. Run: `python pipeline.py compounds/X.yaml [--start-from STEP]`

---

## ANHARMONIC PHONONS (SSCHA)

For compounds with imaginary harmonic modes (CaH6, LaH10, YH6), use SSCHA to renormalize phonons before EPW.

**Script**: `Simulations/pipeline/sscha_relax.py`

```bash
python sscha_relax.py compounds/CaH6_200GPa.yaml --nconfigs 50 --niter 20 --temp 300
```

**How it works**:
1. Loads harmonic dyn files from `runs/<prefix>/2_phonons/`
2. Calls `ForcePositiveDefinite()` to replace imaginary modes with small positive values (starting point only)
3. Generates random displaced supercell configurations (thermal sampling)
4. Runs pw.x force calculations on each configuration
5. Updates trial force constants via SSCHA gradient minimization
6. Repeats until gradient < 1e-4 or `--niter` reached
7. Saves renormalized dyn files to `runs/<prefix>/sscha/`

**Key parameters**:
- `--nconfigs 10`: PoC (fast, ~20 min/iter on Mac); `50`: production quality
- `--temp`: temperature in K — run at several T to check anharmonic T-dependence
- Supercell is auto-detected from dyn0 q-grid (must match phonon q-grid)

**After SSCHA**: copy `sscha/<prefix>_sscha_irred.dyn*` over harmonic dyn files, rerun q2r → matdyn → EPW.

**Runtime estimate**: 50 configs × ~2 min/config × 10 iterations = ~17 hours on Mac for CaH6 (448-atom supercell). Use AWS for production.

---

## VALIDATED RESULTS

| Compound | Pressure | Our Tc | Literature |
|----------|----------|--------|------------|
| H3S      | 200 GPa  | 216 K  | ~200–220 K ✓ |
| H3Se     | 150 GPa  | 97 K   | ~100–130 K ✓ |
| CaH6     | 200 GPa  | 22.6 K (harmonic only — SSCHA required) | ~220–235 K |
| NbN      | 0 GPa    | FAILED — requires nbndsub=13 for correct Wannierization (see NbN notes below) | ~12–16 K |

---

## NbN WANNIER NOTES (d-band system pitfalls)

NbN (Fm-3m rock salt) is a prototype d-band superconductor. It is much harder to Wannierize correctly than H₃X compounds:

**Band structure after excluding bands 1-4 (N-2s at -35 eV, N-2p at -11 eV):**
- Band 5: Nb-5s/N-s at ~3 eV
- Bands 6-8: Nb-t2g at ~17 eV
- Bands 9-11: Nb-eg at ~18.75 eV
- Bands 12-13: crossing Ef at ~20.6 eV (coarse Ef=20.57, fine Ef=20.85 eV)
- Bands 14+: above 21 eV

**The nbndsub=9 problem**: To get egnv written (needed for Eliashberg), dis_froz_max must be > Ef_fine=20.85 eV. But at dis_froz_max≥19.74, max bands=11 at some k-points (>9). So nbndsub=9 cannot have a frozen window spanning Ef.

**With no frozen window (dis_froz_max=dis_froz_min=3.0)**: Wannier90 minimizes spread but distorts the Fermi-surface bands. a2F run gives λ=50 (completely unphysical; should be ~0.8-1.0). ω_log=1 meV (should be ~30-50 meV).

**Correct setup for NbN**: nbndsub=13 with projections Nb:s(1)+Nb:p(3)+Nb:d(5)+N:s(1)+N:p(3)=13. Then dis_froz_max can be set to ~20.74 eV (where max bands = exactly 13 at worst k-point). This fully spans the Fermi-crossing manifold.

**For pipeline PoC on Mac**: NbN is not recommended. Choose H₃X compounds (clean s+p physics) or wait until AWS for d-band systems requiring larger WF subspaces.

---

## PLATFORM LIMITATIONS

### Apple Silicon (ARM64) — QE 7.4.1 ph.x hexagonal crash
ph.x segfaults during irreducible representation computation for hexagonal space groups (P6/mmm, P6_3/mmc, etc.) on Apple Silicon. Crash occurs after printing the q-point list, before printing the irreps table. Signal: `SIGSEGV` with no useful backtrace.

**Affected compounds**: MgB2, all hexagonal structures.
**Not affected**: cubic (Fm-3m, Im-3m, F-43m), tetragonal — all compounds run so far work.
**Workaround**: avoid hexagonal compounds on Mac. Use NbN (Fm-3m) or other cubic structures instead. For hexagonal, run on Linux/AWS.

---

## ALL KNOWN BUGS (fixed in current code)

| Bug | Fix |
|-----|-----|
| QE resolves `outdir=./tmp/` from invocation dir | `cwd=str(dirs["step"])` in all `run_qe()` calls |
| NSCF "inconsistent k-grids" | `gen_nscf()` uses explicit `K_POINTS {crystal}` N³ list |
| pp.py "No dyn0 file" | `_run_pp_py()` builds `_ph0/` symlinks before calling pp.py |
| EPW pass1 falsely FAILED | `job_done()` accepts `"Total program execution"` (EPW) + `"JOB DONE"` (pw/ph) |
| `A: null` crash | Generator raises clear error — always provide numeric A |
| EPW `write_freq` error 108 | Pre-create `tmp/{prefix}.ephmat/` before pass2 in `run_epw2()` |
| EOF reading .epmatwp | Pass1 and pass2 must have identical Wannier windows; different windows → delete .epmatwp and rerun pass1 |
| Stale `.ephmat/egnv` = 0 bytes | Delete `tmp/{prefix}.ephmat/` and `tmp/{prefix}.epb*` when changing fsthick |
| ibrav=3 wrong sym ops for non-0/0.5 positions | Use ibrav=0 + explicit CELL_PARAMETERS |
| Two compounds sharing prefix | Each compound needs unique prefix in YAML |

---

## WHICH COMPOUNDS NEED SSCHA

CaH6, LaH10, and other clathrate hydrides are **quantum anharmonic crystals** — their H-cage phonons require nuclear quantum effects to stabilize. In harmonic PBE:
- CaH6@150GPa: ~100 imaginary modes → dynamically "unstable"
- CaH6@200GPa: 98 imaginary modes → still unstable (all pressures in harmonic PBE)
- EPW harmonic result: λ≈1.1, Allen-Dynes Tc≈22 K (real Tc ~220 K)

The imaginary modes are NOT a pipeline error. They are real physics. The harmonic result gives correct λ from real modes, but Tc is suppressed because most coupling modes are imaginary. SSCHA renormalizes the phonon spectrum and recovers the correct Tc.

**For harmonic-only compounds** (H3S, H3Se, H3Te): pipeline gives accurate Tc directly.
**For clathrate hydrides** (CaH6, LaH10, YH6): pipeline gives λ + Allen-Dynes estimate; SSCHA needed for Tc.

---

## QUICK DIAGNOSIS

| Symptom | Cause | Fix |
|---------|-------|-----|
| 4 sym ops instead of 48 | ibrav=3 with 0.25/0.75 coords | Switch to ibrav=0 |
| EPW NaN coupling | Frozen window wrong | Compute dis_froz_max from NSCF eigenvalues |
| `egnv` file = 0 bytes | No k-points in Fermi window (nkfs=0) | Widen dis_froz_max; check fine-grid Ef shift |
| EOF reading .epmatwp | Windows changed between pass1/pass2 | Delete .epmatwp, rerun pass1 with new windows |
| `write_freq` error 108 | `tmp/{prefix}.ephmat/` dir missing | Pre-create dir before running EPW (pipeline does this automatically now; for manual single-pass runs: `mkdir -p tmp/{prefix}.ephmat/`) |
| ph.x SIGSEGV before irreps table | Hexagonal space group on Apple Silicon (QE 7.4.1 ARM64 bug) | Switch to cubic compound or run on Linux |
| Imaginary phonon modes | Pressure too low or wrong structure | Try higher pressure; check positions |
| Fine-grid Ef shift > 0.5 eV | Poor Wannier interpolation | Change projections (e.g. sp3 → s+p explicit) |
| 0 bands in Fermi window | Over-compressed compound | Try lower pressure (e.g. 200→150 GPa) |
| Eliashberg Broyden `factorization` crash | All temps > Allen-Dynes Tc (Δ→0, matrix singular) | Lower `temperatures_K` to below Tc — check "Estimated Allen-Dynes Tc" printed before the crash |
| ephmat1 = 0 bytes when temps < Tc | EPW `liso=.true.` skips ephmat write at low T (code path bug) | Run pass2 with original high temps to get coupling numbers; accept Broyden crash as expected |
| CaH6/LaH10 harmonic Tc << literature | Quantum anharmonic crystal — imaginary modes dominate | Harmonic EPW gives λ but Tc is unphysical; SSCHA required for correct phonon spectrum |
| `mem_size_eliashberg` error: exceeds max_memlt | Default max_memlt=2.85 GB too small for large fine grids | Set `max_memlt = 16.0` in pass2 input (already added to generator) |
| egnv = 0 bytes in two-pass EPW | `epwread=.true.` unconditionally skips the q-loop that writes egnv — unfixable with flags (confirmed: even with pre-created ephmat dir, kmaps=.true., ephwrite=.true. — egnv is still 0 with epwread=.true.) | Use single-pass EPW (`epwwrite=.true.`, `epwread=.false.`, `eliashberg=.true.` all in one run). Two-pass only needed for >500 MB .epmatwp. |
| EPW requires epwwrite or epwread = .true. | EPW validation check | Always set `epwwrite=.true.` in pass1, `epwread=.true.` in pass2; for single-pass use `epwwrite=.true.` |
| `efermig` cannot bracket Ef | dis_froz_max below Ef — Wannier subspace doesn't span Fermi level | Run eigenvalue scan; ensure dis_froz_max > Ef at ALL k-points |
| too few projection functions | nbndsub > number of projection orbitals | Count projections: Nb:d=5, N:p=3 → 8 total; nbndsub must equal projection count |
| egnv=0 despite single-pass, froz below Ef | egnv is only written when dis_froz_max spans Ef; d-band systems where froz window can't reach Ef will always have egnv=0 | Use a2f=.true. + eliashberg=.false. to get Allen-Dynes Tc from spectral function (no egnv needed) |
| λ=50 or other unphysical value from a2f | No frozen window → Wannier90 distorts Fermi-surface bands → wrong e-ph coupling; common when dis_froz_max < Ef | Frozen window MUST span Ef for physical results; if nbndsub constraint prevents this, increase nbndsub to match the full band manifold at Ef |
