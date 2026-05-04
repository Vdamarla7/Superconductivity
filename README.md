# Hydride Superconductor Pipeline

A Python pipeline for computing superconducting critical temperatures (Tc) of hydride compounds from first principles, using Quantum ESPRESSO + EPW + Wannier90.

**Current status**: 6 compounds validated across 3 compound classes (hydrides, d-metal nitrides/elemental, A15 ternaries). Nb₃Sn in progress.

---

## What this does

Given a crystal structure and pressure, the pipeline runs:

```
pw.x (vc-relax) → pw.x (SCF) → ph.x (phonons) → q2r + matdyn
    → pw.x (NSCF) → epw.x (Wannierize) → epw.x (Eliashberg Tc)
```

Each step is checkpointed — if you restart, completed steps are skipped automatically.

---

## Validated results

| Compound | Pressure | Our Tc | Literature | Notes |
|---|---|---|---|---|
| H₃S | 200 GPa | **216 K** | ~200–220 K (exp: 203 K) | ✓ Full EPW pipeline |
| H₃Se | 150 GPa | **97 K** | ~100–130 K | ✓ Full EPW pipeline |
| CaH₆ | 200 GPa | 22.6 K | ~220–235 K | ⚠ Harmonic only — SSCHA required |
| Nb | 0 GPa | **11.4 K** | 9.3 K (exp) | ✓ d-band, a2f path |
| NbN | 0 GPa | **16.6 K** (μ\*=0.18) | 12–16 K | ✓ d-band, nbndsub=13 required |
| MgCNi3 | 0 GPa | **7.9 K** | 7–8 K (exp) | ✓ EPW fails (d-band NaN); ph.x workaround |
| Nb₃Sn | 0 GPa | — | 18.3 K (exp) | 🔄 In progress (phonons running) |
| PdH | 0 GPa | — | ~9 K (exp) | ⚠ Imaginary phonons at X-point (harmonic DFT) |

See [results/](results/) for per-compound summaries and methodology notes.

---

## Repository layout

```
pipeline/
  pipeline.py       # Orchestrator — runs all 6 steps, checkpointing, arg parsing
  generators.py     # QE/EPW input file generators (one function per step)
  postprocess.py    # Output parsers: Fermi energy, Tc, sym ops, imaginary modes
  kpoints.py        # K-point grid utilities (explicit crystal list for NSCF/EPW)
  sscha_relax.py    # Anharmonic phonon renormalization via SSCHA (for CaH6, LaH10)

compounds/
  template.yaml       # Annotated template — start here for a new compound
  H3S_200GPa.yaml     # H₃X family (sp physics)
  H3Se_150GPa.yaml
  CaH6_200GPa.yaml    # Clathrate hydride (SSCHA required)
  Nb_0GPa.yaml        # Elemental d-metal (minimal Wannier basis)
  NbN_0GPa.yaml       # Rock-salt d-metal nitride (nbndsub=13)
  MgCNi3_0GPa.yaml    # Ternary anti-perovskite (EPW fails; ph.x workaround)
  PdH_0GPa.yaml       # Interstitial hydride (imaginary phonons in harmonic DFT)
  Nb3Sn_0GPa.yaml     # A15 ternary (44 WFs; in progress)

results/
  H3S_200GPa/summary.md
  H3Se_150GPa/summary.md
  CaH6_200GPa/summary.md   # includes SSCHA explanation

.claude/skills/simulate-superconductor/SKILL.md
  # Claude skill with full pipeline know-how — load with /simulate-superconductor
```

---

## Dependencies

**External (not included):**
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) 7.x compiled with EPW and Wannier90
- [pslibrary](https://github.com/dalcorso/pslibrary) pseudopotentials (PBE USPP)

**Python:**
```bash
pip install pyyaml
# Optional for SSCHA: pip install python-sscha cellconstructor
```

Set paths in your compound YAML:
```yaml
paths:
  qe_bin:     "/path/to/qe-7.x/bin"
  pseudo_dir: "/path/to/pslibrary/pbe/PSEUDOPOTENTIALS/"
```

---

## Running

```bash
# Run full pipeline for H3S at 200 GPa
python pipeline/pipeline.py compounds/H3S_200GPa.yaml

# Resume from a specific step
python pipeline/pipeline.py compounds/H3S_200GPa.yaml --start-from epw1

# Check what inputs would be generated without running anything
python pipeline/pipeline.py compounds/H3S_200GPa.yaml --dry-run
```

Output goes to `runs/<prefix>/` with subdirectories `1_vcrelax/`, `2_phonons/`, `3_epw/`.

---

## Adding a new compound

1. Copy `compounds/template.yaml` → `compounds/MyCompound.yaml`
2. Fill in structure, pseudopotentials, Wannier parameters, and paths
3. Run `--dry-run` to check generated inputs
4. Verify sym ops after SCF: `grep "Sym. Ops." runs/<prefix>/2_phonons/scf.out` (expect 48 for Im-3m)
5. Compute correct `dis_froz_max` from NSCF eigenvalues before committing to EPW

The template YAML has inline comments for every parameter. See also the [Claude skill](.claude/skills/simulate-superconductor/SKILL.md) for detailed Wannier window rules, common error fixes, and compound-specific guidance.

---

## Key design decisions

**ibrav=0 for everything except simple H₃X**  
`ibrav=3` only works when all atomic positions use 0 and 0.5 fractional coordinates. CaH₆ has H at ±0.25 — using ibrav=3 gives 4 sym ops instead of 48, producing wrong phonons and garbage Tc. Use ibrav=0 + explicit CELL_PARAMETERS for all other structures; QE does full symmetry detection from coordinates.

**Explicit K_POINTS crystal for NSCF**  
EPW requires the full N³ k-grid, not the symmetry-reduced set. `K_POINTS {automatic}` causes "inconsistent k-grids" in Wannierization. `kpoints.py` generates the full explicit list.

**Two-pass EPW**  
Pass 1 (Wannierize) writes `.epmatwp` — the e-ph matrix in the Wannier basis (~64–500 MB). Pass 2 reads it and does fine-grid interpolation + Eliashberg. If you change Wannier windows between passes, delete `.epmatwp` and rerun pass 1.

**Pre-create `tmp/{prefix}.ephmat/`**  
EPW's `write_freq()` opens a file in this directory before creating it. The pipeline pre-creates the directory; for manual runs do `mkdir -p runs/<prefix>/3_epw/tmp/<prefix>.ephmat/`.

---

## CaH₆ and SSCHA

CaH₆ (and LaH₁₀, YH₆) are quantum anharmonic crystals: their H-cage phonons have imaginary harmonic frequencies at all pressures in PBE. This is real physics — nuclear quantum effects stabilize them. The harmonic EPW result (Tc≈23 K) is unphysical. SSCHA renormalizes the phonon spectrum and recovers Tc≈220 K.

```bash
# Run SSCHA on CaH6 (requires python-sscha + cellconstructor)
python pipeline/sscha_relax.py compounds/CaH6_200GPa.yaml \
    --nconfigs 50 --niter 20 --temp 300

# Then copy renormalized dyn files and rerun from phonons step
python pipeline/pipeline.py compounds/CaH6_200GPa.yaml --start-from nscf
```

Runtime estimate on Mac: ~17 hours (50 configs × 10 iters, 448-atom supercell).  
Use AWS for production runs.

---

## Platform notes

**Apple Silicon (M-series Mac):**  
QE 7.4.1 ph.x segfaults for hexagonal space groups (P6/mmm, P6₃/mmc) on ARM64. All cubic structures (Im-3m, Fm-3m, F-43m) work fine. Affected: MgB₂ and all hexagonal compounds. Use Linux/AWS for hexagonal.

---

## Phase 2: novel compound search (planned, AWS)

The goal is to find new ternary hydrides that are:
- Superconducting with Tc > 20 K
- Dynamically stable at ambient pressure (no imaginary phonons at 0 GPa)
- Synthesizable: metastable after decompression from 20–50 GPa

**Target chemical space:** A₂MH₆ / A₃MH₆ families where M ∈ {Ir, Rh, Os, Ru, Pt, Pd} and A ∈ {Mg, Ca, Al, Li, Na}. These form stable octahedral [MH₆]ⁿ⁻ complexes — the key structural motif for ambient-pressure metastable superconductors (prototype: Mg₄Pt₃H₆, Tc=2.9 K).

**Screening pipeline (per candidate):**
```
AIRSS random structure search at 0 + 20–50 GPa
  → pw.x fast relaxation
  → ph.x on 2×2×2 (dynamic stability filter — discard if imaginary at 0 GPa)
  → pymatgen + MP API (convex hull distance — flag if > 150 meV/atom)
  → NEB decomposition barrier (kinetic stability — keep if > 0.3 eV/atom)
  → Full EPW Tc (finalists only)
```

**AWS budget estimate:** 100 candidates screened in ~1 week wall time (fully parallel), full EPW for top 5–10 candidates in ~2–5 days each.

---

## References

- Drozdov et al., Nature **525**, 73 (2015) — H₃S experimental Tc=203 K
- Wang et al., Phys. Rev. B **92**, 224507 (2015) — CaH₆ theoretical prediction
- Poncé et al., Comput. Phys. Commun. **209**, 116 (2016) — EPW
- Giustino et al., Phys. Rev. B **76**, 165108 (2007) — EPW method
- Monacelli et al., J. Phys.: Condens. Matter **33**, 363001 (2021) — SSCHA review
