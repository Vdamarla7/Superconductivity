"""Screen a crystal structure: superconductivity likelihood, Tc estimate, thermodynamic stability.

  python ml/predict.py structure.cif [more.cif ...]
  python ml/predict.py compounds/CaH6_200GPa.yaml          # scflow configs work too
  python ml/predict.py --json POSCAR

These are fast statistical screens trained on ambient-pressure experimental (SuperCon via 3DSC)
and DFT (Materials Project) data. They prioritise candidates for scflow; they do not replace it.
"""

import argparse
import json
import sys
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from pymatgen.core import Composition, Structure

sys.path.insert(0, str(Path(__file__).parent))
from features import featurize, structure_from_scflow_yaml  # noqa: E402

MODELS = Path(__file__).parent / "models"
DOMAIN_KEYS = ["volume_per_atom", "density", "packing_fraction", "MagpieData mean AtomicWeight",
               "MagpieData mean Electronegativity", "frac_H", "nsites"]


class Model:
    def __init__(self, name):
        meta = json.loads((MODELS / f"{name}.json").read_text())
        self.booster = lgb.Booster(model_file=str(MODELS / f"{name}.txt"))
        self.features, self.domain, self.metrics = meta["features"], meta["domain"], meta["metrics"]

    def __call__(self, row):
        x = pd.DataFrame([row]).reindex(columns=self.features).astype(float)
        return float(self.booster.predict(x)[0])

    def out_of_domain(self, row):
        flags = []
        for k in DOMAIN_KEYS:
            lo, hi = self.domain.get(k, (None, None))
            v = row.get(k)
            if lo is not None and v is not None and not (lo <= v <= hi):
                flags.append(f"{k}={v:.3g} outside training range [{lo:.3g}, {hi:.3g}]")
        return flags


def load_structure(path):
    p = Path(path)
    if p.suffix in (".yaml", ".yml"):
        s, cfg = structure_from_scflow_yaml(p)
        return s, float(cfg["compound"].get("pressure_GPa", 0.0))
    return Structure.from_file(p), 0.0


def known_label(comp):
    known = MODELS / "sc_known.csv"
    if not known.exists():
        return None
    df = pd.read_csv(known)
    target = comp.reduced_composition
    for _, r in df.iterrows():
        try:
            if Composition(r["formula_sc"]).reduced_composition.almost_equals(target, rtol=0.01):
                return {"formula": r["formula_sc"], "tc_K": r["tc"], "class": r["sc_class"]}
        except Exception:
            continue
    return None


def predict(path, subset="conventional"):
    s, pressure = load_structure(path)
    row = featurize(s)
    clf, reg = Model(f"sc_classifier_{subset}"), Model(f"tc_regressor_{subset}")
    hull = Model("ehull_regressor") if (MODELS / "ehull_regressor.txt").exists() else None
    p_sc = clf(row)
    tc = float(np.expm1(reg(row)))
    ehull = max(0.0, hull(row)) if hull else None
    warnings = clf.out_of_domain(row)
    if hull:
        warnings += [w for w in hull.out_of_domain(row) if w not in warnings]
    if pressure > 5:
        warnings.insert(0, f"structure is at {pressure:g} GPa; all training data are at ambient pressure, "
                           f"so the stability estimate refers to 0 GPa chemistry and is not meaningful here")
    return {
        "input": str(path),
        "formula": s.composition.reduced_formula,
        "p_superconductor_given_tested": round(p_sc, 3),
        "tc_estimate_K": round(tc, 1),
        "e_above_hull_estimate_eV_per_atom": round(ehull, 3) if ehull is not None else None,
        "known_in_supercon": known_label(s.composition),
        "applicability_warnings": warnings,
        "model_cv": {
            "sc_roc_auc_grouped": round(clf.metrics["grouped"]["roc_auc"], 3),
            "tc_mae_K_grouped": round(reg.metrics["grouped"]["mae_K"], 1),
            "ehull_mae_eV_grouped": round(hull.metrics["grouped"]["mae_eV_per_atom"], 3) if hull else None,
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="+")
    ap.add_argument("--subset", choices=["conventional", "all"], default="conventional",
                    help="conventional: trained on non-cuprate/iron/heavy-fermion classes (what scflow can verify)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    results = [predict(p, args.subset) for p in args.structures]
    if args.json:
        print(json.dumps(results, indent=2))
        return
    for r in results:
        print(f"\n  {r['formula']}  ({r['input']})")
        print(f"    P(superconductor | tested)  {r['p_superconductor_given_tested']:.2f}"
              f"    (grouped-CV ROC-AUC {r['model_cv']['sc_roc_auc_grouped']})")
        print(f"    Tc estimate                 {r['tc_estimate_K']:.1f} K"
              f"    (grouped-CV MAE {r['model_cv']['tc_mae_K_grouped']} K)")
        if r["e_above_hull_estimate_eV_per_atom"] is not None:
            print(f"    E above hull estimate       {r['e_above_hull_estimate_eV_per_atom']:.3f} eV/atom"
                  f"    (grouped-CV MAE {r['model_cv']['ehull_mae_eV_grouped']} eV/atom)")
        if r["known_in_supercon"]:
            k = r["known_in_supercon"]
            print(f"    SuperCon entry              {k['formula']}: Tc = {k['tc_K']} K ({k['class']})")
        for w in r["applicability_warnings"]:
            print(f"    ⚠ {w}")


if __name__ == "__main__":
    main()
