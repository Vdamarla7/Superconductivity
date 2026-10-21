"""Rank Materials Project compounds as superconductor candidates for scflow.

  python ml/screen.py                       # writes ml/models/candidates.csv
  python ml/screen.py --max-sites 8 --top 30

Filters (all from Materials Project DFT, not from a model): metallic, within --max-ehull of the
convex hull, at most --max-sites atoms per cell (scflow cost grows steeply with cell size), and
not already in SuperCon. The remainder is ranked by the conventional-superconductor classifier and
Tc regressor. The output is a shortlist for first-principles calculation, not a prediction of
superconductivity: the classifier's P(SC | tested) overstates the base rate for random compounds.
"""

import argparse
import json
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from pymatgen.core import Composition

ML = Path(__file__).parent
DATA, MODELS = ML / "data", ML / "models"
RADIOACTIVE = {"Tc", "Pm", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu"}


def reduced(formula):
    try:
        return Composition(formula).reduced_formula
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-ehull", type=float, default=0.025, help="eV/atom, from MP DFT")
    ap.add_argument("--max-sites", type=int, default=12)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--allow-radioactive", action="store_true")
    args = ap.parse_args()

    mp = pd.read_parquet(DATA / "mp.parquet", columns=["material_id", "formula", "e_above_hull", "is_metal",
                                                        "nsites", "spacegroup", "theoretical"])
    X = pd.read_parquet(DATA / "mp_features.parquet")
    dsc = pd.read_csv(DATA / "3dsc_mp.csv", comment="#", usecols=["formula_sc", "pretty_formula_2", "material_id_2"])
    # 3DSC often lists only doped variants, so exclude parent compositions and parent MP ids too.
    known = {reduced(f) for f in pd.concat([dsc["formula_sc"], dsc["pretty_formula_2"]])}
    known_ids = set(dsc["material_id_2"])

    keep = (mp["is_metal"].astype("boolean").fillna(False) & (mp["e_above_hull"] <= args.max_ehull)
            & (mp["nsites"] <= args.max_sites) & X.notna().any(axis=1)).to_numpy()
    mp, X = mp[keep].reset_index(drop=True), X[keep].reset_index(drop=True)
    mp["reduced"] = mp["formula"].map(reduced)
    novel = ~mp["reduced"].isin(known) & ~mp["material_id"].isin(known_ids)
    if not args.allow_radioactive:
        novel &= ~mp["formula"].map(lambda f: bool(RADIOACTIVE & {e.symbol for e in Composition(f).elements}))
    mp, X = mp[novel].reset_index(drop=True), X[novel].reset_index(drop=True)

    clf = lgb.Booster(model_file=str(MODELS / "sc_classifier_conventional.txt"))
    reg = lgb.Booster(model_file=str(MODELS / "tc_regressor_conventional.txt"))
    feats = json.loads((MODELS / "sc_classifier_conventional.json").read_text())["features"]
    mp["p_sc"] = clf.predict(X[feats].astype(float))
    mp["tc_est_K"] = np.expm1(reg.predict(X[feats].astype(float)))
    mp["score"] = mp["p_sc"] * mp["tc_est_K"]

    out = mp.sort_values("score", ascending=False).head(args.top)
    cols = ["material_id", "formula", "spacegroup", "nsites", "e_above_hull", "theoretical", "p_sc", "tc_est_K"]
    out[cols].to_csv(MODELS / "candidates.csv", index=False)
    print(f"{keep.sum()} metallic near-hull MP structures with ≤{args.max_sites} sites; "
          f"{len(mp)} not in SuperCon; top {args.top} → {MODELS / 'candidates.csv'}\n")
    print(out[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
