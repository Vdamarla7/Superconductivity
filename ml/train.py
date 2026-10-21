"""Train and evaluate the screening models.

  python ml/train.py sc          # superconductivity: classifier + Tc regressor (3DSC_MP)
  python ml/train.py stability   # thermodynamic stability: E_above_hull regressor (Materials Project)

Evaluation uses 5-fold CV grouped by chemical system (of the undoped parent structure), so a
compound's doped variants and polymorphs never sit on both sides of a split. Random-split CV is
reported alongside it to show how much a naive split overstates performance.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from pymatgen.core import Composition, Structure
from sklearn.metrics import (average_precision_score, balanced_accuracy_score, mean_absolute_error, r2_score,
                             roc_auc_score)
from sklearn.model_selection import GroupKFold, KFold

sys.path.insert(0, str(Path(__file__).parent))
from features import chemsys, featurize_many  # noqa: E402

ML = Path(__file__).parent
DATA, MODELS = ML / "data", ML / "models"
CONVENTIONAL = {"Other", "Chevrel", "Carbon"}
# 3DSC builds doped entries as disordered supercells, and 81% of those are superconductors vs 52% of
# ordered entries. These columns encode how the dataset was made, not the material, so they are dropped.
DATASET_ARTIFACTS = ["is_ordered", "nsites"]
SEED = 0
JOBS = int(os.environ.get("ML_JOBS", "8"))

CLS_PARAMS = dict(objective="binary", n_estimators=600, learning_rate=0.03, num_leaves=31, min_child_samples=20,
                  subsample=0.8, subsample_freq=1, colsample_bytree=0.7, reg_lambda=1.0, verbose=-1, random_state=SEED, n_jobs=JOBS)
REG_PARAMS = dict(CLS_PARAMS, objective="regression_l1")


def _cached_features(name, structures):
    path = DATA / f"{name}_features.parquet"
    if path.exists():
        return pd.read_parquet(path)
    t = time.time()
    X = featurize_many(structures, n_jobs=JOBS)
    X.to_parquet(path)
    print(f"  featurized {len(X)} structures in {time.time() - t:.0f} s → {path}")
    return X


def _cv(model_fn, X, y, groups, kind):
    """Out-of-fold predictions for grouped and random 5-fold CV."""
    out = {}
    for label, splitter in (("grouped", GroupKFold(5)), ("random", KFold(5, shuffle=True, random_state=SEED))):
        pred = np.zeros(len(y))
        split = splitter.split(X, y, groups) if label == "grouped" else splitter.split(X, y)
        for tr, te in split:
            m = model_fn()
            m.fit(X.iloc[tr], y[tr])
            pred[te] = m.predict_proba(X.iloc[te])[:, 1] if kind == "cls" else m.predict(X.iloc[te])
        out[label] = pred
    return out


def _cls_metrics(y, p):
    return {"roc_auc": roc_auc_score(y, p), "pr_auc": average_precision_score(y, p),
            "balanced_acc@0.5": balanced_accuracy_score(y, p > 0.5), "positive_rate": float(y.mean())}


def _domain(X):
    return {c: [float(np.nanpercentile(X[c], 1)), float(np.nanpercentile(X[c], 99))] for c in X.columns}


def _save(name, model, features, domain, metrics):
    MODELS.mkdir(parents=True, exist_ok=True)
    model.booster_.save_model(MODELS / f"{name}.txt")
    (MODELS / f"{name}.json").write_text(json.dumps({"features": features, "domain": domain, "metrics": metrics},
                                                     indent=2, default=float))


def train_sc():
    df = pd.read_csv(DATA / "3dsc_mp.csv", comment="#")
    root = DATA / "3dsc_repo" / "superconductors_3D"
    structures = [Structure.from_file(root / c) for c in df["cif"]]
    X = _cached_features("3dsc", structures)
    ok = X.notna().any(axis=1).to_numpy()
    df, X = df[ok].reset_index(drop=True), X[ok].drop(columns=DATASET_ARTIFACTS).reset_index(drop=True)
    groups = df["pretty_formula_2"].map(lambda f: chemsys(Composition(f))).to_numpy()
    feats = list(X.columns)
    report = {"n": len(df), "n_groups": int(len(set(groups)))}

    for subset, mask in (("all", np.ones(len(df), bool)), ("conventional", df["sc_class"].isin(CONVENTIONAL).to_numpy())):
        Xs, ds, gs = X[mask].reset_index(drop=True), df[mask].reset_index(drop=True), groups[mask]
        y_cls = (ds["tc"] > 0).to_numpy().astype(int)
        oof = _cv(lambda: lgb.LGBMClassifier(**CLS_PARAMS), Xs, y_cls, gs, "cls")
        cls_m = {k: _cls_metrics(y_cls, v) for k, v in oof.items()}
        ordered = (~ds["synth_doped"]).to_numpy()
        cls_m["grouped_ordered_only"] = _cls_metrics(y_cls[ordered], oof["grouped"][ordered])

        pos = y_cls == 1
        y_reg = np.log1p(ds.loc[pos, "tc"].to_numpy())
        oof_r = _cv(lambda: lgb.LGBMRegressor(**REG_PARAMS), Xs[pos].reset_index(drop=True), y_reg, gs[pos], "reg")
        tc_true = np.expm1(y_reg)
        reg_m = {k: {"mae_K": mean_absolute_error(tc_true, np.expm1(v)), "r2_logTc": r2_score(y_reg, v),
                     "baseline_mae_K": mean_absolute_error(tc_true, np.full_like(tc_true, np.median(tc_true)))}
                 for k, v in oof_r.items()}
        o = ordered[pos]
        reg_m["grouped_ordered_only"] = {
            "mae_K": mean_absolute_error(tc_true[o], np.expm1(oof_r["grouped"][o])),
            "r2_logTc": r2_score(y_reg[o], oof_r["grouped"][o]),
            "baseline_mae_K": mean_absolute_error(tc_true[o], np.full(o.sum(), np.median(tc_true)))}
        report[subset] = {"n": int(mask.sum()), "classifier": cls_m, "tc_regressor": reg_m}

        clf = lgb.LGBMClassifier(**CLS_PARAMS).fit(Xs, y_cls)
        reg = lgb.LGBMRegressor(**REG_PARAMS).fit(Xs[pos], y_reg)
        _save(f"sc_classifier_{subset}", clf, feats, _domain(Xs), cls_m)
        _save(f"tc_regressor_{subset}", reg, feats, _domain(Xs), reg_m)

        imp = pd.Series(clf.booster_.feature_importance("gain"), index=feats).sort_values(ascending=False)
        report[subset]["top_features"] = list(imp.index[:10])

    known = df[["formula_sc", "tc", "sc_class", "material_id_2"]]
    known.to_csv(MODELS / "sc_known.csv", index=False)
    (MODELS / "sc_report.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(report, indent=2, default=lambda x: round(float(x), 3)))


def train_stability():
    df = pd.read_parquet(DATA / "mp.parquet")
    structures = [Structure.from_str(s, fmt="json") for s in df["structure"]]
    X = _cached_features("mp", structures)
    ok = X.notna().any(axis=1).to_numpy() & df["e_above_hull"].notna().to_numpy()
    df, X = df[ok].reset_index(drop=True), X[ok].reset_index(drop=True)
    groups = df["formula"].map(lambda f: chemsys(Composition(f))).to_numpy()
    y = np.clip(df["e_above_hull"].to_numpy(), 0, 1.0)
    oof = _cv(lambda: lgb.LGBMRegressor(**dict(REG_PARAMS, n_estimators=1500, num_leaves=127)), X, y, groups, "reg")
    stable = (df["e_above_hull"] <= 1e-6).to_numpy()
    near = (df["e_above_hull"] <= 0.1).to_numpy()
    metrics = {}
    for k, p in oof.items():
        metrics[k] = {"mae_eV_per_atom": mean_absolute_error(y, p),
                      "baseline_mae": mean_absolute_error(y, np.full_like(y, np.median(y))),
                      "stable_roc_auc": roc_auc_score(stable, -p),
                      "within_0.1eV_roc_auc": roc_auc_score(near, -p)}
    report = {"n": len(df), "n_groups": int(len(set(groups))), "stable_fraction": float(stable.mean()),
              "metrics": metrics}
    model = lgb.LGBMRegressor(**dict(REG_PARAMS, n_estimators=1500, num_leaves=127)).fit(X, y)
    _save("ehull_regressor", model, list(X.columns), _domain(X), metrics)
    (MODELS / "stability_report.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(report, indent=2, default=lambda x: round(float(x), 3)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("task", choices=["sc", "stability"])
    args = ap.parse_args()
    train_sc() if args.task == "sc" else train_stability()
