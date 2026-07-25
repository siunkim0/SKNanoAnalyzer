#!/usr/bin/env python3
"""Train the XGBoost W-vs-QCD taggers (nominal + mass-decorrelated).

  python3 train_bdt.py            # both variants
  python3 train_bdt.py nom|md     # one variant

Outputs models/bdt_<variant>.json and models/bdt_<variant>_importance.json.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import xgboost as xgb

BASE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(BASE, "models"), exist_ok=True)

PARAMS = dict(
    objective="binary:logistic", eval_metric="auc", tree_method="hist",
    max_depth=6, eta=0.08, subsample=0.8, colsample_bytree=0.8,
    min_child_weight=50, device="cuda",
)
NROUNDS, EARLY = 2000, 50


def train(variant):
    meta = json.load(open(os.path.join(BASE, "data", "meta.json")))
    feats = meta["features"] if variant == "nom" else meta["features_md"]
    wcol = "w_flatpt" if variant == "nom" else "w_flatptm"

    df = pd.read_parquet(os.path.join(BASE, "data", "dataset.parquet"))
    df = df[df.y >= 0]          # drop t->bqq spectator rows
    tr, va = df[df.split == "train"], df[df.split == "val"]
    dtr = xgb.DMatrix(tr[feats], label=tr.y, weight=tr[wcol])
    dva = xgb.DMatrix(va[feats], label=va.y, weight=va[wcol])

    params = dict(PARAMS)
    try:
        bst = xgb.train(params, dtr, NROUNDS, evals=[(dva, "val")],
                        early_stopping_rounds=EARLY, verbose_eval=100)
    except xgb.core.XGBoostError:      # no GPU -> CPU fallback
        params["device"] = "cpu"
        bst = xgb.train(params, dtr, NROUNDS, evals=[(dva, "val")],
                        early_stopping_rounds=EARLY, verbose_eval=100)

    out = os.path.join(BASE, "models", f"bdt_{variant}.json")
    bst.save_model(out)
    imp = bst.get_score(importance_type="total_gain")
    with open(out.replace(".json", "_importance.json"), "w") as f:
        json.dump(dict(sorted(imp.items(), key=lambda kv: -kv[1])), f, indent=2)
    print(f"[{variant}] best_iter={bst.best_iteration} val_auc={bst.best_score:.5f} -> {out}")


if __name__ == "__main__":
    variants = sys.argv[1:] or ["nom", "md"]
    for v in variants:
        train(v)
