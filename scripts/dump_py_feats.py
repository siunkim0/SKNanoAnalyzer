"""Dump the BDT Python pipeline's 23 features + score for the C++ parity check.

Reads one of the BDT-project *_features.parquet files (Phase 2 output) and
runs the v5 xgboost model on it. Writes ``py_feats_dump.csv`` with the same
column order the C++ side uses (HiggsBDT.cc temporary dump block).

Run from the BDT project root, e.g.:

    cd /data6/Users/snuintern2/BDT
    python /data6/Users/snuintern2/folder/sknano/new/SKNanoAnalyzer/scripts/dump_py_feats.py \\
        --features data/ntuples/ggH_ZZ_4l_features.parquet \\
        --model    data/models/bdt_v5.json \\
        --out      /data6/Users/snuintern2/folder/sknano/new/SKNanoAnalyzer/py_feats_dump.csv \\
        --n        500

Match this run to the SAME MC sample you produce the C++ dump from. Path 1
training applied the SR mask m4l in [105, 140] at training time, so for an
apples-to-apples parity check with the C++ side (which still sees the full
training window pre-Task-2) we DO NOT pre-mask here — the C++ dump fires for
any event that reaches the ONNX call after Phase 5's m4l > m4l_min cut.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import pandas as pd
import xgboost as xgb

FEATURES = [
    "m4l", "mZ1", "mZ2",
    "pt4l", "eta4l",
    "pt_Z1", "pt_Z2", "dR_Z1Z2",
    "pt_mu1", "pt_mu2", "pt_mu3", "pt_mu4",
    "eta_mu1", "eta_mu2", "eta_mu3", "eta_mu4",
    "pt_Z1_over_m4l", "pt_Z2_over_m4l",
    "cos_theta_star", "cos_theta1", "cos_theta2", "Phi", "Phi1",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", required=True, type=Path,
                   help="Phase-2 parquet, e.g. data/ntuples/ggH_ZZ_4l_features.parquet")
    p.add_argument("--model", required=True, type=Path,
                   help="xgboost JSON model, e.g. data/models/bdt_v5.json")
    p.add_argument("--out", required=True, type=Path,
                   help="Output CSV (py_feats_dump.csv).")
    p.add_argument("--n", type=int, default=500,
                   help="Cap on number of events to dump (default 500).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.features.exists():
        sys.exit(f"feature parquet not found: {args.features}")
    if not args.model.exists():
        sys.exit(f"model not found: {args.model}")

    arr = ak.from_parquet(args.features)
    missing = [c for c in FEATURES if c not in arr.fields]
    if missing:
        sys.exit(f"feature parquet missing columns: {missing}")

    df = pd.DataFrame({c: ak.to_numpy(arr[c]) for c in FEATURES})
    df = df.head(args.n).reset_index(drop=True)

    X = df[FEATURES].to_numpy(dtype=np.float32)
    model = xgb.XGBClassifier()
    model.load_model(args.model)
    p_sig = model.predict_proba(X)[:, 1]

    df["py_score"] = p_sig
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} events, {len(df.columns)} cols)")


if __name__ == "__main__":
    main()
