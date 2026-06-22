"""Compare per-event C++ vs Python feature dumps for HiggsBDT.

Inputs (default: in cwd):
    cpp_feats_dump.csv  (produced by HiggsBDT.cc temporary dump block)
    py_feats_dump.csv   (produced by scripts/dump_py_feats.py)

Matching: events are paired by rounding (m4l, mZ1, mZ2) to 3 d.p. — the same
key is used on both sides. The C++ and Python skims do NOT preserve event
order, so we cannot rely on position.

Output: PASS/FAIL table per feature with max |Δ| and Pearson r.

Critical checks called out separately:
  * Helicity slots (cos_theta_star..Phi1) — sign-flip detector.
  * Phi, Phi1 — also mirror check (correlation ~ -1 means a plane-normal
    cross-product was reversed somewhere).

The 23 mechanical features (slots 0..17) should match to ~1e-3 once the
mechanical reconstruction is right; helicity angles 18..22 are the actual
test of the C++ port. Threshold for PASS is configurable.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FEATURES = [
    "m4l", "mZ1", "mZ2",
    "pt4l", "eta4l",
    "pt_Z1", "pt_Z2", "dR_Z1Z2",
    "pt_mu1", "pt_mu2", "pt_mu3", "pt_mu4",
    "eta_mu1", "eta_mu2", "eta_mu3", "eta_mu4",
    "pt_Z1_over_m4l", "pt_Z2_over_m4l",
    "cos_theta_star", "cos_theta1", "cos_theta2", "Phi", "Phi1",
]
HELI_SLOTS = list(range(18, 23))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpp", type=Path, default=Path("cpp_feats_dump.csv"))
    p.add_argument("--py",  type=Path, default=Path("py_feats_dump.csv"))
    p.add_argument("--tol-mech", type=float, default=1e-2,
                   help="max |Δ| tolerance for mechanical (mass/pT/eta) features")
    p.add_argument("--tol-heli", type=float, default=1e-2,
                   help="max |Δ| tolerance for helicity features")
    p.add_argument("--key-decimals", type=int, default=3)
    return p.parse_args()


def load_with_key(path: Path, key_decimals: int) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"missing file: {path}")
    df = pd.read_csv(path)
    for col in ("m4l", "mZ1", "mZ2"):
        if col not in df.columns:
            sys.exit(f"{path} is missing column {col!r}")
    df["__key"] = (
        df["m4l"].round(key_decimals).astype(str) + "|" +
        df["mZ1"].round(key_decimals).astype(str) + "|" +
        df["mZ2"].round(key_decimals).astype(str)
    )
    df = df.drop_duplicates("__key", keep="first").reset_index(drop=True)
    return df


def main() -> None:
    args = parse_args()
    cpp = load_with_key(args.cpp, args.key_decimals)
    py  = load_with_key(args.py,  args.key_decimals)

    merged = cpp.merge(py, on="__key", suffixes=("_cpp", "_py"))
    n_cpp, n_py, n_match = len(cpp), len(py), len(merged)
    print(f"events: cpp={n_cpp}  py={n_py}  matched={n_match}")
    if n_match == 0:
        sys.exit("no overlap — re-dump from the same MC sample")
    if n_match < 20:
        print(f"WARN: only {n_match} matched events; statistics will be noisy")

    rows = []
    overall_pass = True
    for i, feat in enumerate(FEATURES):
        c_col = f"{feat}_cpp"
        p_col = f"{feat}_py"
        if c_col not in merged.columns or p_col not in merged.columns:
            print(f"  skipping {feat}: not present in both dumps")
            continue
        x_cpp = merged[c_col].to_numpy(dtype=np.float64)
        x_py  = merged[p_col].to_numpy(dtype=np.float64)
        diff = x_cpp - x_py
        max_abs = float(np.max(np.abs(diff)))
        # Pearson correlation; if either side is constant, fall back to nan.
        if np.std(x_cpp) > 0 and np.std(x_py) > 0:
            r = float(np.corrcoef(x_cpp, x_py)[0, 1])
        else:
            r = float("nan")
        tol = args.tol_heli if i in HELI_SLOTS else args.tol_mech
        ok = max_abs < tol
        flag = "PASS" if ok else "FAIL"
        # Sign-flip diagnostic: correlation near -1 with similar magnitude.
        signflip = (r < -0.9)
        mirror_note = " (sign-flip?)" if signflip else ""
        rows.append((feat, max_abs, r, tol, flag + mirror_note))
        if not ok:
            overall_pass = False

    # Pretty-print.
    print(f"\n{'feature':18s} {'max|Δ|':>12s} {'corr':>8s} {'tol':>10s}  status")
    print("-" * 64)
    for feat, max_abs, r, tol, status in rows:
        print(f"{feat:18s} {max_abs:12.4e} {r:8.4f} {tol:10.0e}  {status}")
    print("-" * 64)
    print("HELICITY block (slots 18-22) — actual port test:")
    for feat, max_abs, r, tol, status in rows[18:23]:
        print(f"  {feat:18s} {max_abs:12.4e} corr={r:+.4f}  {status}")

    if overall_pass:
        print("\nOVERALL: PASS — C++ feature pipeline matches Python.")
        sys.exit(0)
    else:
        print("\nOVERALL: FAIL — fix the failing features in HiggsBDT.cc and re-dump.")
        sys.exit(1)


if __name__ == "__main__":
    main()
