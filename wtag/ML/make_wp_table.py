#!/usr/bin/env python3
"""Working-point table for the v2 W taggers (DeepAK8/ParticleNet convention:
WPs at fixed inclusive QCD mistag rate on the test split).

For each tagger and eps_b in {5%, 2.5%, 1%, 0.5%}: score threshold, inclusive
signal efficiency, and per-pT-bin signal efficiency / actual mistag.

Outputs results/wp_table.md, results/wp_table.tex, results/wp_table.json.
"""
import json
import os

import numpy as np
import pandas as pd

import evaluate as ev

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, "results")
WPS = [(0.05, "5\\%"), (0.025, "2.5\\%"), (0.01, "1\\%"), (0.005, "0.5\\%")]
PT_BINS = ev.PT_BINS
TAGGERS = ["part_nom", "bdt_nom", "part_md", "bdt_md", "pnet_wvsqcd"]
LABELS = {"part_nom": "ParT-lite", "bdt_nom": "BDT",
          "part_md": "ParT-lite MD", "bdt_md": "BDT MD",
          "pnet_wvsqcd": "ParticleNet WvsQCD"}


def main():
    meta = json.load(open(os.path.join(BASE, "data", "meta.json")))
    df = pd.read_parquet(os.path.join(BASE, "data", "dataset.parquet"))
    df = df[(df.split == "test") & (df.y >= 0)].copy()
    w = df.weight.abs().values
    y = df.y.values

    scores = {"bdt_nom": ev.score_bdt(df, "nom", meta),
              "bdt_md": ev.score_bdt(df, "md", meta),
              "part_nom": ev.score_part(df, "nom"),
              "part_md": ev.score_part(df, "md"),
              "pnet_wvsqcd": df.pnet_wvsqcd.values}

    out = {}
    md = ["| tagger | WP (mistag) | threshold | eps_S incl | "
          + " | ".join(f"eps_S {lo}-{hi}" for lo, hi in PT_BINS) + " |",
          "|---" * (4 + len(PT_BINS)) + "|"]
    tex = [
        "\\begin{tabular}{llcc" + "c" * len(PT_BINS) + "}",
        "\\hline",
        "tagger & WP ($\\varepsilon_B$) & threshold & $\\varepsilon_S$ (incl.) & "
        + " & ".join(f"$\\varepsilon_S^{{{lo}\\text{{--}}{hi}}}$" for lo, hi in PT_BINS)
        + " \\\\", "\\hline",
    ]
    for name in TAGGERS:
        s = scores[name]
        sb = s[y == 0]
        wb = w[y == 0]
        order = np.argsort(-sb)
        cw = np.cumsum(wb[order]) / wb.sum()
        out[name] = {}
        for eb, eb_lab in WPS:
            thr = float(sb[order][np.searchsorted(cw, eb)])
            row = {"threshold": thr,
                   "eps_s_incl": float(np.average(s[y == 1] > thr, weights=w[y == 1]))}
            for lo, hi in PT_BINS:
                m = (df.pt > lo) & (df.pt < hi)
                ms, mb = m.values & (y == 1), m.values & (y == 0)
                row[f"eps_s_{lo}_{hi}"] = float(np.average(s[ms] > thr, weights=w[ms]))
                row[f"eps_b_{lo}_{hi}"] = float(np.average(s[mb] > thr, weights=w[mb]))
            out[name][str(eb)] = row
            md.append(f"| {LABELS[name]} | {eb:.1%} | {thr:.4f} | "
                      f"{row['eps_s_incl']:.3f} | "
                      + " | ".join(f"{row[f'eps_s_{lo}_{hi}']:.3f}" for lo, hi in PT_BINS)
                      + " |")
            tex.append(f"{LABELS[name]} & {eb_lab} & {thr:.3f} & "
                       f"{row['eps_s_incl']:.3f} & "
                       + " & ".join(f"{row[f'eps_s_{lo}_{hi}']:.3f}" for lo, hi in PT_BINS)
                       + " \\\\")
        md.append("| | | | | " + " | ".join("" for _ in PT_BINS) + " |")
        tex.append("\\hline")
    tex.append("\\end{tabular}")

    with open(os.path.join(RES, "wp_table.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(RES, "wp_table.tex"), "w") as f:
        f.write("\n".join(tex) + "\n")
    with open(os.path.join(RES, "wp_table.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\n".join(md))


if __name__ == "__main__":
    main()
