#!/usr/bin/env python3
"""Evaluate all taggers on the test split and produce the v2 deliverables.

Discriminants compared:
  bdt_nom / bdt_md          — XGBoost (models/)
  part_nom / part_md        — ParT-lite (models/)
  tau21 (inverted), tau21+mass-window, N2DDT (quantile-map scan)
  pnet_wvsqcd (nominal ceiling), pnet_xqqvsqcd (MD ceiling)

Outputs in results/:
  roc_inclusive.png, roc_ptbins.png     ROC curves (bkg rej vs sig eff)
  rejection_table.md                    bkg rejection @ eps_s = 0.3/0.5/0.7
  eff_vs_pt.png                         sig eff & mistag vs pT at fixed-mistag WPs
  sculpting.png, jsd_vs_effb.png        QCD m_SD shapes pass/fail + JSD scan
  summary.json                          all numbers
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

PT_BINS = [(200, 300), (300, 500), (500, 800), (800, 1200)]
MASS_WIN = (65.0, 105.0)
EPS_S_POINTS = (0.3, 0.5, 0.7)


# ---------------------------------------------------------------- scoring
def score_bdt(df, variant, meta):
    import xgboost as xgb
    feats = meta["features"] if variant == "nom" else meta["features_md"]
    bst = xgb.Booster()
    bst.load_model(os.path.join(BASE, "models", f"bdt_{variant}.json"))
    return bst.predict(xgb.DMatrix(df[feats]))


def score_part(df, variant):
    import torch
    from models import ParTLite
    from train_partlite import build_tokens
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    ck = torch.load(os.path.join(BASE, "models", f"part_{variant}.pt"),
                    map_location=dev, weights_only=False)
    model = ParTLite(ck["n_token_feat"], len(ck["gfeats"]))
    model.load_state_dict(ck["state_dict"])
    model.to(dev).eval()
    tok, mask = build_tokens(df, ck["md"])
    glob = ((df[ck["gfeats"]].values - ck["norm_mean"]) / ck["norm_std"]).astype(np.float32)
    glob = np.nan_to_num(glob, nan=0.0, posinf=0.0, neginf=0.0)  # mirror make_loader
    out = np.empty(len(df), dtype=np.float32)
    import torch as T
    with T.no_grad():
        for i in range(0, len(df), 65536):
            s = slice(i, i + 65536)
            out[s] = T.sigmoid(model(T.from_numpy(tok[s]).to(dev),
                                     T.from_numpy(mask[s]).to(dev),
                                     T.from_numpy(glob[s]).to(dev))).cpu().numpy()
    return out


def n2ddt_score(df, quantile_map):
    """N2 - map(rho, pt); tagged if < 0.  Returns -(n2 - map) so that
    'greater = more W-like' matches the other discriminants."""
    rho, pt = df["rho"].values, df["pt"].values
    ir = np.clip(np.digitize(rho, quantile_map["rho_edges"]) - 1, 0,
                 len(quantile_map["rho_edges"]) - 2)
    ip = np.clip(np.digitize(pt, quantile_map["pt_edges"]) - 1, 0,
                 len(quantile_map["pt_edges"]) - 2)
    return -(df["n2b1"].values - quantile_map["map"][ir, ip])


def build_n2_map(qcd, q):
    """q-quantile of N2 in (rho, pt) bins on the QCD training sample."""
    rho_edges = np.linspace(-6.0, -1.0, 15)
    pt_edges = np.array([200, 300, 400, 500, 650, 800, 1000, 1200], dtype=float)
    m = np.full((len(rho_edges) - 1, len(pt_edges) - 1), np.nan)
    ir = np.clip(np.digitize(qcd["rho"], rho_edges) - 1, 0, len(rho_edges) - 2)
    ip = np.clip(np.digitize(qcd["pt"], pt_edges) - 1, 0, len(pt_edges) - 2)
    for i in range(len(rho_edges) - 1):
        for j in range(len(pt_edges) - 1):
            v = qcd["n2b1"].values[(ir == i) & (ip == j)]
            if len(v) > 50:
                m[i, j] = np.quantile(v, q)
    # fill holes with global quantile
    m[np.isnan(m)] = np.quantile(qcd["n2b1"], q)
    return {"rho_edges": rho_edges, "pt_edges": pt_edges, "map": m}


# ---------------------------------------------------------------- metrics
def roc(y, score, w):
    order = np.argsort(-score)
    y, w = y[order], w[order]
    cs = np.cumsum(w * y) / np.sum(w * y)
    cb = np.cumsum(w * (1 - y)) / np.sum(w * (1 - y))
    return cs, cb  # eps_s, eps_b along threshold scan


def rej_at(eps_s_target, cs, cb):
    i = np.searchsorted(cs, eps_s_target)
    i = min(i, len(cb) - 1)
    return (1.0 / cb[i]) if cb[i] > 0 else np.inf


def jsd(p, q):
    p, q = p / p.sum(), q / q.sum()
    m = 0.5 * (p + q)
    def kl(a, b):
        mask = a > 0
        return np.sum(a[mask] * np.log2(a[mask] / b[mask]))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def main():
    meta = json.load(open(os.path.join(BASE, "data", "meta.json")))
    df = pd.read_parquet(os.path.join(BASE, "data", "dataset.parquet"))
    df = df[(df.split == "test") & (df.y >= 0)].copy()
    df["rho"] = np.log(np.maximum(df.sdmass, 1e-3) ** 2 / df.pt ** 2)
    df["tau21_inv"] = -df["tau21"]        # greater = more W-like
    mass_ok = (df.sdmass > MASS_WIN[0]) & (df.sdmass < MASS_WIN[1])
    df["tau21_masswin"] = np.where(mass_ok, -df["tau21"], -1e9)
    w = df.weight.abs().values            # physics (xsec) weights for evaluation
    y = df.y.values

    scores = {}
    for v in ("nom", "md"):
        p = os.path.join(BASE, "models", f"bdt_{v}.json")
        if os.path.exists(p):
            scores[f"bdt_{v}"] = score_bdt(df, v, meta)
        p = os.path.join(BASE, "models", f"part_{v}.pt")
        if os.path.exists(p):
            scores[f"part_{v}"] = score_part(df, v)
    scores["tau21"] = df["tau21_inv"].values
    scores["tau21+masswin"] = df["tau21_masswin"].values
    qcd_tr_mask = (df.y == 0)
    n2map = build_n2_map(df[qcd_tr_mask], 0.05)
    scores["n2ddt(5%)"] = n2ddt_score(df, n2map)
    scores["pnet_wvsqcd"] = df.pnet_wvsqcd.values
    scores["pnet_xqqvsqcd"] = df.pnet_xqqvsqcd.values

    summary = {"rejection": {}, "pt_bins": {}}

    # inclusive ROC + rejection table
    plt.figure(figsize=(7, 6))
    lines = ["| tagger | " + " | ".join(f"rej@eps_s={e}" for e in EPS_S_POINTS) + " |",
             "|---" * (1 + len(EPS_S_POINTS)) + "|"]
    for name, s in scores.items():
        cs, cb = roc(y, s, w)
        rejs = [rej_at(e, cs, cb) for e in EPS_S_POINTS]
        summary["rejection"][name] = dict(zip(map(str, EPS_S_POINTS), rejs))
        lines.append("| " + name + " | " + " | ".join(f"{r:,.0f}" for r in rejs) + " |")
        sel = cs > 0.05
        plt.plot(cs[sel], 1.0 / np.maximum(cb[sel], 1e-8), label=name)
    plt.yscale("log")
    plt.xlabel("signal efficiency")
    plt.ylabel("background rejection (1/eps_b)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.title(f"W vs QCD, AK8 pT 200-1200, {'+'.join(meta['eras'])}")
    plt.savefig(os.path.join(RES, "roc_inclusive.png"), dpi=140, bbox_inches="tight")
    plt.close()
    with open(os.path.join(RES, "rejection_table.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # per-pT-bin rejection at eps_s = 0.5
    for lo, hi in PT_BINS:
        m = (df.pt > lo) & (df.pt < hi)
        summary["pt_bins"][f"{lo}-{hi}"] = {}
        for name, s in scores.items():
            cs, cb = roc(y[m.values], s[m.values], w[m.values])
            summary["pt_bins"][f"{lo}-{hi}"][name] = rej_at(0.5, cs, cb)

    # eff vs pT at fixed overall mistag
    plt.figure(figsize=(7, 5))
    pts = 0.5 * (np.array([b[0] for b in PT_BINS]) + np.array([b[1] for b in PT_BINS]))
    for name in [k for k in ("bdt_nom", "part_nom", "pnet_wvsqcd", "tau21+masswin") if k in scores]:
        s = scores[name]
        bkg = s[(y == 0)]
        thr = np.quantile(bkg[np.isfinite(bkg)], 0.99)  # 1% overall mistag
        effs = [np.average((s > thr)[m.values & (y == 1)],
                           weights=w[m.values & (y == 1)])
                for lo, hi in PT_BINS
                for m in [(df.pt > lo) & (df.pt < hi)]]
        plt.plot(pts, effs, "o-", label=f"{name} (1% mistag WP)")
    plt.xlabel("AK8 jet pT [GeV]")
    plt.ylabel("signal efficiency")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(RES, "eff_vs_pt.png"), dpi=140, bbox_inches="tight")
    plt.close()

    # mass sculpting: QCD m_SD passing vs all, JSD vs eps_b
    edges = np.linspace(20, 250, 47)
    qcd = df[y == 0]
    base_hist, _ = np.histogram(qcd.sdmass, bins=edges, weights=w[y == 0])
    plt.figure(figsize=(7, 5))
    plt.stairs(base_hist / base_hist.sum(), edges, color="k", label="QCD inclusive")
    jsd_scan = {}
    for name in [k for k in ("bdt_nom", "bdt_md", "part_nom", "part_md",
                             "n2ddt(5%)", "pnet_wvsqcd", "pnet_xqqvsqcd") if k in scores]:
        s = scores[name][y == 0]
        thr = np.quantile(s[np.isfinite(s)], 0.95)   # eps_b = 5%
        hp, _ = np.histogram(qcd.sdmass[s > thr], bins=edges, weights=w[y == 0][s > thr])
        if hp.sum() > 0:
            plt.stairs(hp / hp.sum(), edges, label=f"{name} pass (5% eps_b)")
        jsd_scan[name] = []
        for eb in (0.20, 0.10, 0.05, 0.02, 0.01):
            thr = np.quantile(s[np.isfinite(s)], 1 - eb)
            hp, _ = np.histogram(qcd.sdmass[s > thr], bins=edges, weights=w[y == 0][s > thr])
            jsd_scan[name].append(jsd(hp + 1e-12, base_hist + 1e-12) if hp.sum() > 0 else np.nan)
    plt.xlabel("QCD jet m_SD [GeV]")
    plt.ylabel("normalized")
    plt.legend(fontsize=7)
    plt.savefig(os.path.join(RES, "sculpting.png"), dpi=140, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7, 5))
    for name, vals in jsd_scan.items():
        plt.plot([0.20, 0.10, 0.05, 0.02, 0.01], vals, "o-", label=name)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("background efficiency")
    plt.ylabel("JSD( QCD m_SD pass || inclusive )")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(RES, "jsd_vs_effb.png"), dpi=140, bbox_inches="tight")
    plt.close()
    summary["jsd_vs_effb"] = jsd_scan

    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(open(os.path.join(RES, "rejection_table.md")).read())
    print("results in", RES)


if __name__ == "__main__":
    main()
