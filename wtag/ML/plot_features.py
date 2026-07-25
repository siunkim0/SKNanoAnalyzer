#!/usr/bin/env python3
"""Per-input-variable separation check: merged W jet (signal, y==1) vs normal
QCD jet (background, y==0), for every tagger input in meta.json.

Answers "is each variable different enough to be worth tagging on?" two ways:
  * physics (xsec) weights  -> the honest, as-observed separation
  * flat-pT weights         -> intrinsic power with the pT-spectrum advantage
                               removed (both classes share a flat pT spectrum,
                               matching how the tagger is trained)
For each feature both are drawn, and the TMVA separation S in [0,1] is printed
in every panel (0 = identical shapes, 1 = no overlap).

Run on tamsa2 with the LCG el9 view (pandas/matplotlib):
  source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
  python3 plot_features.py

Outputs:
  results/physics/<feature>.png    one file per input, physics-weighted, sig vs bkg
  results/flat_pt/<feature>.png    one file per input, flat-pT-weighted, sig vs bkg
  results/separation_ranking.png   S per feature, physics vs flat, sorted
  results/separation_table.md      the same numbers as a table
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, "results")
PHYSDIR = os.path.join(RES, "physics")
FLATDIR = os.path.join(RES, "flat_pt")
os.makedirs(PHYSDIR, exist_ok=True)
os.makedirs(FLATDIR, exist_ok=True)

NBINS = 50
PRANGE = (0.5, 99.5)   # percentile clip for the histogram range (kills outliers)

# pretty axis labels; anything missing falls back to the raw column name
LABELS = {
    "sdmass": "m_SD [GeV]", "mass": "m [GeV]", "mreg": "m_reg [GeV]",
    "tau21": "tau_21", "tau32": "tau_32", "tau43": "tau_43",
    "n2b1": "N2^(1)", "n3b1": "N3^(1)", "lsf3": "LSF_3",
    "nconst": "N_constituents", "chmult": "N_charged", "nemult": "N_neutral",
    "chhef": "charged had. E frac", "nehef": "neutral had. E frac",
    "chemef": "charged em E frac", "neemef": "neutral em E frac",
    "muef": "muon E frac",
    "sj1_mass": "subjet1 mass [GeV]", "sj2_mass": "subjet2 mass [GeV]",
    "sj1_btag": "subjet1 b-tag", "sj2_btag": "subjet2 b-tag",
    "sj_dr": "dR(subjets)", "sj_z": "subjet z",
    "nsv": "N_SV", "sv_summass": "sum SV mass [GeV]",
    "sv_maxdlensig": "max SV d_len sig", "sv_sumntracks": "sum SV N_tracks",
    "pt": "AK8 pT [GeV]", "eta": "AK8 eta",
}
# context kinematics: shown for reference, not tagger inputs
EXTRA = ["pt", "eta"]

SIG_COL, BKG_COL = "#c0392b", "#2c6fbb"   # red = merged W, blue = QCD


def separation(s, b):
    """TMVA separation of two *normalized* histograms; 0=identical, 1=disjoint."""
    den = s + b
    m = den > 0
    return 0.5 * np.sum((s[m] - b[m]) ** 2 / den[m])


def edges_for(x):
    # nan-robust: a handful of NaN rows (~1e-6) would otherwise poison percentile
    lo, hi = np.nanpercentile(x, PRANGE)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x)) + 1e-6
    return np.linspace(lo, hi, NBINS + 1)


def norm_hist(x, w, edges):
    h, _ = np.histogram(x, bins=edges, weights=w)
    tot = h.sum()
    return h / tot if tot > 0 else h


def plot_one(f, xs, ws, xb, wb, edges, sub, outdir):
    """One PNG, one single-panel plot: signal vs background for feature f."""
    hs = norm_hist(xs, ws, edges)
    hb = norm_hist(xb, wb, edges)
    S = separation(hs, hb)
    c = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.step(c, hb, where="mid", color=BKG_COL, lw=1.6, label="QCD (normal)")
    ax.step(c, hs, where="mid", color=SIG_COL, lw=1.6, label="merged W")
    ax.fill_between(c, hb, step="mid", color=BKG_COL, alpha=0.15)
    ax.fill_between(c, hs, step="mid", color=SIG_COL, alpha=0.15)
    ax.set_xlabel(LABELS.get(f, f))
    ax.set_ylabel("normalized")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", frameon=False)
    ax.text(0.03, 0.95, f"S = {S:.3f}", transform=ax.transAxes,
            fontsize=11, va="top", ha="left", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.85))
    extra = "   [kinematic]" if f in EXTRA else ""
    ax.set_title(f"merged W vs normal QCD jet - {sub}{extra}", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"{f}.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    return S


def plot_feature(f, data, sig, bkg):
    """Two single-panel PNGs for feature f: physics/ and flat_pt/."""
    x = data[f]
    edges = edges_for(x[sig | bkg])
    s_phys = plot_one(f, x[sig], data["weight"][sig], x[bkg], data["weight"][bkg],
                      edges, "physics (xsec) weights", PHYSDIR)
    s_flat = plot_one(f, x[sig], data["w_flatpt"][sig], x[bkg], data["w_flatpt"][bkg],
                      edges, "flat-pT weights", FLATDIR)
    return s_phys, s_flat


def ranking(sep_phys, sep_flat, feats):
    order = sorted(feats, key=lambda f: max(sep_phys[f], sep_flat[f]))
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(order) + 1.2))
    ax.barh(y + 0.19, [sep_phys[f] for f in order], height=0.36,
            color="#7f8c8d", label="physics weights")
    ax.barh(y - 0.19, [sep_flat[f] for f in order], height=0.36,
            color="#e67e22", label="flat-pT weights")
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(f, f) for f in order], fontsize=8)
    ax.set_xlabel("separation S  (0 = identical, 1 = disjoint)")
    ax.set_title("merged W vs normal QCD jet: per-variable separation")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    ax.set_xlim(0, None)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "separation_ranking.png"),
                dpi=140, bbox_inches="tight")
    plt.close(fig)

    order_desc = order[::-1]
    lines = ["| feature | S (physics) | S (flat-pT) |", "|---|---|---|"]
    for f in order_desc:
        lines.append(f"| {LABELS.get(f, f)} | {sep_phys[f]:.3f} | {sep_flat[f]:.3f} |")
    with open(os.path.join(RES, "separation_table.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return order_desc


def main():
    meta = json.load(open(os.path.join(BASE, "data", "meta.json")))
    feats = list(meta["features"]) + EXTRA
    cols = feats + ["y", "weight", "w_flatpt"]

    print(f"reading {len(cols)} columns ...")
    # y >= 0 keeps signal (1) and QCD (0), drops t->bqq spectators (-1)
    tbl = pq.read_table(os.path.join(BASE, "data", "dataset.parquet"),
                        columns=cols, filters=[("y", ">=", 0)])
    data = {c: tbl.column(c).to_numpy() for c in cols}
    for c in feats:
        data[c] = data[c].astype(np.float32)
    data["weight"] = np.abs(data["weight"]).astype(np.float64)   # shapes: |w|
    data["w_flatpt"] = data["w_flatpt"].astype(np.float64)
    sig = data["y"] == 1
    bkg = data["y"] == 0
    print(f"merged W (signal): {int(sig.sum()):,}   normal QCD (bkg): {int(bkg.sum()):,}")

    sep_phys, sep_flat = {}, {}
    for f in feats:
        sep_phys[f], sep_flat[f] = plot_feature(f, data, sig, bkg)
        print(f"  {LABELS.get(f, f):22s} -> {f}.png")
    order = ranking(sep_phys, sep_flat, feats)

    print("\n--- separation S (sorted, flat-pT weighted) ---")
    print(f"  {'feature':22s} {'physics':>8s} {'flat-pT':>8s}")
    for f in order:
        print(f"  {LABELS.get(f, f):22s} {sep_phys[f]:8.3f} {sep_flat[f]:8.3f}")
    print("\nresults in", RES)


if __name__ == "__main__":
    main()
