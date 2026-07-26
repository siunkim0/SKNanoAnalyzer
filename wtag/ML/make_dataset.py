#!/usr/bin/env python3
"""Build the training dataset directly from the per-sample v2 outputs on /gv0
(870M stored jets — far more than training needs, so samples are skimmed with
per-sample downsampling; the keep-probability is divided back into `weight`
so cross-section normalization stays exact).

Stage 1 (PyROOT RDataFrame, multithreaded): selection + downsampling
         → data/skim/<era>_<sample>.root
Stage 2 (uproot/pandas): derived features, class labels, splits, flattening
         weights → data/dataset.parquet + data/meta.json

Run on tamsa2 with the LCG el9 view:
  python3 make_dataset.py            # both stages
  python3 make_dataset.py skim|pack  # one stage
"""
import glob
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = "/gv0/Users/snuintern2/SKNanoOutput/Wtag"
ERAS = ["2022_v2", "2022EE_v2"]
SKIM = os.path.join(BASE, "data", "skim")
OUTDIR = os.path.join(BASE, "data")

PT_MIN, PT_MAX = 200.0, 1200.0
MSD_MIN, MSD_MAX = 20.0, 250.0
SELECTION = (f"pt > {PT_MIN} && pt < {PT_MAX} && "
             f"sdmass > {MSD_MIN} && sdmass < {MSD_MAX}")

SIG_SAMPLES = ["TTLJ_powheg", "TTJJ_powheg", "WW_pythia", "WZ_pythia"]
BKG_PREFIX = "QCD_Pt-"
# cap on *selected* jets per (sample, era) after cuts
SIG_CAP = 3_000_000        # label==1 is a small fraction; caps mostly the spectators
QCD_CAP = 800_000          # 18 sample-eras -> <= 14.4M background

MASS_FEATURES = ["sdmass", "mass", "mreg", "sj1_mass", "sj2_mass", "sj_dr"]
FEATURES = [
    "sdmass", "mass", "mreg",
    # tau31/tau41 are the two best non-tagger variables measured in v2c
    # (rejection 18.0 / 15.9 at eps_S=0.5, vs 14.4 for tau21). They are not
    # reachable from the adjacent ratios: tau31 = tau21 * tau32 is a product,
    # which a tree ensemble cannot form, and no raw tau_N is in this list.
    "tau21", "tau32", "tau43", "tau31", "tau41",
    "n2b1", "n3b1", "lsf3",
    "nconst", "chmult", "nemult",
    "chhef", "nehef", "chemef", "neemef", "muef",
    "sj1_mass", "sj2_mass", "sj1_btag", "sj2_btag", "sj_dr", "sj_z",
    "nsv", "sv_summass", "sv_maxdlensig", "sv_sumntracks",
]
TOKEN_COLS = ["sj1_pt", "sj2_pt"] + [
    f"sv{k}_{v}" for k in (1, 2, 3) for v in ("pt", "mass", "dlensig", "ntracks")]
# Benchmarks, deliberately NOT features -- a model given these would just learn
# to copy them. pnet_qcd/xcc/xbb are here so the *combined* mass-decorrelated
# discriminant (Xqq+Xcc+Xbb)/(...+QCD) can be reconstructed: the stored ratios
# are r = X/(X+QCD), so X = QCD*r/(1-r) needs the raw QCD probability.
EXTRAS = ["pt", "eta", "weight", "label", "pnet_wvsqcd", "pnet_tvsqcd",
          "pnet_xqqvsqcd", "pnet_xggvsqcd", "pnet_qcd", "pnet_xccvsqcd",
          "pnet_xbbvsqcd", "gen_wpt", "gen_drqq", "gen_hasb"]


def skim():
    import ROOT
    ROOT.EnableImplicitMT(8)
    os.makedirs(SKIM, exist_ok=True)
    for era in ERAS:
        for path in sorted(glob.glob(f"{SRC}/{era}/*.root")):
            sample = os.path.basename(path).replace(".root", "")
            is_qcd = sample.startswith(BKG_PREFIX)
            if not is_qcd and sample not in SIG_SAMPLES:
                continue
            cap = QCD_CAP if is_qcd else SIG_CAP
            out = f"{SKIM}/{era}_{sample}.root"
            if os.path.exists(out):
                print(f"[{era}] {sample:22s} skim exists — skipped")
                continue
            df = ROOT.RDataFrame("jets", path)
            # exact selected count (reads only pt/sdmass) sets the keep-probability
            n_sel = df.Filter(SELECTION).Count().GetValue()
            keep = min(1.0, cap / max(n_sel, 1))
            d = (df.Filter(SELECTION)
                   .Define("rnd", "gRandom->Rndm()")
                   .Filter(f"rnd < {keep}")
                   .Define("weight_c", f"weight / {keep}"))
            cols = ROOT.std.vector("string")()
            for c in (FEATURES + TOKEN_COLS + EXTRAS):
                if c == "weight":
                    continue
                if c in ("tau21", "tau32", "tau43", "tau31", "tau41", "mreg"):
                    continue  # derived later
                cols.push_back(c)
            for c in ("tau1", "tau2", "tau3", "tau4", "masscorr", "weight_c"):
                cols.push_back(c)
            d.Snapshot("jets", out, cols)
            n_out = ROOT.RDataFrame("jets", out).Count().GetValue()
            print(f"[{era}] {sample:22s} sel={n_sel:>11,} keep={keep:.4f} -> {n_out:,}",
                  flush=True)


def pack():
    import pandas as pd
    import uproot

    frames = []
    for path in sorted(glob.glob(f"{SKIM}/*.root")):
        name = os.path.basename(path).replace(".root", "")
        era, sample = name.split("_v2_", 1)   # "<era>_v2_<sample>"
        with uproot.open(path + ":jets") as t:
            df = t.arrays(library="pd")
        df = df.rename(columns={"weight_c": "weight"})
        is_qcd = sample.startswith(BKG_PREFIX)
        if is_qcd:
            df = df[df.label == 0].copy()
            df["y"] = 0
        else:
            df = df[df.label.isin([1, 2])].copy()
            # y: 1 = clean W, -1 = t->bqq spectator (never train/eval background)
            df["y"] = np.where(df.label == 1, 1, -1)
        df["proc"], df["era"] = sample, era
        frames.append(df)
        print(f"{name:34s} kept {len(df):>9,}")
    data = pd.concat(frames, ignore_index=True)

    # drop rare corrupted rows (NaN in any feature; O(10) per 17M jets)
    n0 = len(data)
    data = data.dropna(subset=[c for c in FEATURES + TOKEN_COLS if c in data]).copy()
    if len(data) < n0:
        print(f"dropped {n0 - len(data)} NaN rows")
    # sanitize sentinels (N2/N3 defined only for raw pT>250; PNet uses -10)
    for c in ("n2b1", "n3b1", "lsf3"):
        data[c] = np.where(data[c] > -1.0, data[c], -1.0)
    # The vs-QCD ratios use -10 for "not evaluated"; fold that onto -1 so a
    # single test (< 0) means undefined. pnet_qcd is a raw probability with no
    # sentinel, so it is left alone -- clamping it would corrupt the X = QCD*r/(1-r)
    # inversion that the combined MD benchmark depends on.
    for c in ("pnet_wvsqcd", "pnet_tvsqcd", "pnet_xqqvsqcd", "pnet_xggvsqcd",
              "pnet_xccvsqcd", "pnet_xbbvsqcd"):
        data[c] = np.where(data[c] >= 0.0, data[c], -1.0)
    # derived features
    data["mreg"] = data.masscorr * data.mass
    for num, den, nm in ((2, 1, "tau21"), (3, 2, "tau32"), (4, 3, "tau43"),
                         (3, 1, "tau31"), (4, 1, "tau41")):
        d = data[f"tau{den}"]
        data[nm] = np.where(d > 0, data[f"tau{num}"] / d, -1.0)
    # rare NaNs (a handful of jets with NaN subjet-btag / SV kinematics) -> sentinel
    # (after the derived columns exist)
    data[FEATURES + TOKEN_COLS] = data[FEATURES + TOKEN_COLS].fillna(-1.0)

    # train/val/test split 60/20/20, reproducible
    rng = np.random.default_rng(20260717)
    r = rng.random(len(data))
    data["split"] = np.select([r < 0.6, r < 0.8], ["train", "val"], default="test")

    # per-class flattening weights (training rows: label 1 vs QCD label 0)
    trainable = data.y >= 0
    pt_edges = np.linspace(PT_MIN, PT_MAX, 51)
    msd_edges = np.linspace(MSD_MIN, MSD_MAX, 24)
    data["w_flatpt"] = 0.0
    data["w_flatptm"] = 0.0
    for y in (0, 1):
        m = (data.y == y) & trainable
        pts = data.loc[m, "pt"].values
        hist, _ = np.histogram(pts, bins=pt_edges)
        ipt = np.clip(np.digitize(pts, pt_edges) - 1, 0, 49)
        w = np.where(hist[ipt] > 0, 1.0 / hist[ipt], 0.0)
        data.loc[m, "w_flatpt"] = w * m.sum() / w.sum()
        im = np.clip(np.digitize(data.loc[m, "sdmass"], msd_edges) - 1, 0, 22)
        comb = ipt * 23 + im
        h2 = np.bincount(comb, minlength=50 * 23).astype(float)
        w2 = np.where(h2[comb] > 0, 1.0 / h2[comb], 0.0)
        data.loc[m, "w_flatptm"] = w2 * m.sum() / w2.sum()

    cols = list(dict.fromkeys(
        FEATURES + TOKEN_COLS + EXTRAS
        + ["proc", "era", "y", "split", "w_flatpt", "w_flatptm"]))
    out = os.path.join(OUTDIR, "dataset.parquet")
    data[cols].to_parquet(out + ".tmp", row_group_size=1_000_000)
    os.replace(out + ".tmp", out)   # atomic — readers never see a partial file

    meta = {
        "eras": ERAS,
        "features": FEATURES,
        "features_md": [f for f in FEATURES if f not in MASS_FEATURES],
        "selection": {"pt": [PT_MIN, PT_MAX], "sdmass": [MSD_MIN, MSD_MAX]},
        "n": {"sig_label1": int(((data.y == 1) & (data.label == 1)).sum()),
              "spectator_label2": int((data.label == 2).sum()),
              "qcd": int((data.y == 0).sum())},
    }
    with open(os.path.join(OUTDIR, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print("wrote dataset.parquet:", meta["n"])


if __name__ == "__main__":
    stages = sys.argv[1:] or ["skim", "pack"]
    if "skim" in stages:
        skim()
    if "pack" in stages:
        pack()
