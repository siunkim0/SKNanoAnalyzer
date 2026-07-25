#!/usr/bin/env python3
"""Register Wtag-v2 samples (TTto4Q + inclusive QCD_PT bins) in the SKNano
Run3_v13_Run2_v9 sample DB for 2022 and 2022EE.

Writes data/<db>/<era>/Sample/ForSNU/<alias>.json (path list) and inserts the
alias into CommonSampleInfo.json with exact nmc/sumW from the files' Runs
trees.  sumSign: pythia QCD has all-positive weights (sumSign = nmc); for
powheg TTto4Q it is summed from Events.genWeight (slow, one pass).

Run inside an env with ROOT (LCG el9 view).  Idempotent — re-running
overwrites the same entries.
"""
import glob
import json
import os
import sys

import ROOT

ROOT.EnableImplicitMT(8)

REPO = "/data6/Users/snuintern2/wtag/SKNanoAnalyzer"
DB = "Run3_v13_Run2_v9"
STORE = "/gv0/DATA/SKNano/Run3NanoAODv13p1"
ERAS = ["2022", "2022EE"]

# alias -> (PD dir name, xsec [pb], ref, needs_sign_loop)
# QCD xsecs copied from the group's Run3_v12 DB (XSDB, 13.6 TeV).
# TTto4Q: sigma_tt(NNLO) * B(W->qq)^2 = 923.6 * 0.6741^2.
SAMPLES = {
    "TTJJ_powheg":      ("TTto4Q_TuneCP5_13p6TeV_powheg-pythia8", 419.68,
                         "https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO", True),
    "QCD_Pt-30to50":    ("QCD_PT-30to50_TuneCP5_13p6TeV_pythia8",   113300000.0, "XSDB", False),
    "QCD_Pt-50to80":    ("QCD_PT-50to80_TuneCP5_13p6TeV_pythia8",    16760000.0, "XSDB", False),
    "QCD_Pt-80to120":   ("QCD_PT-80to120_TuneCP5_13p6TeV_pythia8",   2534000.0, "XSDB", False),
    "QCD_Pt-120to170":  ("QCD_PT-120to170_TuneCP5_13p6TeV_pythia8",   445800.0,  "XSDB", False),
    "QCD_Pt-170to300":  ("QCD_PT-170to300_TuneCP5_13p6TeV_pythia8",  113700.0,  "XSDB", False),
    "QCD_Pt-300to470":  ("QCD_PT-300to470_TuneCP5_13p6TeV_pythia8",  7589.0,    "XSDB", False),
    "QCD_Pt-470to600":  ("QCD_PT-470to600_TuneCP5_13p6TeV_pythia8",  626.4,     "XSDB", False),
    "QCD_Pt-600to800":  ("QCD_PT-600to800_TuneCP5_13p6TeV_pythia8",  178.6,     "XSDB", False),
    "QCD_Pt-800to1000": ("QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8", 30.57,     "XSDB", False),
    "QCD_Pt-1000to1400": ("QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8", 8.92,    "XSDB", False),
    "QCD_Pt-1400to1800": ("QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8", 0.8103,  "XSDB", False),
    "QCD_Pt-1800to2400": ("QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8", 0.1148,  "XSDB", False),
    "QCD_Pt-2400to3200": ("QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8", 0.007542, "XSDB", False),
}


def natural_key(path):
    base = os.path.basename(path)
    num = "".join(c for c in base if c.isdigit())
    return (int(num) if num else 0, base)


def register(era, alias, pd, xsec, ref, sign_loop):
    files = sorted(glob.glob(f"{STORE}/{era}/{pd}/*/*/NANOAOD_*.root"), key=natural_key)
    if not files:
        print(f"[{era}] {alias}: NO FILES under {STORE}/{era}/{pd} — skipped")
        return None

    runs = ROOT.RDataFrame("Runs", files)
    nmc = runs.Sum("genEventCount").GetValue()
    sumw = runs.Sum("genEventSumw").GetValue()

    if sign_loop:
        ev = ROOT.RDataFrame("Events", files)
        sumsign = ev.Define("s", "genWeight > 0 ? 1. : (genWeight < 0 ? -1. : 0.)").Sum("s").GetValue()
    else:
        sumsign = float(nmc)

    forsnu = {
        "name": alias, "isMC": 1, "PD": pd, "xsec": xsec,
        "sumsign": -1, "sumW": -1, "nmc": -1, "path": files,
    }
    forsnu_path = f"{REPO}/data/{DB}/{era}/Sample/ForSNU/{alias}.json"
    with open(forsnu_path, "w") as f:
        json.dump(forsnu, f, indent=4)

    common_path = f"{REPO}/data/{DB}/{era}/Sample/CommonSampleInfo.json"
    with open(common_path) as f:
        common = json.load(f)
    common[alias] = {
        "isMC": 1, "PD": pd, "xsec": xsec,
        "nmc": float(nmc), "sumsign": float(sumsign), "sumW": float(sumw),
        "ref": ref,
    }
    with open(common_path, "w") as f:
        json.dump(common, f, indent=4)

    print(f"[{era}] {alias}: {len(files)} files  nmc={nmc:.0f}  sumsign={sumsign:.0f}  sumW={sumw:.6g}")
    return nmc


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for era in ERAS:
        for alias, (pd, xsec, ref, sign_loop) in SAMPLES.items():
            if only and only not in alias:
                continue
            register(era, alias, pd, xsec, ref, sign_loop)
