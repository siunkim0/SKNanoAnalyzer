#!/usr/bin/env python3
"""Wtag v2c post-processing: sum the branch-comparison histograms into two files.

Only the `Branch/` and `BranchNoMSD/` directories are read. The condor output is
~76 GB, almost all of it the 20M-entry `jets` tree, and plain `hadd` would drag
that through memory (it is what pushes the hadd nodes past their 8 GB cgroup) for
a result the branch study never looks at. Reading just the histograms turns the
whole merge into a few hundred KB and a few seconds.

Class assignment is not by label alone: `label == 0` also occurs in the signal
samples, but the ML dataset deliberately takes the QCD class from the QCD samples
only, so the split has to happen at the file level here too.

  Branch/       m_SD in [20,250] -- the deployment regime
  BranchNoMSD/  no mass window   -- the unconditional survey

Run: source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
     python3 merge_v2c.py
"""
import glob
import os
import re
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)          # we own the histograms, not gDirectory
ROOT.EnableImplicitMT(8)              # the tau pass reads ~30 GB of `jets` tree

SRC = "/gv0/Users/snuintern2/SKNanoOutput/Wtag"
DST = "/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/v2c"
ERAS = ("2022_v2c", "2022EE_v2c")
DIRS = ("Branch", "BranchNoMSD")

SIG = ["TTLJ_powheg", "TTJJ_powheg", "WW_pythia", "WZ_pythia"]
QCD = ["QCD_Pt-%s" % b for b in
       ("170to300", "300to470", "470to600", "600to800", "800to1000",
        "1000to1400", "1400to1800", "1800to2400", "2400to3200")]

RUNLOG = "/gv0/Users/snuintern2/SKNanoRunlog"


def rescue_map():
    """{(era, sample): [per-job files]} for samples whose hadd node died.

    The hadd nodes are hardcoded to request_memory 8192 (python/SKNano.py) and a
    few of them need ~9.8 GB, because hadd drags the 20M-entry `jets` tree
    through memory. When one dies its output file is left truncated -- readable,
    silently short, exactly the kind of thing that turns into a wrong number
    rather than an error.

    The tell is exact rather than heuristic: hadd.sh removes output/hists_*.root
    only after hadd returns, so per-job files still sitting there mean that
    sample's hadd did not finish. Merge from those instead, which also skips the
    tree entirely and so cannot hit the memory limit in the first place.
    """
    out = {}
    for hs in glob.glob(os.path.join(RUNLOG, "*_Wtag", "*", "*", "hadd.sh")):
        jobs = sorted(glob.glob(os.path.join(os.path.dirname(hs), "output",
                                             "hists_*.root")))
        if not jobs:
            continue                      # hadd finished and cleaned up
        m = re.search(r"([^/\s]+)/([^/\s]+)\.root\s", open(hs).read())
        if not m:
            continue
        era_tag, sample = m.group(1), m.group(2)
        if era_tag in ERAS:
            out[(era_tag, sample)] = jobs
    return out


def add_file(path, acc):
    """Add one file's Branch/ + BranchNoMSD/ histograms into acc. False if unusable."""
    # A file still being written (or truncated by a killed hadd) is present but
    # unreadable, and TFile.Open *raises* rather than returning null.
    try:
        f = ROOT.TFile.Open(path)
    except OSError:
        return False
    if not f or f.IsZombie():
        return False
    for d in DIRS:
        tdir = f.Get(d)
        if not tdir:
            continue
        for key in tdir.GetListOfKeys():
            h = key.ReadObj()
            if not h.InheritsFrom("TH1"):
                continue
            full = "%s/%s" % (d, key.GetName())
            if full in acc:
                acc[full].Add(h)
            else:
                acc[full] = h.Clone()
    f.Close()
    return True


def file_list(samples, rescue):
    """[(path, is_rescue)] for every era x sample, per-job files where hadd died."""
    out = []
    for era in ERAS:
        for s in samples:
            if (era, s) in rescue:             # hadd died; use the per-job files
                out += [(p, True) for p in rescue[(era, s)]]
            else:
                out.append((os.path.join(SRC, era, s + ".root"), False))
    return out


def collect(samples, out_name, rescue):
    """Sum Branch/ and BranchNoMSD/ over every era x sample into one file."""
    acc = {}                                   # "Dir/hist" -> TH1
    missing = []
    for era in ERAS:
        for s in samples:
            if (era, s) in rescue:             # hadd died; use the per-job files
                jobs = rescue[(era, s)]
                ok = sum(add_file(p, acc) for p in jobs)
                print("   [rescue] %s/%s: hadd did not finish, summed %d/%d "
                      "per-job files" % (era, s, ok, len(jobs)))
                if ok < len(jobs):
                    missing.append("%s/%s (%d/%d job files readable)"
                                   % (era, s, ok, len(jobs)))
                continue
            path = os.path.join(SRC, era, s + ".root")
            if not os.path.exists(path):
                missing.append("%s/%s" % (era, s))
                continue
            if not add_file(path, acc):
                missing.append("%s/%s (unreadable)" % (era, s))
                continue

    out = ROOT.TFile(os.path.join(DST, out_name), "RECREATE")
    for d in DIRS:
        out.mkdir(d)
    for full, h in sorted(acc.items()):
        d, name = full.split("/", 1)
        out.cd(d)
        h.Write(name)
    out.Close()
    return acc, missing


# tau_i/tau_j pairs built here rather than in the analyzer. NEW are the three
# the v2c jobs never filled; CHECK are the three that exist both ways and so
# measure how faithful this offline route is.
TAU_NEW = ((3, 1), (4, 1), (4, 2))
TAU_CHECK = ((2, 1), (3, 2), (4, 3))


def sep(hs, hb):
    """TMVA separation of two histograms; same definition as wtag_plots.py."""
    a, b = hs.Clone(), hb.Clone()
    a.Scale(1. / a.Integral()); b.Scale(1. / b.Integral())
    return 0.5 * sum((a.GetBinContent(i) - b.GetBinContent(i)) ** 2
                     / (a.GetBinContent(i) + b.GetBinContent(i))
                     for i in range(1, a.GetNbinsX() + 1)
                     if a.GetBinContent(i) + b.GetBinContent(i) > 0)


def tau_class(samples, label, rescue):
    """One class's tau_i/tau_j histograms, {(dir, name): TH1}, from the jets tree."""
    files = [p for p, _ in file_list(samples, rescue) if os.path.exists(p)]
    df = ROOT.RDataFrame("jets", files).Filter("label == %d" % label) \
                                       .Filter("pt > 200. && pt < 1200.")
    sel = {"BranchNoMSD": df, "Branch": df.Filter("sdmass > 20. && sdmass < 250.")}

    booked = {}
    for i, j in TAU_NEW + TAU_CHECK:
        name = "tau%d%d" % (i, j)
        for d, node in sel.items():
            # -1 when the denominator is undefined, then bucket out-of-range into
            # the padding bins, exactly as the analyzer's S() lambda does.
            v, x = "v_%s_%s" % (d, name), "x_%s_%s" % (d, name)
            n = node.Define(v, "tau%d > 0.f ? tau%d / tau%d : -1.f" % (j, i, j)) \
                    .Define(x, "{0} < -1.f ? -1.05f : ({0} > 1.f ? 1.05f : {0})"
                               .format(v))
            booked[(d, name)] = n.Histo1D((d + name + str(label), "",
                                           220, -1.1, 1.1), x, "weight")
    return {k: h.GetValue().Clone() for k, h in booked.items()}


def tau_ratios(rescue):
    """Append tau31/tau41/tau42 to the merged files, built from the `jets` tree.

    The analyzer filled only the adjacent ratios (tau21/32/43), but tau1..tau4
    are each written to the tree, so the skip-a-step ratios are recoverable
    without resubmitting 26 DAGs. Same axis and same two selections as the
    analyzer's S() lambda, so the numbers drop straight into the same ranking.

    One caveat, measured rather than assumed: the tree stores floats truncated
    to ~11 mantissa bits (p_T lands on 0.125 GeV steps), so jets sitting within
    a rounding step of a window edge fall on the other side of it -- 0.02% of W
    and 0.03% of QCD. tau21/32/43 are rebuilt alongside and their S compared
    against the analyzer's own histograms, so the cost of this route is printed
    rather than assumed; a diff above a few times 0.001 would mean the offline
    numbers are not good enough to put in the same table as the rest.
    """
    tree = {}
    for samples, label, cls in ((SIG, 1, "W_"), (QCD, 0, "QCD_")):
        tree[cls] = tau_class(samples, label, rescue)

    print("   fidelity check against the analyzer-filled adjacent ratios:")
    worst = 0.
    for d in DIRS:
        for i, j in TAU_CHECK:
            name = "tau%d%d" % (i, j)
            ref = {}
            for cls, fn in (("W_", "SIG.root"), ("QCD_", "QCD.root")):
                f = ROOT.TFile.Open(os.path.join(DST, fn))
                ref[cls] = f.Get("%s/%s%s" % (d, cls, name))
                f.Close()
            if not all(ref.values()):
                continue
            a = sep(ref["W_"], ref["QCD_"])
            b = sep(tree["W_"][(d, name)], tree["QCD_"][(d, name)])
            worst = max(worst, abs(b - a))
            print("      %-12s %-6s analyzer S=%.4f   tree S=%.4f   diff %+.4f"
                  % (d, name, a, b, b - a))
    print("   worst |diff| = %.4f%s" % (worst, "" if worst < 0.005 else
                                        "   <-- TOO LARGE, do not quote these"))

    for cls, fn in (("W_", "SIG.root"), ("QCD_", "QCD.root")):
        out = ROOT.TFile(os.path.join(DST, fn), "UPDATE")
        for i, j in TAU_NEW:
            name = "tau%d%d" % (i, j)
            for d in DIRS:
                out.cd(d)
                tree[cls][(d, name)].Write(cls + name)
        out.Close()
    print("   wrote %d new histograms (%s) x 2 classes x %d sets"
          % (len(TAU_NEW), ", ".join("tau%d%d" % p for p in TAU_NEW), len(DIRS)))


def main():
    os.makedirs(DST, exist_ok=True)
    rescue = rescue_map()
    if rescue:
        print("%d sample(s) with an unfinished hadd, merging from per-job files:"
              % len(rescue))
        for era, s in sorted(rescue):
            print("   %s/%s" % (era, s))
        print()
    bad = []
    for samples, out_name, tag in ((SIG, "SIG.root", "W (ttbar+diboson)"),
                                   (QCD, "QCD.root", "QCD")):
        acc, missing = collect(samples, out_name, rescue)
        n = {d: sum(1 for k in acc if k.startswith(d + "/")) for d in DIRS}
        print("%-10s -> %s   %s"
              % (tag, out_name, "  ".join("%s: %d hists" % (d, n[d]) for d in DIRS)))
        if missing:
            bad += missing
            print("   MISSING (%d): %s" % (len(missing), ", ".join(missing)))
        if n["Branch"] != n["BranchNoMSD"]:
            print("   WARN: the two sets disagree in size -- expected identical "
                  "variable lists, got %d vs %d" % (n["Branch"], n["BranchNoMSD"]))

    if bad:
        print("\nSkipping the tau-ratio pass -- it reads the same files.")
    else:
        print("\n=== tau31 / tau41 / tau42 from the jets tree ===")
        tau_ratios(rescue)

    print("\n=== merged into %s ===" % DST)
    for fn in sorted(os.listdir(DST)):
        print("   %-12s %8.1f KB"
              % (fn, os.path.getsize(os.path.join(DST, fn)) / 1024.))
    if bad:
        print("\nINCOMPLETE: %d sample(s) missing -- the merge is not "
              "trustworthy until they land." % len(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
