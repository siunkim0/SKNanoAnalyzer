#!/usr/bin/env python3
# Wtag v1 plots. Run: rot && python3 wtag_plots.py   (after merge_v1.sh)
#   Stage 1 (gen):  fraction dR(qq') < 0.4 / 0.8  vs  W pT;  <dRqq> & <dRjj>
#                   profiles vs W pT with the 2 m_W / pT expectation overlaid
#   Stage 2 (reco): fraction merged / resolved / overlap(AK4) / captured(AK8)
#                   vs W pT, + printed 50%-crossing thresholds
#   Stage 3 (jet):  <chMult>/<nConst>/<neMult>/<nSV>/<mass> vs jet pT for
#                   mergedW vs singleW vs lightuds vs gluon, + ChMult shape
#                   comparison in fixed jet-pT slices (pT-matched hypothesis test)
#   AK8:            softdrop mass / tau21 / ParticleNet WvsQCD / NConst
#                   capturedW vs other (shapes + profiles)
#
# Histos summed over TTLJ+WW+WZ (already xsec-weighted via MCweight). Pick one
# with SAMPLE=TTLJ|WW|WZ. Reads /data6/.../wtag/SKNanoAnalyzer/wtag/<ERA>_<TAG>/<sample>.root.
import ROOT, os
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

ERA    = os.environ.get('ERA', '2022')
TAG    = os.environ.get('TAG', 'v1')
SAMPLE = os.environ.get('SAMPLE', 'all')
INDIR  = f"/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/{ERA}_{TAG}"
OUTDIR = f"/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/Tools/Plots/{TAG}"
os.makedirs(OUTDIR, exist_ok=True)
FILES  = ["TTLJ.root", "WW.root", "WZ.root"] if SAMPLE == 'all' else [f"{SAMPLE}.root"]
WPT_REBIN = 5   # 10 GeV -> 50 GeV bins for the fraction curves

_files = []
def get(hname):
    """Sum a histogram across the sample files; returns a detached clone (or None)."""
    tot = None
    for fn in FILES:
        path = os.path.join(INDIR, fn)
        if not os.path.exists(path):
            continue
        f = ROOT.TFile.Open(path)
        _files.append(f)
        h = f.Get(hname)
        if not h:
            continue
        if tot is None:
            tot = h.Clone(hname.replace('/', '_') + "_sum")
            tot.SetDirectory(0)
        else:
            tot.Add(h)
    return tot

def frac(num, den, rebin=WPT_REBIN):
    """Binomial fraction num/den vs pT (both 1D in W pT)."""
    if not num or not den:
        return None
    n, d = num.Clone(), den.Clone()
    if rebin > 1:
        n.Rebin(rebin); d.Rebin(rebin)
    n.Divide(n, d, 1., 1., "B")
    return n

def crossing(h, level=0.5, rising=True):
    """First pT where the fraction curve crosses `level` (linear interpolation)."""
    if not h:
        return None
    prev_x, prev_y = None, None
    for i in range(1, h.GetNbinsX() + 1):
        x, y = h.GetBinCenter(i), h.GetBinContent(i)
        if h.GetBinError(i) == 0 and y == 0 and prev_y is None:
            continue                       # skip empty leading bins
        if prev_y is not None:
            up   = rising and prev_y < level <= y
            down = (not rising) and prev_y > level >= y
            if up or down:
                return prev_x + (level - prev_y) * (x - prev_x) / (y - prev_y)
        prev_x, prev_y = x, y
    return None

# ---- Stage 1 & 2: fraction-vs-WpT curves -----------------------------------
def plot_fractions():
    den = get("Gen/W_Pt_all")
    if not den:
        print("[skip] Gen/W_Pt_all missing"); return
    curves = [
        ("Gen/W_Pt_dRqq_lt04",     "#DeltaR(qq') < 0.4 (gen)",    ROOT.kBlue+1,  20),
        ("Gen/W_Pt_dRqq_lt08",     "#DeltaR(qq') < 0.8 (gen)",    ROOT.kAzure+7, 24),
        ("Reco/W_Pt_ak4merged",    "AK4 merged (1 jet)",          ROOT.kRed+1,   21),
        ("Reco/W_Pt_ak4resolved",  "AK4 resolved (2 jets)",       ROOT.kMagenta+1, 25),
        ("Reco/W_Pt_ak4overlap",   "AK4 jets overlap (#DeltaR<0.8)", ROOT.kOrange+7, 22),
        ("Reco/W_Pt_ak8captured",  "AK8 captured (1 fatjet)",     ROOT.kGreen+2, 23),
    ]
    c = ROOT.TCanvas("cfrac", "", 900, 700)
    leg = ROOT.TLegend(0.50, 0.20, 0.88, 0.47)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    frame = ROOT.TH2F("frame", ";generator W p_{T} [GeV];fraction of hadronic W",
                      10, 0., 1000., 10, 0., 1.05)
    frame.Draw()
    keep = []
    print(f"--- 50% crossings ({SAMPLE}) ---")
    for hname, label, col, mk in curves:
        h = frac(get(hname), den)
        if not h:
            continue
        h.SetLineColor(col); h.SetMarkerColor(col); h.SetMarkerStyle(mk); h.SetMarkerSize(0.8)
        h.Draw("PE SAME")
        leg.AddEntry(h, label, "PE")
        keep.append(h)
        rising = "resolved" not in hname
        x50 = crossing(h, 0.5, rising)
        print(f"  {hname:26s} 50% at W pT ~ {x50:.0f} GeV" if x50 else
              f"  {hname:26s} no 50% crossing in range")
    for x in (200., 400.):   # naive 2mW/pT boundaries
        ln = ROOT.TLine(x, 0., x, 1.05); ln.SetLineStyle(2); ln.SetLineColor(ROOT.kGray+2)
        ln.Draw(); keep.append(ln)
    leg.Draw()
    c.SaveAs(f"{OUTDIR}/efficiency_vs_WpT_{SAMPLE}.pdf")

# ---- Stage 2: exclusive AK4 topology, stacked so fractions sum to 1 --------
def plot_fraction_stack():
    """Every hadronic W falls in exactly one category: merged (1 AK4 jet) /
    resolved (2 AK4 jets) / not both matched -> stack sums to 1 per pT bin."""
    den  = get("Gen/W_Pt_all")
    both = get("Reco/W_Pt_ak4bothmatched")
    if not den or not both:
        print("[skip] fraction stack: missing histograms"); return
    unmatched = den.Clone("W_Pt_unmatched"); unmatched.SetDirectory(0)
    unmatched.Add(both, -1.)
    cats = [   # bottom of the stack first
        (get("Reco/W_Pt_ak4merged"),   "AK4 merged (1 jet)",      ROOT.kRed-7),
        (get("Reco/W_Pt_ak4resolved"), "AK4 resolved (2 jets)",   ROOT.kAzure-9),
        (unmatched,                    "not both quarks matched", ROOT.kGray+1),
    ]
    d = den.Clone(); d.Rebin(WPT_REBIN)
    stack = ROOT.THStack("stack", ";generator W p_{T} [GeV];fraction of hadronic W")
    keep = [stack, d]
    leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    entries = []
    for num, label, col in cats:
        if not num:
            print("[skip] fraction stack: missing category"); return
        n = num.Clone(); n.Rebin(WPT_REBIN)
        n.Divide(n, d)                       # plain ratio: stack sums to 1
        n.SetFillColor(col); n.SetLineColor(ROOT.kBlack); n.SetLineWidth(1)
        stack.Add(n, "HIST"); keep.append(n)
        entries.append((n, label))
    c = ROOT.TCanvas("cstack", "", 900, 700)
    frame = ROOT.TH2F("framestack", ";generator W p_{T} [GeV];fraction of hadronic W",
                      10, 0., 1000., 10, 0., 1.05)
    frame.Draw(); keep.append(frame)
    stack.Draw("HIST SAME")
    for n, label in reversed(entries):       # legend top-to-bottom = stack top-to-bottom
        leg.AddEntry(n, label, "F")
    for x in (200., 400.):                   # naive 2mW/pT boundaries
        ln = ROOT.TLine(x, 0., x, 1.05); ln.SetLineStyle(2); ln.SetLineColor(ROOT.kGray+3)
        ln.Draw(); keep.append(ln)
    leg.Draw()
    ROOT.gPad.RedrawAxis()
    c.SaveAs(f"{OUTDIR}/fraction_vs_WpT_{SAMPLE}.pdf")

# ---- Stage 1: <dR> profiles vs W pT with 2mW/pT overlay --------------------
def plot_dr_profiles():
    c = ROOT.TCanvas("cdr", "", 900, 700)
    leg = ROOT.TLegend(0.45, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    frame = ROOT.TH2F("frdr", ";generator W p_{T} [GeV];#LT#DeltaR#GT",
                      10, 0., 1000., 10, 0., 3.5)
    frame.Draw()
    keep = []
    for hn, label, col, mk in [("Gen/dRqq_vs_WPt",  "#DeltaR(q,q') gen",         ROOT.kBlue+1, 20),
                               ("Reco/dRjj_vs_WPt", "#DeltaR(j_{1},j_{2}) reco", ROOT.kRed+1,  21)]:
        h2 = get(hn)
        if not h2:
            continue
        prof = h2.ProfileX("px_" + hn.replace('/', '_'))
        prof.Rebin(2)
        prof.SetLineColor(col); prof.SetMarkerColor(col); prof.SetMarkerStyle(mk); prof.SetMarkerSize(0.8)
        prof.Draw("PE SAME")
        leg.AddEntry(prof, label, "PE")
        keep.append(prof)
    fexp = ROOT.TF1("fexp", "2*80.4/x", 80., 1000.)
    fexp.SetLineColor(ROOT.kGray+2); fexp.SetLineStyle(2)
    fexp.Draw("SAME")
    leg.AddEntry(fexp, "2 m_{W} / p_{T}", "L")
    for y in (0.4, 0.8):
        ln = ROOT.TLine(0., y, 1000., y); ln.SetLineStyle(3); ln.SetLineColor(ROOT.kGray+1)
        ln.Draw(); keep.append(ln)
    leg.Draw()
    c.SaveAs(f"{OUTDIR}/dR_vs_WpT_{SAMPLE}.pdf")

    # 2D map of dRqq vs W pT
    h2 = get("Gen/dRqq_vs_WPt")
    if h2:
        c2 = ROOT.TCanvas("cdr2d", "", 900, 700)
        c2.SetRightMargin(0.13); c2.SetLogz()
        h2.SetTitle(";generator W p_{T} [GeV];#DeltaR(q,q')")
        h2.Draw("COLZ")
        fexp.Draw("SAME")
        c2.SaveAs(f"{OUTDIR}/dRqq_vs_WpT_2D_{SAMPLE}.pdf")

# ---- <observable> vs jet pT for Stage 3 categories -------------------------
def plot_profile(stem, ytitle, ymax, cats):
    c = ROOT.TCanvas("c_" + stem, "", 900, 700)
    leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    frame = ROOT.TH2F("fr_" + stem, f";jet p_{{T}} [GeV];{ytitle}",
                      10, 0., 1500., 10, 0., ymax)
    frame.Draw()
    keep = []
    for suf, label, col, mk in cats:
        h2 = get(f"Jet/{stem}_vs_Pt_{suf}")
        if not h2:
            continue
        prof = h2.ProfileX(f"px_{stem}_{suf}")
        prof.Rebin(5)   # 10 GeV -> 50 GeV bins
        prof.SetLineColor(col); prof.SetMarkerColor(col); prof.SetMarkerStyle(mk); prof.SetMarkerSize(0.8)
        prof.Draw("PE SAME")
        leg.AddEntry(prof, label, "PE")
        keep.append(prof)
    leg.Draw()
    c.SaveAs(f"{OUTDIR}/{stem}_vs_Pt_{SAMPLE}.pdf")

STAGE3_CATS = [
    ("mergedW",  "merged W jet",   ROOT.kRed+1,     21),
    ("singleW",  "single-q W jet", ROOT.kBlue+1,    20),
    ("lightuds", "light u/d/s jet", ROOT.kAzure+7,  24),
    ("gluon",    "gluon jet",      ROOT.kGreen+2,   23),
]

def plot_stage3():
    plot_profile("ChMult", "#LT N_{charged} #GT per jet",      40., STAGE3_CATS)
    plot_profile("NConst", "#LT N_{constituents} #GT per jet", 60., STAGE3_CATS)
    plot_profile("NeMult", "#LT N_{neutral} #GT per jet",      30., STAGE3_CATS)
    plot_profile("NSV",    "#LT N_{SV} #GT per jet",           1.2, STAGE3_CATS)
    plot_profile("Mass",   "#LT jet mass #GT [GeV]",          120., STAGE3_CATS)

# ---- Stage 3: shape comparison in fixed jet-pT slices (pT-matched test) ----
def plot_slices(stem, xtitle, ptlos=(200., 400.), pthis=(400., 600.), xmax=60.):
    for ptlo, pthi in zip(ptlos, pthis):
        c = ROOT.TCanvas(f"cs_{stem}_{ptlo:.0f}", "", 900, 700)
        leg = ROOT.TLegend(0.55, 0.62, 0.88, 0.88)
        leg.SetBorderSize(0); leg.SetFillStyle(0)
        leg.SetHeader(f"{ptlo:.0f} < jet p_{{T}} < {pthi:.0f} GeV")
        keep, means = [], []
        for suf, label, col, mk in STAGE3_CATS:
            h2 = get(f"Jet/{stem}_vs_Pt_{suf}")
            if not h2:
                continue
            ax = h2.GetXaxis()
            h = h2.ProjectionY(f"py_{stem}_{suf}_{ptlo:.0f}",
                               ax.FindBin(ptlo + 1e-3), ax.FindBin(pthi - 1e-3))
            if h.Integral() <= 0:
                continue
            means.append((label, h.GetMean(), h.GetMeanError()))
            h.Scale(1. / h.Integral())
            h.SetLineColor(col); h.SetLineWidth(2)
            h.SetTitle(f";{xtitle};a.u.")
            h.GetXaxis().SetRangeUser(0., xmax)
            h.Draw("HIST SAME" if keep else "HIST")
            leg.AddEntry(h, label, "L")
            keep.append(h)
        if keep:
            m = max(h.GetMaximum() for h in keep)
            for h in keep:
                h.SetMaximum(1.3 * m)
        leg.Draw()
        c.SaveAs(f"{OUTDIR}/{stem}_slice{ptlo:.0f}to{pthi:.0f}_{SAMPLE}.pdf")
        print(f"--- <{stem}> in {ptlo:.0f}-{pthi:.0f} GeV ({SAMPLE}) ---")
        for label, mean, err in means:
            print(f"  {label:16s} {mean:6.2f} +- {err:.2f}")

# ---- AK8 shapes: capturedW vs other -----------------------------------------
def plot_ak8_shape(stem, xtitle, xmax, rebin=1):
    c = ROOT.TCanvas("c8_" + stem, "", 900, 700)
    leg = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    keep = []
    for suf, label, col in [("capturedW", "captured W", ROOT.kRed+1),
                            ("other",     "other AK8",  ROOT.kGray+2)]:
        h = get(f"FatJet/{stem}_{suf}")
        if not h:
            continue
        if rebin > 1:
            h.Rebin(rebin)
        if h.Integral() > 0:
            h.Scale(1. / h.Integral())
        h.SetLineColor(col); h.SetLineWidth(2)
        h.SetTitle(f";{xtitle};a.u.")
        h.GetXaxis().SetRangeUser(0., xmax)
        h.Draw("HIST SAME" if keep else "HIST")
        leg.AddEntry(h, label, "L")
        keep.append(h)
    if keep:
        m = max(h.GetMaximum() for h in keep)
        for h in keep:
            h.SetMaximum(1.3 * m)
    leg.Draw()
    c.SaveAs(f"{OUTDIR}/ak8_{stem.lower()}_{SAMPLE}.pdf")

def plot_ak8():
    plot_ak8_shape("SDMass", "AK8 softdrop mass [GeV]",       300., rebin=2)
    plot_ak8_shape("Tau21",  "AK8 #tau_{2}/#tau_{1}",           1., rebin=2)
    plot_ak8_shape("WvsQCD", "ParticleNet W vs QCD",            1., rebin=2)
    plot_ak8_shape("NConst", "AK8 N_{constituents}",          150., rebin=2)

# ---- v2b: raw FatJet branch comparison, merged W vs normal QCD jet ---------
# One panel per NanoAOD FatJet branch, both classes area-normalized, with the
# TMVA separation S printed in the panel (same definition as ML/plot_features.py).
# The W class comes from the ttbar/diboson merge and the QCD class from the QCD
# merge: `label == 0` also exists in the signal samples, and the ML dataset
# deliberately takes the QCD class from the QCD samples only.
BR_INDIR = os.environ.get('BR_INDIR', "/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/v2b")
BR_OUT   = "/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/Tools/Plots/v2b/branches"

# Full spelled-out names: an axis that says "charged particle multiplicity" is
# readable by someone who has never seen the NanoAOD branch called chMultiplicity.
BR_LABELS = {
    "pt":        "AK8 jet transverse momentum p_{T} [GeV]",
    "eta":       "AK8 jet pseudorapidity #eta",
    "phi":       "AK8 jet azimuthal angle #phi",
    "mass":      "AK8 jet invariant mass [GeV]",
    "msoftdrop": "soft-drop mass (PUPPI, corrected) [GeV]",
    "mreg":      "ParticleNet regressed mass [GeV]",
    "masscorr":  "ParticleNet mass-regression correction factor",
    "area":      "jet catchment area",
    "tau1": "N-subjettiness #tau_{1} (1 axis)",
    "tau2": "N-subjettiness #tau_{2} (2 axes)",
    "tau3": "N-subjettiness #tau_{3} (3 axes)",
    "tau4": "N-subjettiness #tau_{4} (4 axes)",
    "tau21": "N-subjettiness ratio #tau_{2}/#tau_{1} (2-prong)",
    "tau32": "N-subjettiness ratio #tau_{3}/#tau_{2} (3-prong)",
    "tau43": "N-subjettiness ratio #tau_{4}/#tau_{3} (4-prong)",
    "n2b1": "energy-correlation ratio N_{2} (#beta=1)",
    "n3b1": "energy-correlation ratio N_{3} (#beta=1)",
    "lsf3": "lepton subjet fraction (3 subjets)",
    "nconst": "number of particle-flow constituents",
    "chmult": "charged particle multiplicity (PUPPI-weighted)",
    "nemult": "neutral particle multiplicity (PUPPI-weighted)",
    "nbhad":  "number of b hadrons (ghost-clustered)",
    "nchad":  "number of c hadrons (ghost-clustered)",
    "chhef":  "charged hadron energy fraction",
    "nehef":  "neutral hadron energy fraction",
    "chemef": "charged electromagnetic energy fraction",
    "neemef": "neutral electromagnetic energy fraction",
    "muef":   "muon energy fraction",
    "ddbvl":      "DeepDoubleX H#rightarrowbb vs QCD (V2, mass-decorrelated)",
    "ddcvb":      "DeepDoubleX H#rightarrowcc vs H#rightarrowbb (V2)",
    "ddcvl":      "DeepDoubleX H#rightarrowcc vs QCD (V2)",
    "btag_deepb": "DeepCSV b+bb discriminator",
    "btag_hbb":   "Higgs#rightarrowbb tagger discriminator",
}
# ParticleNet score suffix -> spelled-out physics process
PNET_PROC = {
    "h4qvsqcd": "H#rightarrowVV#rightarrow4q vs QCD",
    "hbbvsqcd": "H#rightarrowbb vs QCD", "hccvsqcd": "H#rightarrowcc vs QCD",
    "wvsqcd":   "W#rightarrowqq vs QCD", "zvsqcd": "Z#rightarrowqq vs QCD",
    "tvsqcd":   "top#rightarrowbqq vs QCD",
    "qcd":      "QCD score (sum)",
    "qcd0hf":   "QCD with 0 heavy-flavour hadrons",
    "qcd1hf":   "QCD with 1 heavy-flavour hadron",
    "qcd2hf":   "QCD with 2 heavy-flavour hadrons",
    "xbbvsqcd": "X#rightarrowbb vs QCD", "xccvsqcd": "X#rightarrowcc vs QCD",
    "xqqvsqcd": "X#rightarrowqq (uds) vs QCD", "xggvsqcd": "X#rightarrowgg vs QCD",
    "xtevsqcd": "X#rightarrowe#tau_{h} vs QCD",
    "xtmvsqcd": "X#rightarrow#mu#tau_{h} vs QCD",
    "xttvsqcd": "X#rightarrow#tau_{h}#tau_{h} vs QCD",
}
def br_label(n):
    if n in BR_LABELS:
        return BR_LABELS[n]
    if n.startswith("PNetM_"):
        return "ParticleNet (with mass): " + PNET_PROC.get(n[6:], n[6:])
    if n.startswith("PNet_"):
        return "ParticleNet (mass-decorrelated): " + PNET_PROC.get(n[5:], n[5:])
    return n

def br_get(fname, hname):
    """Fetch one histogram from a merged class file; returns a detached clone."""
    path = os.path.join(BR_INDIR, fname)
    if not os.path.exists(path):
        return None
    f = ROOT.TFile.Open(path)
    _files.append(f)
    h = f.Get(hname)
    if not h:
        return None
    h = h.Clone(hname.replace('/', '_') + "_" + fname.replace('.root', ''))
    h.SetDirectory(0)
    return h

def separation(hs, hb):
    """TMVA separation of two normalized histograms; 0 = identical, 1 = disjoint."""
    S = 0.
    for i in range(1, hs.GetNbinsX() + 1):
        s, b = hs.GetBinContent(i), hb.GetBinContent(i)
        if s + b > 0:
            S += (s - b) ** 2 / (s + b)
    return 0.5 * S

def _zoom(hists):
    """X-range covering the populated bins of both classes (kills empty margins,
    e.g. the unused negative half of a [-1,1] score axis when no sentinel fires)."""
    lo, hi = None, None
    for h in hists:
        for i in range(1, h.GetNbinsX() + 1):
            if h.GetBinContent(i) > 0:
                lo = h.GetBinLowEdge(i) if lo is None else min(lo, h.GetBinLowEdge(i))
                break
        for i in range(h.GetNbinsX(), 0, -1):
            if h.GetBinContent(i) > 0:
                e = h.GetBinLowEdge(i) + h.GetBinWidth(i)
                hi = e if hi is None else max(hi, e)
                break
    return lo, hi

def plot_branches():
    os.makedirs(BR_OUT, exist_ok=True)
    sig = ROOT.TFile.Open(os.path.join(BR_INDIR, "SIG.root"))
    d = sig.Get("Branch")
    if not d:
        print(f"[skip] no Branch/ directory in {BR_INDIR}/SIG.root"); return
    names = sorted(k.GetName()[2:] for k in d.GetListOfKeys()
                   if k.GetName().startswith("W_"))
    print(f"{len(names)} branches from {BR_INDIR}")

    seps = {}
    for n in names:
        hs = br_get("SIG.root", f"Branch/W_{n}")
        hb = br_get("QCD.root", f"Branch/QCD_{n}")
        if not hs or not hb or hs.Integral() <= 0 or hb.Integral() <= 0:
            print(f"  [skip] {n}: empty"); continue
        hs.Scale(1. / hs.Integral()); hb.Scale(1. / hb.Integral())
        seps[n] = separation(hs, hb)

        c = ROOT.TCanvas("cbr_" + n, "", 900, 620)
        c.SetBottomMargin(0.14)
        for h, col in ((hb, ROOT.kAzure + 2), (hs, ROOT.kRed + 1)):
            h.SetLineColor(col); h.SetLineWidth(2)
            h.SetFillColorAlpha(col, 0.18)
        hb.SetTitle(f";{br_label(n)};normalized")
        # spelled-out names are long: shrink and re-seat the axis title so it fits
        hb.GetXaxis().SetTitleSize(0.036)
        hb.GetXaxis().SetTitleOffset(1.25)
        lo, hi = _zoom([hs, hb])
        if lo is not None and hi > lo:
            hb.GetXaxis().SetRangeUser(lo, hi)
        hb.SetMaximum(1.35 * max(hs.GetMaximum(), hb.GetMaximum()))
        hb.SetMinimum(0.)
        hb.Draw("HIST"); hs.Draw("HIST SAME")

        leg = ROOT.TLegend(0.60, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0); leg.SetFillStyle(0)
        leg.AddEntry(hs, "merged W", "L")
        leg.AddEntry(hb, "QCD (normal)", "L")
        leg.Draw()
        tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.042)
        tx.DrawLatex(0.16, 0.85, f"S = {seps[n]:.3f}")
        c.SaveAs(f"{BR_OUT}/{n}.png")

    order = sorted(seps, key=seps.get, reverse=True)

    def plain(n):   # ROOT latex -> readable markdown
        for a, b in (("#rightarrow", "->"), ("#tau", "tau"), ("#mu", "mu"),
                     ("#beta", "beta"), ("#eta", "eta"), ("#phi", "phi"),
                     ("_{T}", "T"), ("_{h}", "h"), ("_{SD}", "SD"),
                     ("_{2}", "2"), ("_{3}", "3"), ("_{1}", "1"), ("_{4}", "4"),
                     ("{", ""), ("}", ""), ("#", "")):
            n = n.replace(a, b)
        return n

    with open(f"{BR_OUT}/separation_table.md", "w") as fh:
        fh.write("| variable | histogram | S (physics weights) |\n|---|---|---|\n")
        for n in order:
            fh.write(f"| {plain(br_label(n))} | `{n}` | {seps[n]:.3f} |\n")

    h = ROOT.TH1F("brrank", ";;separation S", len(order), 0, len(order))
    for i, n in enumerate(order[::-1]):
        h.SetBinContent(i + 1, seps[n])
        h.GetXaxis().SetBinLabel(i + 1, br_label(n))
    c = ROOT.TCanvas("crank", "", 1250, 1450)
    c.SetLeftMargin(0.50); c.SetGridx()   # wide margin: names are spelled out
    c.SetTopMargin(0.02); c.SetRightMargin(0.03); c.SetBottomMargin(0.06)
    h.SetFillColor(ROOT.kOrange + 1); h.SetBarWidth(0.8); h.SetBarOffset(0.1)
    h.GetXaxis().SetLabelSize(0.016)          # category names
    h.GetYaxis().SetLabelSize(0.016)          # the S scale (drawn horizontally)
    h.GetYaxis().SetTitleSize(0.020); h.GetYaxis().SetTitleOffset(1.30)
    h.SetStats(0)
    h.Draw("hbar")
    c.SaveAs(f"{BR_OUT}/separation_ranking.png")

    print(f"\n  {'branch':28s} {'S':>7s}")
    for n in order:
        print(f"  {n:28s} {seps[n]:7.3f}")
    print(f"\n[done] {len(order)} panels -> {BR_OUT}")

if __name__ == "__main__":
    if os.environ.get('MODE') == 'branches':
        plot_branches()
    else:
        plot_fractions()
        plot_fraction_stack()
        plot_dr_profiles()
        plot_stage3()
        plot_slices("ChMult", "N_{charged} per jet", xmax=50.)
        plot_slices("NConst", "N_{constituents} per jet", xmax=70.)
        plot_ak8()
        print(f"[done] plots -> {OUTDIR}")
