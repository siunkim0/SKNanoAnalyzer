"""Working-point scan: S/B/Asimov-Z vs BDT cut, with cut-based reference.

Reads the per-sample ROOT outputs produced by HiggsBDT and Higgs analyzers
on 2018 MC + data, integrates the m4l histograms in the [120, 130] GeV
window (same window for BDT and cut-based, so the comparison is fair), and
prints + plots Asimov Z vs the BDT working point.

Edit the BASE path and SIGNAL/BKG sample-to-file maps below before running.

Output:
  - stdout: cut-based S/B/Z + per-WP table + best-WP summary
  - plots/wp_scan_2018.png : Z(WP) curve with horizontal cut-based line

Definitions:
  S = sum over {ggH, VBFH, VH}
  B = sum over {ZZTo4L, DY, TTLL, TTLJ}
  Z = sqrt(2 * ((S+B)*ln(1+S/B) - S))      -- Asimov, no syst
  Best WP: argmax over WPs with B >= 0.5  -- avoids low-stat Z blow-up
"""
from __future__ import annotations

import ctypes
from pathlib import Path

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

# --------------------------------------------------------------------------
# CONFIG — edit these paths after the Task 2c production runs finish.
# --------------------------------------------------------------------------
ERA  = "2018_v5_3"
BASE = Path("/data6/Users/snuintern2/Higgs") / ERA

# sample tag -> ROOT filename (under BASE).
SIGNAL = {
    "ggH":  "ggHToZZTo4L.root",
    "VBFH": "VBFHToZZTo4L.root",
    "VH":   "VHToNonbb.root",
}
BKG = {
    "ZZ":   "ZZTo4L_powheg.root",
    "DY":   "DY.root",
    "TTLL": "TTLL_powheg.root",
    "TTLJ": "TTLJ_powheg.root",
}

SYST_DIR = "Central"                                # nominal histogram dir
WPS      = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95]
WINDOW   = (120.0, 130.0)
B_MIN    = 0.5                                      # guard on best-WP pick
OUT_PNG  = Path(f"plots/wp_scan_{ERA}.png")

# Cut-based reference histogram. The Higgs.cc analyzer wasn't run, so we read
# H_Mass from the SAME HiggsBDT output files — that's the BDT analyzer's full
# selection with no BDT cut applied. In the [120,130] counting window this is
# equivalent to Higgs.cc's H_Mass_final (which only adds an m4l > 70 cut),
# and using the same files keeps the comparison apples-to-apples.
CUT_HIST = "H_Mass"
CUT_FILE_OVERRIDE: dict[str, str] = {}  # if cut-based files ever live elsewhere


# --------------------------------------------------------------------------
def yield_window(h: ROOT.TH1, lo: float, hi: float) -> tuple[float, float]:
    ax = h.GetXaxis()
    b1 = ax.FindBin(lo + 1e-6)
    b2 = ax.FindBin(hi - 1e-6)
    err = ctypes.c_double(0.0)
    val = h.IntegralAndError(b1, b2, err)
    return val, err.value


def sum_group(files: dict[str, str], histname: str,
              lo: float, hi: float) -> tuple[float, float]:
    """Sum yield_window over a sample group; quadrature on the error."""
    tot, var = 0.0, 0.0
    for tag, fname in files.items():
        path = CUT_FILE_OVERRIDE.get(tag, BASE / fname)
        f = ROOT.TFile.Open(str(path))
        if not f or f.IsZombie():
            print(f"  WARN cannot open {path}")
            if f: f.Close()
            continue
        h = f.Get(f"{SYST_DIR}/{histname}")
        if not h:
            print(f"  WARN missing {SYST_DIR}/{histname} in {path}")
            f.Close()
            continue
        v, e = yield_window(h, lo, hi)
        tot += v
        var += e * e
        f.Close()
    return tot, float(np.sqrt(var))


def asimov_Z(S: float, B: float) -> float:
    if B <= 0 or not np.isfinite(S) or not np.isfinite(B):
        return float("nan")
    return float(np.sqrt(2.0 * ((S + B) * np.log(1.0 + S / B) - S)))


def main() -> None:
    lo, hi = WINDOW
    print(f"era={ERA}  window m4l in [{lo:.1f}, {hi:.1f}] GeV\n")

    # cut-based reference (one histogram, same window)
    Sc, Sce = sum_group(SIGNAL, CUT_HIST, lo, hi)
    Bc, Bce = sum_group(BKG,    CUT_HIST, lo, hi)
    Zc = asimov_Z(Sc, Bc)
    print(f"cut-based ({CUT_HIST}):")
    print(f"  S = {Sc:8.3f} +- {Sce:6.3f}")
    print(f"  B = {Bc:8.3f} +- {Bce:6.3f}")
    print(f"  Z = {Zc:8.3f}\n")

    print(f"{'WP':>5} {'S':>9} {'B':>9} {'S_err':>8} {'B_err':>8} {'Z_asimov':>10}")
    print("-" * 56)
    rows: list[tuple[int, float, float, float, float, float]] = []
    best: tuple[int | None, float] = (None, -1.0)
    for wp in WPS:
        hn = f"H_Mass_BDT{wp:03d}"
        S, Se = sum_group(SIGNAL, hn, lo, hi)
        B, Be = sum_group(BKG,    hn, lo, hi)
        Z = asimov_Z(S, B)
        rows.append((wp, S, B, Se, Be, Z))
        print(f"{wp/100:5.2f} {S:9.3f} {B:9.3f} {Se:8.3f} {Be:8.3f} {Z:10.3f}")
        if np.isfinite(Z) and B >= B_MIN and Z > best[1]:
            best = (wp, Z)

    if best[0] is not None:
        print(f"\nbest WP = {best[0]/100:.2f}  (Z = {best[1]:.3f}, B >= {B_MIN})")
        if Zc > 0 and np.isfinite(Zc):
            ratio = best[1] / Zc
            delta_pct = (ratio - 1) * 100
            verdict = "improves over" if best[1] > Zc else "does NOT beat"
            print(f"  vs cut-based Z = {Zc:.3f}  ->  BDT {verdict} cut-based "
                  f"by {delta_pct:+.1f}%")
    else:
        print(f"\nNo WP with B >= {B_MIN}; cannot pick a best WP. Inspect table above.")

    # Plot Z vs WP, with cut-based horizontal reference.
    xs = np.array([r[0] / 100 for r in rows], dtype=np.float64)
    zs = np.array([r[5] if np.isfinite(r[5]) else 0.0 for r in rows],
                  dtype=np.float64)
    g = ROOT.TGraph(len(xs), xs, zs)
    g.SetTitle(";BDT score cut;Asimov Z (m4l #in [120,130] GeV)")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)
    g.SetLineWidth(2)

    c = ROOT.TCanvas("c_wp", "wp scan", 800, 600)
    g.Draw("APL")
    if np.isfinite(Zc) and Zc > 0:
        line = ROOT.TLine(float(xs.min()), Zc, float(xs.max()), Zc)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        tex = ROOT.TLatex()
        tex.SetNDC(False)
        tex.SetTextColor(ROOT.kRed)
        tex.SetTextSize(0.035)
        tex.DrawLatex(float(xs.max()) - 0.10, Zc * 1.02,
                      f"cut-based Z = {Zc:.2f}")
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    c.SaveAs(str(OUT_PNG))
    print(f"\nwrote {OUT_PNG}")


if __name__ == "__main__":
    main()
