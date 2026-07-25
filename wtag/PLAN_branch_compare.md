# All FatJet branches: merged W vs QCD — plan and results

**Status: done 2026-07-25.** 50 branch panels + ranking + table in
`Tools/Plots/v2b/branches/`. Analyzer histograms in tag `v2b`; the v2 `jets`
tree, `ML/data/dataset.parquet` and the v2 figures were not touched.

Goal: for every raw NanoAOD **FatJet (AK8)** branch in `branches.json`, draw the
merged-W (signal) vs QCD (background) shape comparison with the TMVA separation
`S` printed per panel, plus a ranking table.

## Decisions (locked 2026-07-25)

- **Collection: AK8 / FatJet only.** A merged W only exists as a large-R jet
  (ΔR(qq′) ≈ 2m_W/p_T ≈ 0.8 at 200 GeV > R=0.4), so the AK4 `Jet` block in
  `branches.json` is out of scope.
- **Selection: reconstructed AK8 jet p_T ∈ [200, 1200], m_SD ∈ [20, 250].**
  Cut on *jet* p_T (both classes have it), not gen W p_T (QCD has none). Same
  boundary as the v2 ML dataset, so the numbers stay comparable with v2.
- **merged W = label 1** (FatJet gen-matched to a hadronic W, both daughters
  ΔR<0.8, extra hard b vetoed). **QCD = label 0.** label 2 (W+b / top) is a
  spectator and is not drawn.
- **Physics (xsec) weights only. No flat-p_T reweighting** — the technique was
  already tried on the H→ZZ→4μ BDT and did not pay off, and this is a
  descriptive shape comparison, not a training-input study.
- **Full statistics**, no `--reduction`.

## Approach: histograms in the analyzer

The comparison histograms are filled **inside `Wtag.cc`** and read back by
`Tools/wtag_plots.py`, like the v1 `FatJet/*` histograms already there.

`Wtag.cc` also writes a `jets` TTree (one entry per AK8 jet, 55 branches) — the
ML training dataset that `ML/make_dataset.py` skims into `dataset.parquet`.
**That tree was left exactly as it was.** Adding the new branches there instead
would have meant re-skimming and re-packing a second parquet just to draw
shapes, with `dataset.parquet` — the inputs the trained v2 BDT/ParT models
expect — in the blast radius.

### No `DataFormats/` work was needed

- `Analyzers/include/SKNanoLoader.h:489-541` already declares every FatJet
  branch used here; `AnalyzerCore::GetAllFatJets()` already fills them into the
  `FatJet` object. Only the public *getters* in `DataFormats/include/FatJet.h`
  are missing — and they aren't needed, because `Wtag.cc` reads loader arrays
  directly by fatjet index (the `FatJet_jetId[ifj]` precedent at line 221).
  The index is safe: `GetAllFatJets()` loops `i = 0..nFatJet` with no `continue`,
  so RVec index == NanoAOD index.
- ⇒ Analyzers-only change, `rebuild.sh`, no LinkDef or dictionary regeneration.

## What was done

1. **`Wtag.cc`** — added a `Branch/{W,QCD}_*` block in the fatjet loop after
   `label` is computed: 51 histograms per class. Binning: tagger scores 200 bins
   over [−1,1] (keeps the "undefined" sentinel visible instead of clipped),
   energy fractions 100×[0,1], multiplicities 150×[0,150], masses 120×[0,300],
   τ_i and τ ratios 100×[0,1], `area` 100×[0,4], hadron counts 10×[0,10].
   ParticleNet fills are gated on `Run == 3` (those loader arrays are only
   resized for Run 3). Prefixes: **`PNetM_` = with-mass, `PNet_` =
   mass-decorrelated**, which also resolves the `particleNetWithMass_QCD` vs
   `particleNet_QCD` name collision.
2. **Syntax check + `rebuild.sh`** in the el9 sif — clean.
3. **Local validation** (20k events, TTLJ + QCD_Pt-470to600): m_SD peaks at the
   W mass, `area` = 2.11 ≈ πR² = 2.01, `PNetM_wvsqcd` 0.92 (W) vs 0.19 (QCD),
   b-veto working, no under/overflow.
4. **Condor** (tag `v2b`): 13 samples × 2 eras = 26 DAGs, 999 nodes, all with
   nodes. 4 `hadd` nodes held on the 8 GB cgroup limit (measured 9.8 GB) —
   released at 12–20 GB. **The memory pressure is from hadding the 20M-entry
   `jets` tree, not the histograms.** Next time, raise `--memory` at submit.
5. **Merge** — summed only the `Branch/` histograms across samples rather than
   hadding the trees: 76 GB of input → `v2b/SIG.root` (170 KB) + `v2b/QCD.root`
   (93 KB). Signal class from TTLJ/TTJJ/WW/WZ, QCD class from the QCD bins only:
   `label == 0` also exists in the signal samples, and `make_dataset.py`
   deliberately takes the QCD class from QCD samples only.
6. **Plots** — `plot_branches()` added to `Tools/wtag_plots.py`
   (`MODE=branches`). Axis labels and the table are spelled out in full
   ("charged particle multiplicity", not `chmult`); the raw histogram name is
   kept in the table's second column for traceability.

### Reproduce

```bash
# 1. produce (from the SKNanoAnalyzer repo root, tamsa2)
source setup.sh >/dev/null 2>&1 && bash Wtag.sh          # tag v2b
# 2. merge Branch/ histograms -> v2b/SIG.root + v2b/QCD.root
# 3. plot
cd Tools && MODE=branches python3 wtag_plots.py
```

## Results

Verification: v2b reproduces v2 exactly — TTLJ 2022 gives **20,735,955** jets in
both. All 26 outputs opened clean (no zombie/recovered files), QCD samples carry
zero `W_` histograms as they must.

Top separators (physics-weighted TMVA `S`; full list in
`Tools/Plots/v2b/branches/separation_table.md`):

| variable | histogram | S |
|---|---|---|
| ParticleNet (with mass): W→qq vs QCD | `PNetM_wvsqcd` | 0.674 |
| ParticleNet (with mass): QCD score (sum) | `PNetM_qcd` | 0.648 |
| ParticleNet (with mass): Z→qq vs QCD | `PNetM_zvsqcd` | 0.642 |
| ParticleNet (with mass): top→bqq vs QCD | `PNetM_tvsqcd` | 0.500 |
| ParticleNet (with mass): H→cc vs QCD | `PNetM_hccvsqcd` | 0.495 |
| soft-drop mass | `msoftdrop` | 0.463 |
| ParticleNet regressed mass | `mreg` | 0.442 |
| AK8 jet invariant mass | `mass` | 0.426 |
| N-subjettiness ratio τ₂/τ₁ | `tau21` | 0.396 |
| energy-correlation ratio N₂ (β=1) | `n2b1` | 0.325 |
| ParticleNet (with mass): H→VV→4q vs QCD | `PNetM_h4qvsqcd` | 0.301 |
| N-subjettiness τ₁ | `tau1` | 0.278 |

The ordering is the expected physics: the dedicated W tagger first, then the
other with-mass ParticleNet heads (same network, correlated), then the mass
variables, then 2-prongness. The mass-decorrelated heads sit far lower — by
construction they have the mass information removed, which is most of the
signal here.

### Sanity checks

- **`phi` → S = 0.000** and **`eta` → 0.007**: flat, exactly as they must be.
- **`pt` → 0.017**: after xsec weighting the two classes have nearly the same
  p_T spectrum (mean 276 vs 264 GeV). A large `pt` separation would have meant
  the weighting was wrong.

### Three branches carry no information in Run3 NanoAODv13

Worth knowing before anyone tries to use them:

- **`btagHbb` is a constant −1000** — the tagger is not computed. It is present
  in the NanoAOD (and in `branches.json`) but has a single distinct value.
  Excluded from the plots → 50 panels, not 51.
- **`n3b1`: 99.93% (W) / 100.00% (QCD) sit at the −1 sentinel.** N3 is not
  stored, though N2 is (`n2b1` has a real distribution, S = 0.325). Its
  S = 0.000 means "both classes are the same delta function", **not** "no
  discrimination" — do not read it as a physics result.
- **`lsf3`: 58% / 64% sentinel**, so its S = 0.011 is mostly dilution.

### Excluded by design

- `jetId` — `Wtag.cc:221` already requires tight (`jetId & 2`), so every jet
  reaching the fill passes and the distribution is a delta function.
- `rawFactor` — JEC bookkeeping (`pT_raw = pT·(1−rawFactor)`), a property of the
  calibration rather than the jet.
- `hadronFlavour` — MC truth, and `nBHadrons`/`nCHadrons` carry the same
  composition information. *(`rawFactor` and `hadronFlavour` are also the only
  two FatJet branches missing from `SKNanoLoader`, so including them would mean
  editing that shared file.)*
- Index branches `subJetIdx1/2`, `electronIdx3SJ`, `muonIdx3SJ`,
  `genJetAK8Idx` — no meaningful distribution.
- `nFatJet` — event-level multiplicity, no per-jet signal/background meaning.
- `nbhad`/`nchad` are drawn but are **MC truth**: QCD-composition context only,
  never tagger inputs.

## Not done

- AK4 `Jet` branches; flat-p_T reweighting; any dataset/parquet rebuild or model
  re-training — this task was descriptive branch comparison only.
- The step-5 merge was run inline rather than saved as a `Tools/merge_v2b.sh`
  alongside `merge_v1.sh`/`merge_v2.sh`.
