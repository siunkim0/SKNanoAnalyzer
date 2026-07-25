# Boosted W Tagger v2 — ML tagger on NanoAOD high-level features

**Builds on v1** (PLAN.md): AK8 captures the W above ~220 GeV; ChMult separates merged-W
from quark jets, (SD-)mass from gluon. v2 turns those findings into an actual trained
tagger with quantified performance.

---

## Direction decision (and why not a full Particle Transformer)

**Verified 2026-07-16:** our NanoAOD (Run3NanoAODv13p1) has **no PF-candidate branches**
— a true ParT/ParticleNet needs per-constituent inputs that only PFNano/CMSSW provides.
Producing PFNano is out of scope. Two honest routes remain, we take both:

1. **XGBoost BDT** on the full high-level feature set — robust, interpretable,
   feature-importance quantifies the v1 ChMult hypothesis.
2. **"ParT-lite"**: a small attention model (ParT-style class-attention pooling,
   ~10⁵ params) over the *token set we do have*: 2 subjets + cone-associated SVs +
   a global feature vector. Keeps the transformer architecture on available inputs.

**Benchmarks** (all computable from stored branches — this is the honest yardstick):
- τ₂₁ cut and τ₂₁ + SD-mass window (CMS Run 2 classic W tag)
- N₂^DDT (mass-decorrelated cut-based, CMS style; quantile map in (ρ, p_T))
- `particleNetWithMass_WvsQCD` (constituent-based ceiling, nominal)
- `particleNet_XqqVsQCD` (mass-decorrelated ceiling)

**CMS/ATLAS grounding:** CMS Run 2/3 uses ParticleNet(-MD)/DeepAK8 and N₂DDT; ATLAS's
standard hadronic-W tagger is 3-variable: m_comb + D₂ + **n_trk** — our v1 ChMult
finding is exactly the ATLAS third variable, and NanoAOD v13 stores FatJet
`chMultiplicity` directly. v2 tests how much of the constituent-tagger gap HLFs close,
and whether ChMult adds rejection on top of mass+τ₂₁ (ATLAS-style).

**Mass decorrelation:** train each model twice — *nominal* (mass features in) and *MD*
(mass features out + QCD reweighted flat in (p_T, m_SD)). Sculpting metric: JSD between
pass/fail QCD m_SD shapes vs background efficiency, compared to N₂DDT.

---

## Samples (era: 2022 + 2022EE, v13 DB — extend to 2023/2023BPix later if stats-limited)

| role | sample | status in DB |
|---|---|---|
| signal (hadronic W) | `TTLJ_powheg`, `WW_pythia`, `WZ_pythia` | registered (v1 set) |
| signal (2 had-W/event) | `TTto4Q` | **on /gv0, must register** |
| background | inclusive `QCD_PT-{170to300 … 2400to3200}` (~9 bins) | **on /gv0, must register** (only Mu/EM-enriched are in DB) |

Registration = add to `data/Run3_v13_Run2_v9/<era>/Sample/CommonSampleInfo.json` +
`ForSNU/<alias>_*.json` per era, xsec from XSDB. Known procedure (zero-node DAG gotcha).

## Jet selection & labels (per-AK8-jet)

- AK8: p_T > 200, |η| < 2.4, tight ID.
- **label = W**: gen W (PID ±24, isLastCopy, hadronic) with ΔR(W, jet) < 0.6 and both
  quark daughters within ΔR(q, jet) < 0.8; **veto** jets that also contain the gen b of
  the same top (top-contaminated → separate label `t`, excluded from W-vs-QCD training).
- **label = QCD**: any selected AK8 jet in the inclusive QCD samples.
- TT/WW/WZ non-matched jets kept with label `other` (spectator, not trained on).

## Stage 0 — Framework extension (build.sh; DataFormats+LinkDef+loader — flag: beyond Analyzers/, same pattern as approved v1 FatJet getters)

- `FatJet.h` + `SKNanoLoader`: add `n2b1, n3b1, chMultiplicity, neMultiplicity,
  chHEF, neHEF, chEmEF, neEmEF, muEF, particleNet_massCorr` (+ XqqVsQCD getter if missing).
- New minimal `SubJet` (pt/eta/phi/mass, btagDeepB, tau1-4, n2b1) and `SV`
  (pt/eta/phi/mass, dlenSig, dxySig, ntracks, chi2) classes + loader arrays.

## Stage 1 — Analyzer: extend `Wtag.cc` in place (v1 hists kept) with a per-jet TTree

Branches: label; era/sample id; event weight (genWeight only — shape study);
kinematics (pt, eta, phi, mass, msoftdrop, msoftdrop×massCorr);
substructure (tau1..4 → τ21/τ32 offline, n2b1, n3b1, lsf3);
multiplicity (nConstituents, chMult, neMult, energy fractions);
subjets (pt1, pt2, m1, m2, ΔR(sj,sj), z = pt2/(pt1+pt2), btagDeepB1/2);
SVs in ΔR<0.8 (nSV, Σmass, max dlenSig);
benchmarks (particleNetWithMass_WvsQCD, particleNet_XqqVsQCD);
gen (W p_T, ΔR(qq̄'), daughter flavours, top-contamination flag).
Output via `NewTree`/`FillTrees` (AnalyzerCore built-ins).

## Stage 2 — Production (ask before submitting; tamsa2; tag `v2`)

Register samples → submit Wtag v2 on {TTLJ, TTto4Q, WW, WZ, QCD×9} × {2022, 2022EE}
→ merge to `hi/wtag/2022_v2/` (hadd per process; keep QCD bins xsec-weighted).

## Stage 3 — Training (`hi/wtag/ML/`; **runs on the user's GPU server** — prepare portable code + merged dataset here, user launches the session there)

- Split 60/20/20 by event; class-balance; **flatten p_T** (both classes) so the tagger
  can't learn the spectrum; MD variants additionally flat in (p_T, m_SD) for QCD.
- Models: XGBoost nominal & MD; ParT-lite (PyTorch CPU) nominal & MD.
- Baselines computed on the same test set.

## Stage 4 — Evaluation / deliverables

1. ROC + table: bkg rejection @ ε_sig = 0.5, in p_T bins 200–300 / 300–500 / 500–800 / 800+.
2. ε_sig and mistag vs p_T at fixed mistag WPs (1%, 2.5%, 5%).
3. Mass-sculpting: QCD m_SD pass/fail shapes + JSD vs ε_bkg; nominal vs MD vs N₂DDT.
4. Feature importance → does ChMult add rejection beyond mass+τ₂₁? (v1 hypothesis, ATLAS n_trk analogue.)
5. Results section appended here; plots in `Tools/Plots/v2/`.

**Success criteria:** beat τ₂₁+mass by ≳×2–3 bkg rejection at ε_s=0.5; quantify remaining
gap to ParticleNet (expected to stay ahead — it sees constituents); MD sculpting ≈ N₂DDT.

## Status (2026-07-17)

- **Stage 0 DONE**: FatJet getters (N2b1/N3b1/Ch-NeMultiplicity/energy fractions/
  PNet MD + massCorr), SubJet & SV loader arrays (raw, no DataFormats class);
  full build OK. ⚠ found framework bug: `FatJet::PassTight()` is a no-op
  (`SetJetID` never fills `j_jetId`) — worked around in Wtag.cc via raw
  `FatJet_jetId[i] & 2`.
- **Stage 1 DONE**: Wtag.cc fills `jets` tree (55 branches incl. 3 dlenSig-ordered
  SV tokens). Validated on 20k TTLJ + QCD470 events: label-1 W jets ⟨m_SD⟩=79 GeV,
  τ21=0.35, ⟨PNet⟩=0.86; label-2 (t→bqq̄) ⟨m_SD⟩=146; QCD all label-0.
  N2/N3 carry a sentinel for raw pT<250 — sanitized in ML/make_dataset.py.
- **Samples REGISTERED** (both eras): TTJJ_powheg (nmc 53.2M/179.2M, sumsign
  ratio 0.9919) + 9 inclusive QCD_Pt bins, exact Runs-tree sums
  (`Tools/register_v2.py`, idempotent).
- **ML package DONE** (`ML/`): make_dataset.py (tamsa2) → parquet; train_bdt.py,
  models.py + train_partlite.py, evaluate.py (GPU server). Wtag.sh rewritten for
  v2 (13 samples × 2 eras); Tools/merge_v2.sh ready.
- **Stage 2 DONE (2026-07-17)**: 26/26 DAGs completed on tamsa2 (870M jets on
  /gv0). 4 final-hadd nodes needed memory bumps (8→16 GB; QCD_Pt-600to800
  2022EE needed 32 GB) via condor_qedit+release. Skipped the /data6 tree hadd
  (Tools/merge_v2.sh unused for trees): ML/make_dataset.py reads /gv0 directly,
  skims with per-sample caps (QCD 800k, sig 3M per sample-era; keep-prob folded
  into `weight`), packs to `ML/data/dataset.parquet` (1.4 GB, atomic write).
  Final: 4.36M W (label 1) / 11.9M QCD / 805k t→bqq̄ spectators (y=−1),
  60/20/20 splits; W ⟨m_SD⟩ = 88.6 GeV. Skims in ML/data/skim/ (24M jets).
- **NEXT**: user starts a session on the GPU server → copy `ML/` (with data/),
  `pip install -r requirements.txt`, run train_bdt.py → train_partlite.py →
  evaluate.py (see ML/README.md).

## Results (v2, test split, AK8 p_T 200–1200, 2022+2022EE)

Background rejection (1/ε_b) at ε_s = 0.5, inclusive (plots + per-p_T-bin
numbers in `ML/results/`):

| tagger | rej@ε_s=0.3 | rej@ε_s=0.5 | rej@ε_s=0.7 | JSD@5%ε_b |
|---|---|---|---|---|
| ParticleNet WvsQCD (ceiling) | 299 | **86** | 30 | 0.38 |
| **ParT-lite nominal** | 136 | **48** | 19 | 0.50 |
| **BDT nominal** | 134 | **46** | 18 | 0.50 |
| τ21 + mass window | 68 | 23 | 9 | — |
| τ21 | 40 | 15 | 7 | — |
| **BDT MD** | 31 | 13 | 6 | 0.046 |
| **ParT-lite MD** | 32 | 12 | 5.5 | **0.011** |
| N2^DDT(5%) | 17 | 7 | 3 | 0.001 |
| ParticleNet XqqVsQCD (MD) | 15 | 5.4 | 2.5 | 0.010 |

**Conclusions.**
1. Success criterion met: nominal HLF taggers beat τ21+mass by **×2.1** at
   ε_s=0.5 (48 vs 23); the remaining gap to constituent-based ParticleNet is
   ×1.8 — NanoAOD HLFs recover most of the classic→ML gain without PFCands.
2. ParT-lite ≈ BDT inclusively but **wins in every p_T bin** (500–800: 42 vs 36;
   800–1200: 42 vs 39) — the inclusive tie is a binning artifact.
3. **Mass decorrelation works**: ParT-lite MD reaches JSD 0.011 at 5% ε_b
   (≈45× less sculpting than nominal), matching the stored MD ParticleNet
   XqqVsQCD sculpting while giving **2.2× its rejection** (12 vs 5.4) and 1.7×
   N2DDT's (12 vs 7). bdt_md sculpts ~4× worse than part_md at same rejection.
4. v1 ChMult hypothesis quantified: **chmult is the #2 BDT feature in both
   variants (16% of total gain)** — behind sdmass (nominal) / τ32 (MD) — the
   ATLAS-n_trk pattern reproduced in CMS NanoAOD.
5. Eval-time fix (2026-07-23): dataset had 8 NaN sv_summass rows → NaN norm in
   ParT checkpoints (training had nan_to_num, eval didn't) → all-NaN scores.
   evaluate.py now mirrors the guard; make_dataset drops NaN rows; norm uses
   nanmean/nanstd. Trained models effectively ignore sv_summass (zeroed column).

**Deliverables (2026-07-23):** working points at fixed mistag (5/2.5/1/0.5%) in
`ML/results/wp_table.{md,tex,json}` (`ML/make_wp_table.py`); full write-up as a
6-page paper: `paper/wtag_v2.tex` → `paper/wtag_v2.pdf` (pdflatex ×2 on
ai-tamsa1; figures pulled from `ML/results/`, WP table `\input` from
`ML/results/wp_table.tex` — regenerate table, recompile, done).

## Decisions (approved 2026-07-16)

1. Direction as planned (BDT + ParT-lite, MD study, full benchmarks). Era 2022+2022EE.
2. Register TTto4Q + inclusive QCD bins in DB — approved.
3. Extend `Wtag.cc` in place with TTree — approved.
4. DataFormats/loader extension — approved.
5. **Training on the user's GPU server** (not this machine); prepare code/dataset here,
   user will start a session there. Still ask before condor submission.
