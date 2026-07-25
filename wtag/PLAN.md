# Boosted W Kinematics → Boosted W Tagger — Plan

**Goal.** Study the hadronic decay W→qq' to build a boosted-W tagger. Concretely:
1. Find at what W p_T the two decay quarks merge into a **single detector jet**.
2. Find at what W p_T the two jets start to **overlap**.
3. Check jet track/substructure kinematics — hypothesis: a **merged (boosted) W jet
   carries more tracks than an ordinary single (light-quark/gluon) jet**.

**Decisions (locked in).**
- Samples: **tt̄ semileptonic** (resolved→merge transition) **+ a high-p_T source**
  (WW / heavy resonance, to populate the fully-merged regime).
- Jet cones: **both AK4 (R=0.4) and AK8 (R=0.8)**.
- Era/campaign: **Run 3 2022**.
- We use MC, so we cut on **generator W p_T** directly (see "Using MC truth" below).

---

## Physics anchor

Two-body opening angle scales as **ΔR ≈ 2·m_W / p_T** (m_W ≈ 80.4 GeV):

| W p_T  | ΔR(q,q') | regime                              |
|--------|----------|-------------------------------------|
| ~160   | ~1.0     | well resolved, two AK4 jets         |
| ~200   | ~0.8     | starts fitting in one **AK8** (R=0.8) |
| ~400   | ~0.4     | starts merging into one **AK4** (R=0.4) |

So "how much p_T makes it one jet" has **two answers** (AK4 ≈ 400, AK8 ≈ 200 GeV);
Stage 2 measures both from full simulation instead of the back-of-envelope formula.

---

## Using MC truth (the "W p_T > n GeV" question)

Because it is MC, the W is available in the `Gen` collection (PID ±24, `isLastCopy`),
so we read its **true p_T** and can select `genW.Pt() > n` directly. Valid and standard.

- A cut **selects, it does not create** events → inclusive tt̄ is starved at high p_T,
  which is why we add a high-p_T sample to *populate* the boosted regime.
- **Gen W p_T** is the clean x-axis for Stages 1–2 (merge/overlap curves).
  **Reco jet p_T** is the matching variable for Stage 3 (what a real tagger sees).
- No generator filter needed at production: scan the gen-W-p_T threshold *in analysis*,
  which gives the "> n GeV" view at every threshold at once.

---

## Framework feasibility (confirmed in SKNano DataFormats)

- **AK4 `Jet`**: `Pt()`, `M()` (mass), `Eta()`, `DeltaR()`;
  `nConstituents()` (total), `chMultiplicity()` (**charged ≈ # tracks**),
  `neMultiplicity()` (neutral), `nSVs()` (# secondary vertices),
  QG discriminators (`btagPNetQvG`, `btagDeepFlavQG`).
- **AK8 `FatJet`**: `SDMass()` (softdrop), `tau1..4` (→ τ₂₁), subjet indices,
  `NConstituents()`, ParticleNet `WvsQCD`.
- **`Gen`**: `PID()`, `MotherIndex()`, status flags (`isLastCopy`, `fromHardProcess`).
- `Analyzers/src/Wtag.cc` exists but is **empty** — this is the skeleton to fill
  (in main repo `folder/sknano/new/SKNanoAnalyzer`). New `Wtag.h` needed too.

> **Observable suite mirrors `hi/track`** (`Analyzers/src/track.cc`): total /
> **charged** / **neutral** multiplicity, per-flavour, in jet-p_T bins.
> The charged & neutral branches fill correctly **only on Run3 2022+ (NanoAODv12+)** —
> we are on that era. `chMultiplicity` is the standard proxy for "# tracks in the jet".

> **Impact parameters (dxy/dz/ip3d) — decision: nSVs only.** This framework has no
> Track / PFCand / SecondaryVertex collection; IP is a *lepton* property, not jet-native
> (in `track.cc` dXY/dZ/IP3D/SIP3D are the muon's). Per-jet/per-track IP is not reachable
> without PFNano. → We use jet-native **`nSVs()`** as the displacement/lifetime proxy for
> every jet; no soft-lepton IP.

---

## Stage 0 — Sample & analyzer setup

- First action on "go": grep the Run3_2022 sample DB
  (`data/<ver>/2022/Sample/CommonSampleInfo.json`) for registered tt̄ + high-p_T
  hadronic-W samples; bring the concrete list back before submitting.
- Fill `Wtag.cc` (+ `Wtag.h`). **DONE.** Samples (v13 DB): `TTLJ_powheg`,
  `WW_pythia`, `WZ_pythia` × {2022, 2022EE}, via `Wtag.sh`, tag `v1`.
- **Output = `FillHist` histograms (not a flat ntuple).** Decision: the whole `hi/track`
  Tools/plotting pipeline is histogram-based and every deliverable (fraction-vs-p_T
  curves, p_T-binned multiplicity) is a 1D ratio or 2D `*_vs_Pt` profile — matches the
  proven `track.cc` pattern and plugs straight into existing Tools. Histograms written:
  - `Gen/`: `W_Pt_all` (denom), `dRqq_vs_WPt`, `W_Pt_dRqq_lt04/_lt08`
  - `Reco/`: `W_Pt_ak4bothmatched/_ak4merged/_ak4resolved/_ak4overlap/_ak8captured`,
    `dRjj_vs_WPt`
  - `Jet/`: `Pt,Mass,NConst,ChMult,NeMult,NSV` (+ `*_vs_Pt`) with category suffix
    `_mergedW / _singleW / _lightuds / _gluon / _c / _b`
  - `FatJet/`: `Pt,SDMass,NConst` (+ `*_vs_Pt`) for `_capturedW` and `_other`
- **AK8 tagger vars (τ₂₁, ParticleNet WvsQCD) DONE** — added `Tau1()..Tau4()` +
  `ParticleNetWithMass_WvsQCD()` getters to `FatJet.h` (build.sh, rootcling regen);
  `FatJet/Tau21_*` and `FatJet/WvsQCD_*` (capturedW + other) now filled.
- Output → `/gv0/Users/snuintern2/SKNanoOutput/Wtag/<era>_v1/`; plots in `hi/wtag`.
- **SUBMITTED (2026-07-16)** all 6 sample×era DAGs (41 nodes each) on tamsa2, running.
  ⚠️ tamsa2 condor schedd IDTOKENS auth is intermittently flaky — 4/6 failed on first
  pass with "Failed to connect to schedd / authenticate using IDTOKENS"; resubmitted the
  failed ones with a retry loop (verify with `condor_q` which runlogs actually queued).
- **COMPLETED (2026-07-16 ~14:00)**: all 6 outputs on /gv0, queue empty; merged via
  `Tools/merge_v1.sh` → `hi/wtag/2022_v1/{TTLJ,WW,WZ}.root` (2022+2022EE hadded).
  `Tools/wtag_plots.py` extended to the full deliverable set (ΔR profiles + 2m_W/p_T
  overlay, all fraction curves + 50%-crossing printout, NeMult/NSV profiles, p_T-sliced
  ChMult/NConst shapes, AK8 τ₂₁/WvsQCD/NConst); plots in `Tools/Plots/v1/`
  (`_all`/`_TTLJ`/`_WW` variants).

## Results (v1, TTLJ+WW+WZ combined; 50% crossings in gen W p_T)

| curve | 50% at | naive 2m_W/p_T |
|---|---|---|
| ΔR(qq') < 0.8 (gen)   | ~227 GeV | 200 |
| AK8 captured (1 fatjet)| ~220 GeV | 200 |
| AK4 jets overlap ΔR<0.8| ~222 GeV | — |
| ΔR(qq') < 0.4 (gen)   | ~462 GeV | 400 |
| AK4 merged (1 jet)    | ~444 GeV | 400 |

Per-sample stable: TTLJ 436/218 GeV, WW 476/241 GeV (AK4-merged / AK8-captured).
Merged fraction plateaus ~90% (AK4, p_T>900), AK8 capture ~95%.

**Hypothesis (Stage 3) confirmed, p_T-matched.** ⟨N_charged⟩ at jet p_T 200–400 GeV:
merged-W **20.7** vs gluon 17.9 vs single-q-W 15.2 vs light-uds 13.6. At 400–600 GeV
merged-W 21.7 ≈ single-q 21.2, gluon 24.0 — the multiplicity edge over quark jets
persists but gluon overtakes at high p_T; ChMult alone separates merged-W from *quark*
jets, mass/SD-mass separates it from gluon (AK8 SD-mass peaks ~85 GeV for captured W
vs falling QCD shape; small t→bqq̄ peak at ~175 GeV from TTLJ full-top capture).

## Stage 1 — Gen truth: boost vs opening angle

- W (PID ±24, `isLastCopy`) → two quark daughters (|PID|≤6, `fromHardProcess`).
- Record genW p_T, ΔR(q₁,q₂), quark p_T.
- **Plots:** 2D ΔR vs genW p_T (+ profile); fraction ΔR<0.4 and ΔR<0.8 vs genW p_T.

## Stage 2 — Detector level: merging & overlap

- Match each gen quark to reco AK4 jets (ΔR<0.4). Classify each W:
  **resolved** (2 distinct jets) / **merged** (both quarks → same jet) / other.
- **Fraction merged-into-1-AK4 vs genW p_T** — "one jet in detector."
- **Overlap:** ΔR between the two matched AK4 jets; fraction ΔR<0.8 vs p_T.
- **AK8:** fraction of W captured in one FatJet (both quarks within ΔR<0.8 of axis) vs p_T.

## Stage 3 — Jet track/substructure kinematics (the hypothesis)

- **Signal jets:** merged-W AK4 jets. **Controls:** single-quark jets from resolved W +
  generic light-q/gluon jets.
- Compare the full `hi/track`-style suite, per flavour, **in matched reco jet-p_T bins**:
  - `nConstituents` (total), `chMultiplicity` (tracks), `neMultiplicity` (neutral)
  - jet p_T, jet mass, ΔR(jet,jet)
  - `nSVs` (secondary-vertex count — displacement proxy)
  - QG disc (`btagPNetQvG`)
- ⚠️ Multiplicity rises with jet p_T (and |η|) → the p_T-binned comparison is essential,
  else the "2-prong" effect is confounded by the merged jet just being higher-p_T.
- **AK8 cross-check:** `NConstituents`, SDMass (peak ~80 GeV), τ₂₁, ParticleNet WvsQCD.

---

## Deliverables

1. ΔR(qq') vs W p_T curve (gen).
2. Fraction-merged (AK4) and fraction-captured (AK8) vs W p_T; AK4 jet-jet overlap vs p_T.
3. Tracks-per-jet comparison (merged-W vs single jet, p_T-matched) — hypothesis test.

## Open items before build

- Analyzer location: fill `Wtag.cc` in main repo, plots in `hi/wtag` — confirm.
- Concrete Run3-2022 sample names (Stage 0 grep).
