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
   kept in the table's second column for traceability. Two rankings are written:
   `separation_ranking.png` (all 50) and `separation_ranking_notagger.png`
   (26 — every network output dropped, see below).

### Reproduce

```bash
# 1. produce (from the SKNanoAnalyzer repo root, tamsa2)
source setup.sh >/dev/null 2>&1 && bash Wtag.sh          # tag v2c
# 2. merge Branch/ + BranchNoMSD/ histograms -> v2c/SIG.root + v2c/QCD.root
#    cd wtag/Tools && python3 merge_v2c.py
#    Reads only the histogram directories, never the `jets` tree, so 76 GB of
#    condor output becomes ~0.5 MB in seconds. It also auto-detects samples whose
#    hadd node died (see below) and sums their per-job files instead.
# 3. plot. `rot` (the LCG *centos7* view) does not load on tamsa2 —
#    `libssl.so.10: cannot open shared object file`. Two working alternatives,
#    both verified 2026-07-25:
#      a) the el9 view — `source /cvmfs/sft.cern.ch/lcg/views/LCG_105/\
#         x86_64-el9-gcc12-opt/setup.sh` (ROOT 6.30/02). On a cold CVMFS cache
#         the first `import ROOT` can take >400 s; that is not a failure, so do
#         not wrap it in a short timeout.
#      b) the el9 container (ROOT 6.34.04), which is what was used here:
singularity exec /data6/Users/snuintern2/private-el9.sif zsh -c '
  export PATH=/opt/conda/bin:$PATH; export MAMBA_ROOT_PREFIX=/opt/conda
  cd /data6/Users/snuintern2/wtag/SKNanoAnalyzer; source setup.sh >/dev/null 2>&1
  unset PYTHONHOME PYTHONPATH GIT_EXEC_PATH
  cd wtag/Tools && MODE=branches python3 wtag_plots.py'
```

## v2c supersedes v2b (2026-07-26)

v2b's numbers were produced with three defects; v2c is the rerun that fixes them
and is the version to quote. Plots in `Tools/Plots/v2c/branches/`.

1. `n3b1` on a score-shaped axis → false S = 0.000 (true value 0.027).
2. The `PNet_x*vsqcd` "not evaluated" sentinel is **−10**, so it fell into
   *underflow* and was silently dropped — and at a class-dependent rate
   (1.5% of W vs 3.3% of QCD), so the two classes were normalised over
   different subsets.
3. Every separation was conditional on the m_SD window with no way to see how
   much of it that window was responsible for.

**Validation — v2c's `Branch/` set reproduces v2b exactly where it should.**
Of the 50 common branches, **only 8 moved**, and they are precisely the 8 the two
binning fixes touch:

| branch | v2b | v2c | |
|---|---|---|---|
| `n3b1` | 0.000 | **0.027** | the whole distribution was in overflow |
| `PNet_xtmvsqcd` | 0.002 | 0.006 | |
| `PNet_xttvsqcd` | 0.003 | 0.006 | |
| `PNet_xtevsqcd` | 0.002 | 0.005 | |
| `PNet_xbbvsqcd` | 0.063 | 0.065 | |
| `PNet_xggvsqcd` | 0.064 | 0.066 | |
| `PNet_xqqvsqcd` | 0.147 | 0.146 | |
| `PNet_xccvsqcd` | 0.237 | 0.235 | |

The other 42 are identical to three decimals. So the sentinel bias was **real but
small** (≤0.004) — worth fixing for correctness, but it never changed a
conclusion. The `n3b1` bug was the consequential one. `btag_hbb` now also appears
(51 branches, not 50): its constant −1000 lands in the −1.05 bucket instead of
vanishing into underflow, so it plots honestly as a single spike at S = 0.000
rather than being dropped as "empty".

Panels whose sentinel spike holds >25% of the events are drawn **log-y** — on a
linear axis the spike flattens the real distribution onto the baseline, which is
what made `n3b1` look empty to the eye even after the axis was fixed.

### The hadd nodes fail, and the failure is silent

4 of 26 hadd nodes died on the 8 GB cgroup limit (they need ~9.8 GB), and
`--memory` does **not** help: `python/SKNano.py:447` hardcodes
`request_memory = 8192` for hadd, independent of the flag. The cost is entirely
the 20M-entry `jets` tree, which the branch study never reads.

A dead hadd leaves a **truncated but readable** output file — the failure mode
that becomes a wrong number rather than an error. `merge_v2c.py` detects it
exactly rather than heuristically: `hadd.sh` deletes `output/hists_*.root` only
after hadd returns, so per-job files still present mean that sample's hadd did
not finish. All 4 were recovered at 40/40 job files, so v2c is complete.

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

### Ranking without the taggers

`separation_ranking_notagger.png` / `separation_table_notagger.md` repeat the
ranking over the **26 branches that are not a network output**, i.e. what the
jet observables achieve on their own. Dropped: all `PNetM_*` / `PNet_*`, the
DeepDoubleX (`ddbvl`, `ddcvb`, `ddcvl`) and DeepCSV (`btag_deepb`) discriminants,
and also **`mreg` and `masscorr`** — those two are ParticleNet regression
outputs, a network's *prediction* of the mass rather than the mass built from
the constituents (`mass`, `msoftdrop`).

| variable | histogram | S |
|---|---|---|
| N-subjettiness ratio τ₄/τ₁ | `tau41` | **0.496** |
| N-subjettiness ratio τ₃/τ₁ | `tau31` | **0.484** |
| soft-drop mass | `msoftdrop` | 0.463 |
| AK8 jet invariant mass | `mass` | 0.426 |
| N-subjettiness ratio τ₂/τ₁ | `tau21` | 0.396 |
| energy-correlation ratio N₂ (β=1) | `n2b1` | 0.325 |
| N-subjettiness τ₁ | `tau1` | 0.278 |
| number of c hadrons *(MC truth)* | `nchad` | 0.195 |
| jet catchment area | `area` | 0.142 |
| charged particle (track) multiplicity | `chmult` | 0.055 |

Mass and prongness are the whole story; everything below `area` is under 0.06.
**`nchad` at 8th place is MC truth**, not a measurable input — it is in the plot
for QCD-composition context only, and the labels now say so.

### τ₃/τ₁ and τ₄/τ₁ beat τ₂/τ₁ — the standard choice is not the best one

The three "skip-a-step" ratios were added on 2026-07-26 (`tau31`, `tau41`,
`tau42`) and **τ₄/τ₁ is now the top non-tagger discriminant**, above
`msoftdrop`:

| | τ₂/τ₁ | τ₃/τ₁ | τ₄/τ₁ | τ₃/τ₂ | τ₄/τ₂ | τ₄/τ₃ |
|---|---|---|---|---|---|---|
| with m_SD | 0.396 | 0.484 | **0.496** | 0.032 | 0.032 | 0.004 |
| no m_SD | 0.598 | 0.661 | **0.662** | 0.056 | 0.021 | 0.009 |

The pattern is sharp and it is entirely about the **denominator**: every ratio
with τ₁ under it separates well, every ratio without it is worthless (≤0.06).
τ_N/τ₁ asks "how much better than *one* axis does this jet fit N axes" — and a
2-prong W improves enormously over 1 axis while a QCD jet, already 1-prong-like,
barely improves at all. Going to 3 or 4 axes keeps accumulating that same
1-prong-vs-not contrast (τ_N is monotonically decreasing), so τ₃/τ₁ and τ₄/τ₁
carry *more* of it than τ₂/τ₁, not less. By contrast τ₃/τ₂ and τ₄/τ₃ compare two
already-good fits, which is a question about 3- and 4-prongness — the right
question for a **top** tagger, and nearly no information for W-vs-QCD.

So `tau21` being the conventional W-tagging variable is convention, not
optimality: on this sample τ₄/τ₁ is worth +0.100 in S over it inside the mass
window. Worth feeding τ₃/τ₁ and τ₄/τ₁ to the v2 ML tagger, which currently gets
τ₁–τ₄ raw and has to build such a ratio itself.

**These three came from the `jets` tree, not from a rerun.** τ₁–τ₄ are each
stored per jet, so the ratios were rebuilt offline in `merge_v2c.py`
(`tau_ratios()`) with the same axis and the same two selections. The tree stores
floats truncated to ~11 mantissa bits (p_T lands on 0.125 GeV steps), which moves
~0.03% of jets across a window edge, so the route is not free by construction —
`tau21`/`tau32`/`tau43` are rebuilt the same way alongside and compared to the
analyzer's own histograms as a fidelity check. **All six agree to 0.0000 in S**,
which is what licenses putting the offline numbers in the table above. `Wtag.cc`
now also fills all six natively, so a future rerun needs no offline step.

### The multiplicity branches

There is no separate "track multiplicity" branch: a track *is* a charged
particle, so `chMultiplicity` **is** the track count, and `nConstituents` is the
total over charged + neutral. Confirmed on the merged histograms — merged W:
26.4 + 16.1 = 42.5 ≈ 43.0; QCD: 30.4 + 16.2 = 46.5 ≈ 47.0 (the small excess is
the PUPPI weighting, which `nConstituents` does not apply). The labels are now
"charged particle (track) multiplicity" and "total particle multiplicity
(charged + neutral constituents)".

Note the separation sits almost entirely in the **charged** half (0.055 vs 0.005
neutral): QCD jets have ~4 more tracks than merged W jets, while the neutral
count is the same in both.

### The m_SD window is a strong conditioning — the ranking is conditional

The selection cuts on `sdmass`, which is itself the top non-tagger discriminant,
so **every S in the `with_msd/` set is conditional on m_SD ∈ [20, 250]**. That is
not a cosmetic caveat.

**The window keeps 98.1% of merged W but only 40.9% of QCD** — it discards 59%
of the background, mostly low-mass, before anything is measured.

As of **v2c the analyzer fills both** (`Branch/` and `BranchNoMSD/`), so this is
now measured directly for all 51 branches — including the ParticleNet and
DeepDoubleX scores, which are *not* in the `jets` tree and so could not be
re-measured offline. Full table in `Tools/Plots/v2c/branches/msd_comparison.md`;
largest shifts:

| variable | with m_SD | no m_SD | Δ |
|---|---|---|---|
| `tau1`          | 0.278 | 0.533 | **+0.254** |
| `mass`          | 0.426 | 0.641 | +0.215 |
| `tau21`         | 0.396 | 0.598 | +0.202 |
| `msoftdrop`     | 0.463 | 0.662 | +0.199 |
| `PNetM_hccvsqcd`| 0.495 | 0.664 | +0.170 |
| `PNetM_h4qvsqcd`| 0.301 | 0.465 | +0.164 |
| `mreg`          | 0.442 | 0.603 | +0.161 |
| `masscorr`      | 0.151 | 0.307 | +0.156 |
| `PNetM_wvsqcd`  | 0.674 | 0.782 | +0.109 |
| `nemult`        | 0.005 | 0.088 | +0.082 |
| `PNet_xqqvsqcd` | 0.146 | 0.082 | −0.064 |
| `n2b1`          | 0.325 | 0.262 | **−0.063** |
| `PNet_xccvsqcd` | 0.235 | 0.105 | **−0.130** |

Shifts of up to ±0.25 are comparable to the separations themselves and they
reorder the ranking: without the window `tau1` overtakes `n2b1`, and `nconst`
overtakes `chmult`.

The **negative** deltas are the interesting ones — those variables are genuinely
*more* discriminating inside the mass window than outside. `n2b1`, `chmult` and,
most sharply, the two mass-decorrelated ParticleNet heads `PNet_xccvsqcd`
(−0.130) and `PNet_xqqvsqcd` (−0.064). That is exactly what mass-decorrelation
is built to do: those heads deliberately discard mass information, so they only
pay off once the mass information has been spent by the window. Outside it they
are competing against variables that still have mass to use, and they lose.

Neither column is "the" answer, and the choice is not cosmetic:

- **no m_SD** is the honest *survey* — "what separates W from QCD". A ranking
  measured inside the window is circular for the mass branches, since
  `msoftdrop` is scored on a sample `msoftdrop` has already purified.
- **with m_SD** is the *deployment* number — "what adds information given a jet
  already in the tagger's working range", which is the operationally relevant
  question for v2, because a real analysis always applies a mass window.

The no-window column is also inflated in a partly trivial way (the discarded QCD
is low-mass and separable on mass alone). Quote whichever answers the question
being asked, and say which one it is.

### Sanity checks

- **`phi` → S = 0.000** and **`eta` → 0.007**: flat, exactly as they must be.
- **`pt` → 0.017**: after xsec weighting the two classes have nearly the same
  p_T spectrum (mean 276 vs 264 GeV). A large `pt` separation would have meant
  the weighting was wrong.

### Branches with sentinel / range problems

- **`btagHbb` is a constant −1000** — the tagger is not computed. It is present
  in the NanoAOD (and in `branches.json`) but has a single distinct value.
  Excluded from the plots → 50 panels, not 51.
- **`lsf3`: 57% / 63% sentinel** (at −1), so its S = 0.011 is mostly dilution.
  It is also **not bounded by 1** despite being a "fraction" — measured max
  **33.5**, ~1.5% of entries above 1, all lost to overflow on the score axis.
- **The `PNet_x*vsqcd` sentinel is −10, not −1 — CORRECTED 2026-07-25.** All the
  ParticleNet vs-QCD *ratios* (`particleNet_X{bb,cc,qq,gg,te,tm,tt}VsQCD` and
  `particleNetWithMass_*vsQCD`) write **−10** when the tagger is not evaluated,
  so on the 200×[−1,1] axis they fell into **underflow** and were silently
  dropped by the normalisation. The dropped fraction is **class-dependent —
  1.5% of W vs 3.3% of QCD** — so the two classes were being normalised over
  different subsets, biasing all seven `PNet_x*` separations. The raw
  `particleNet_QCD*` probabilities have no sentinel (min 0), and
  `PNetM_wvsqcd` is barely affected (0.17%), so the headline S = 0.674 stands.
  `Wtag.cc`'s `S()` lambda is now **220×[−1.1, 1.1]** with out-of-range values
  bucketed into the ±1.05 pad: physical binning unchanged at 0.01/bin, the
  undefined spike visible and separable, and nothing dropped. This also fixes
  `lsf3`'s overflow. **Needs the rerun to take effect.**
- **`n3b1` — CORRECTED 2026-07-25. The v2b value of S = 0.000 was a binning
  bug, not a property of the data.** An earlier version of this document
  claimed "N3 is not stored"; that was **wrong**. N3 *is* stored, with a real
  distribution over **[0.039, 4.211], mean 1.569**. Unlike N2 it is not bounded
  by 1, so on the 200×[−1,1] score axis **~94% of the defined entries fell into
  overflow** — and both `TH1::Integral()` (no args) and `separation()` (bins
  1..N) ignore overflow, leaving only the −1 sentinel spike visible in each
  class. Recomputed from the v2b `jets` trees on a [0,5] axis:

  | | S | W mean | QCD mean |
  |---|---|---|---|
  | `n3b1`, defined entries only | **0.024** | 1.654 | 1.704 |
  | `n3b1`, incl. −1 sentinel | 0.027 | | |
  | `n2b1` control | 0.325 (reproduces the histogram value exactly) | 0.274 | 0.351 |

  So N3 really is a weak discriminant — but at 0.024, not 0.000, which puts it
  mid-pack (between `tau2` 0.026 and `pt` 0.017) rather than last.
  `Wtag.cc` now bins it 250×[−1.25, 5]; **the v2b outputs still carry the old
  binning, so the shipped `n3b1.png` panel is not usable** — rerun to fix it.
  `plot_branches()` now warns for any branch with >0.1% outside its axis.

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
- ~~The step-5 merge was run inline rather than saved as a script.~~ Done in
  v2c: `Tools/merge_v2c.py`.
