# Is the stored ParticleNet already good enough? — go/no-go for the NanoAOD tagger

**Status 2026-07-26: CLOSED — the tagger is not being built.** All three steps
ran (tag `v2d`, 870M jets). Step 3 put the surviving justification in the
"marginal" band; the discovery that NanoAODv15 already ships a transformer
tagger then removed it entirely. See §NanoAODv15 at the bottom — that section is
the decision, everything above it is the evidence trail that led there.

## The question

If `FatJet_particleNetWithMass_WvsQCD` is already stored in every NanoAOD and is
already excellent, there is no point building another W tagger. So: measure how
good it actually is, and find out whether anything is left uncovered.

## What is already measured (v2c + ML/results/summary.json)

QCD rejection (1/mistag) at W efficiency 0.3 / 0.5 / 0.7, physics-weighted,
p_T ∈ (200,1200), m_SD ∈ [20,250]:

| tagger | 0.3 | 0.5 | 0.7 | notes |
|---|---|---|---|---|
| `pnet_wvsqcd` (with-mass) | **299** | **86** | **30** | stored |
| `part_nom` (ours) | 136 | 48 | 19 | |
| `bdt_nom` (ours) | 134 | 46 | 18 | |
| τ₂₁ + mass window | 68 | 23 | 9 | |
| `part_md` (ours) | **32** | **12** | **5.5** | JSD 0.003–0.032 |
| `bdt_md` (ours) | 31 | 13 | 5.9 | JSD 0.019–0.077 |
| `pnet_xqqvsqcd` (stored MD) | 15 | 5.4 | 2.5 | JSD 0.004–0.035 |
| N2DDT(5%) | 17 | 7.1 | 3.2 | |

Cross-checked two ways: the ROC recomputed from the v2c `jets` tree gives
rejection 85.7 for `pnet_wvsqcd` and 5.46 for `pnet_xqqvsqcd`, against 86.3 and
5.43 in the earlier `ML/results/summary.json`. Independent pipelines agree.

**Conclusion so far.** As a general W tagger the project is dominated — 1.6–2.2×
worse at every working point, and structurally so: ParticleNet reads PF
candidates, which NanoAODv13 does not store ([[wtag-v2-status]] verified there is
no PFCands collection). That gap does not close with better features.
**As a mass-decorrelated tagger we win ~2×** at equal sculpting. That is the only
surviving justification, and it rests on one number that has not been checked.

## The load-bearing doubt: we may be benchmarking the wrong head

`XqqVsQCD` is the **light-quark** head. Roughly half of hadronic W decays contain
charm (W→cs), so a large fraction of real W jets are not what that head was
built to select. This is not speculation — measured single-variable rejection at
ε_S = 0.5 inside the mass window:

| head | rejection |
|---|---|
| `PNet_xccvsqcd` | **6.5** |
| `PNet_xqqvsqcd` | 5.4 |

W jets demonstrably populate the Xcc head. So the fair mass-decorrelated
benchmark is the **combined** discriminant, not `Xqq` alone:

```
X = QCD · r / (1 − r)        for each stored r = X/(X + QCD)
score = (Xqq + Xcc + Xbb) / (Xqq + Xcc + Xbb + QCD)
```

**This cannot be computed from what is stored.** The `jets` tree carries only
`wvsqcd`, `tvsqcd`, `xqqvsqcd`, `xggvsqcd` — no Xcc, no raw QCD probability, so
the combination cannot be formed per jet. Hence the rerun below.

## What the reference paper does and does not say

`ref/2604.09809v1-3.pdf` = CMS-JME-25-001, "Particle transformers for
identifying Lorentz-boosted Higgs bosons decaying to a pair of W bosons".

- Its task is **H→WW\*→4q**, i.e. 3- and 4-prong jets — not 2-prong W-vs-QCD.
- The algorithm it beats is **DeepAK8-MD**, which it calls "previously
  state-of-the-art in CMS for H→WW jet identification" (a 2D-image CNN, older
  than ParticleNet), by up to 10× in QCD rejection.
- **ParticleNet is never claimed to be beaten at tagging.** It appears only for
  mass regression (Fig. 4) and inference speed (6% slower than PART).
- Its stated reasons to exist: the H→WW\*→4q topology "has seldom been explored
  at the LHC"; no SM proxy existed to calibrate such a tagger (hence the new
  Lund-jet-plane method); and it enables the first all-hadronic HH→bbWW search.
- Inputs: up to 128 PF candidates + 10 SVs, 2.3M parameters, 37 output classes.

**So the paper does not license "ParticleNet is beatable at W-vs-QCD."** It is
evidence that new taggers get built when a *topology is uncovered*, not when an
existing tagger is merely improvable. If anything, 128 PF candidates is exactly
the input NanoAOD lacks, so it reinforces that a NanoAOD-scale tagger will not
win on raw rejection.

It does, however, hand over the correct framing to test: *is there a
W→cs-inclusive mass-decorrelated score in NanoAOD at all?* That is a gap
argument, not a "we are better" argument, and it is precisely what §Steps
measures.

## Surviving justifications, ranked

1. **A gap in what is stored** — if no combination of stored MD heads covers
   W→cs, a NanoAOD-scale MD W tagger fills a real hole. *Decided by step 3.*
2. **Accessibility** — a NanoAOD-only tagger runs on CMS Open Data and anywhere
   PFNano is unavailable. True regardless, but weak on its own.
3. **Quantifying the cost of NanoAOD scale** — "what do you lose with only
   high-level features" is a number nobody has published; currently ~1.8× at
   ε_S = 0.5. Publishable as a note, not as a tagger.

If step 3 shows the combined stored MD score reaches `part_md` (≈12 at
ε_S = 0.5), justification 1 dies and only 2–3 remain — at which point the
tagger should be dropped and the work rewritten as a benchmarking study.

## Steps

1. **DONE** — `Wtag.cc` gained `pnet_qcd` / `pnet_xccvsqcd` / `pnet_xbbvsqcd`
   on the `jets` tree, guarded `Run == 3` like the histogram block (Run-2 falls
   back to the −10 "not evaluated" sentinel NanoAOD itself uses).
   `ML/make_dataset.py` gained `tau31`/`tau41` as **features** (added to
   `FEATURES`, to the derive-later skip list, and to the ratio loop) and the
   three new ParticleNet columns as **`EXTRAS`, deliberately not features** — a
   model handed the benchmark would only learn to copy it. `pnet_qcd` is
   excluded from the ≥0 sentinel clamp: it is a raw probability with no
   sentinel, and clamping it would corrupt the `X = QCD·r/(1−r)` inversion the
   whole benchmark rests on. **No training run.**
2. **DONE** — rebuilt (`rebuild.sh`, exit 0; all five new literals verified
   present in `libAnalyzers.so`) and submitted as tag `v2d` at 14:11:
   26 DAGs / 999 nodes, zero empty, 586 jobs queued, 0 held at submit.
   Expect the 8 GB hadd limit (hardcoded, `python/SKNano.py:447`) to kill the
   same ~4 big samples again — `merge_v2c.py`'s `rescue_map()` handles it.
3. **Measure the fair MD benchmark.** ROC the combined
   `(Xqq+Xcc+Xbb)/(…+QCD)` score against `part_md`'s 32 / 12 / 5.5, in the same
   weighting and selection, with the JSD sculpting check alongside so the
   comparison is at equal decorrelation.

## Decision rule

| combined stored MD rejection @ ε_S = 0.5 | verdict |
|---|---|
| ≲ 8 | gap is real — continue the MD tagger |
| 8–12 | marginal — continue only with a sharper physics case |
| ≳ 12 | no gap — drop the tagger, rewrite as a benchmarking note |

## RESULT (2026-07-26, v2d)

870M jets, all 26 samples clean, every hadd finished (the 5 held nodes were
fixed with `condor_qedit RequestMemory 24576` + `condor_release` rather than
rescued — `MemoryUsage` 9766 MB confirms 8192 was never enough).

| discriminant | @0.3 | @0.5 | @0.7 | JSD @ ε_b = 20/10/5/2/1% |
|---|---|---|---|---|
| `xqq` only (stored) | 15.2 | 5.5 | 2.5 | 0.006 0.008 0.013 0.024 0.034 |
| **`xqq+xcc` (stored)** | **20.8** | **8.4** | **4.1** | 0.005 0.007 0.009 0.015 0.020 |
| `xqq+xcc+xbb` (stored) | 18.4 | 7.7 | 3.9 | 0.005 0.005 0.006 0.009 0.010 |
| `wvsqcd` with-mass (stored) | 291.7 | 85.7 | 29.4 | 0.257 … 0.495 |
| `part_md` (ours) | 32 | 12 | 5.5 | 0.003 0.006 0.011 0.021 0.032 |

Validation: `xqq` reproduces `summary.json` (15.2 vs 15, 5.5 vs 5.4) and
`wvsqcd` likewise (85.7 vs 86), so the tree-level selection matches the ML
evaluation and the rows are directly comparable. Closure on the inversion:
99.6% of jets have Xqq+Xcc+Xbb+Q ≤ 1 (max 1.207); the small excess is the
tree's ~11-bit float truncation amplified by r/(1−r) as r→1.

**The suspicion was correct — benchmarking against `xqq` alone understated
ParticleNet by 55%.** Adding the charm head takes 5.5 → 8.4 at ε_S = 0.5.
Adding `xbb` on top *hurts* (7.7), as it should: bb is not a W-like final state
and only dilutes.

**Verdict: marginal, at the bottom edge of the middle band.** Our margin over
the best stored MD combination is 1.43× at ε_S = 0.5, not the 2.2× that `xqq`
alone suggested.

**And 1.43× is an overstatement**, because the two are not at equal sculpting:
at ε_b = 1% the stored combination has JSD **0.020** against `part_md`'s
**0.032**, i.e. the stored score is *more* mass-decorrelated. A like-for-like
comparison would cost `part_md` some rejection, and the true margin is somewhere
below 1.43×. Settling that needs `part_md` re-evaluated at matched JSD, which is
an ML job and deliberately out of scope here.

### What this means for the project

Justification 1 (a real gap in what is stored) is **not supported**. There is a
W→cs-inclusive mass-decorrelated discriminant available in NanoAOD today; it
just has to be assembled from two heads instead of read off one. A ~1.3× edge,
at worse decorrelation, is not a reason to build and calibrate a tagger.

Justifications 2 and 3 survive intact, and 3 is now sharper than before: the
cost of NanoAOD scale is **1.8×** in the with-mass regime (48 vs 86) and
**~1.3–1.4×** in the mass-decorrelated one (12 vs 8.4, before the sculpting
correction, i.e. NanoAOD-scale is *ahead* here). That asymmetry is itself the
interesting result — high-level features lose badly when mass is available and
nearly break even once mass is taken away, which says the PF-level information
ParticleNet exploits is largely mass-correlated substructure.

Recommended next move: stop developing the tagger as a product, and write this
up as a benchmarking note answering "what does a NanoAOD-only W tagger cost
you?". The one experiment that could still revive justification 1 is a matched-
JSD comparison, since if `part_md` holds 12 at JSD 0.020 the margin is real.

## Explicitly not in this phase

- No retraining, no GPU-server work. `ML/models/` and `ML/data/dataset.parquet`
  stay as they are; step 1 only edits the feature *list* for a later run.
- No attempt to beat `pnet_wvsqcd` with-mass. That comparison is settled.

## NanoAODv15 closes it (2026-07-26) — DECISION: stop

Step 3 left the project alive but marginal. It is not marginal any more: the
premises it rested on are obsolete in NanoAODv15, which is already in the
central store (`/gv0/DATA/SKNano/NanoAODv15/2024/`, 106 datasets, 2024 only —
v13p1 still covers 2022–2023BPix) and already wired up on `upstream/NanoAODv15`
(`DataFormats/include/FatJet.h`, `SetGloParTResult`).

**GloParT3 is stored per fat jet** — 18 mass-decorrelated heads + 3 with-mass:

```
QCD, TopbWev, TopbWmv, TopbWq, TopbWqq, TopbWtauhv,
WvsQCD, XWW3q, XWW4q, XWWqqev, XWWqqmv,
Xbb, Xcc, Xcs, Xqq, Xtauhtaue, Xtauhtauh, Xtauhtaum
+ withMass{Top,W,Z}vsQCD, massCorrGeneric, massCorrX2p
```

Verified populated in `TTto4Q .../251121_123928/0000/tree_0.root` (48 FatJet
tagger branches; 136,507 leading fat jets in the file):

| branch | min | max | mean |
|---|---|---|---|
| `FatJet_globalParT3_Xcs` | 0.0000 | 0.9541 | 0.0138 |
| `FatJet_globalParT3_Xqq` | 0.0000 | 1.7441 | 0.0885 |
| `FatJet_globalParT3_WvsQCD` | −1.0000 | 0.9995 | 0.1476 |
| `FatJet_particleNetWithMass_WvsQCD` | −10.0000 | 0.9995 | 0.1915 |

(Note the sentinel differs by family: GloParT3 uses −1, ParticleNet −10. And
`Xqq` exceeding 1 says these are not a normalized softmax over the 18 heads —
do not assume Σ = 1 if these are ever used.)

### Why this kills each justification

1. **A gap in what is stored — dead, twice.**
   `Xcs` is a dedicated W→cs head, which was the entire gap argument: §The
   load-bearing doubt was that no *stored* head covers charm, so W→cs had to be
   recovered by hand. It is now trained in. Independently, v15 also adds
   **`FatJet_particleNet_WVsQCD`** (capital V) — a mass-decorrelated ParticleNet
   *W* head, which v13 lacked and which is exactly why step 3 had to assemble
   `Xqq+Xcc` by hand to reach 8.4. A purpose-trained W head will beat a hand-sum
   of two heads never calibrated for that sum, so the real stored MD benchmark
   is now *above* 8.4 and our 12 is not a margin worth defending.
2. **Accessibility (NanoAOD-only) — dead, and this was the sturdier one.**
   The premise was "a transformer needs PF candidates and NanoAOD does not store
   them" ([[wtag-v2-status]] verified that for v13). v15 does not store the
   candidates for the tagger's sake — it stores the tagger's **output**. A
   NanoAOD-scale user now gets PF-level performance with no PF-level inputs, so
   there is nothing our tagger uniquely provides.
3. **Quantifying the cost of NanoAOD scale — survives, but demoted.**
   Still measurable (1.8× with mass, ~1.3–1.4× without, §RESULT), but it is now
   a curiosity rather than a gap: the answer for any actual user is "read the
   branch."

`XWW3q`/`XWW4q` in that head list are the H→WW\*→4q classes of
`ref/2604.09809v1-3.pdf`. GloParT3 is essentially that paper's tagger, shipped
in NanoAOD — which also retires §What the reference paper does and does not say
as a live question.

### v15 does store PF candidates — and it still is not a pivot

`nPFCand`, `PFCand_{pt,eta,phi,mass,pdgId}`, `FatJetPFCand_{jetIdx,pfCandIdx}`;
⟨nPFCand⟩ = 15.0 per event at ⟨nFatJet⟩ = 0.40, i.e. ~37 constituents per AK8
jet. So "true ParT impossible without PFNano" is false for v15 and that line in
[[wtag-v2-status]] is now version-scoped.

It is not an opening. Five kinematic features per candidate is not GloParT's
input: no track impact parameters, no PUPPI weights, no charge, no per-candidate
quality — the displaced-track information that carries heavy-flavour separation
is absent. A reduced-input transformer is trainable on v15; GloParT3 is not
reproducible from it, and it would compete against GloParT3's output sitting in
the same file. Same structural loss as before, one level up.

### What is salvageable

Not a tagger. Two things, in order:

1. **τ₃/τ₁ and τ₄/τ₁ beat the conventional τ₂/τ₁** — 18.0 / 15.9 vs 14.4 at
   ε_S = 0.5 inside the mass window, and both beat `msoftdrop`. Independent of
   taggers and of NanoAOD version; the reason is clean (the denominator sets the
   question — τ_N/τ₁ asks "better than one axis", which *is* the W-vs-QCD
   question, while τ₃/τ₂ and τ₄/τ₃ ask the top-tagger question). Details in
   [[wtag-v2b-branch-comparison]] / `PLAN_branch_compare.md`. Short note.
2. **The §RESULT benchmark table** as supporting material — a clean, finished
   measurement of what high-level features cost against a PF-level tagger.

### Not done, deliberately

No v2e. No matched-JSD re-evaluation of `part_md` — it was the one experiment
that could have revived justification 1, and it is moot now that a trained W
head exists in v15. `ML/models/`, `ML/data/dataset.parquet` and the v2d outputs
on /gv0 are left intact. The `tau31`/`tau41` feature-list edits in
`ML/make_dataset.py` stay: they cost nothing and they are what salvage item 1
would want if it is ever written up.
