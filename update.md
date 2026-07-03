# Review Note — SKNanoLoader: TTreeReader migration + Run3 branch-type bug fixes

For the SKNanoAnalyzer developer / server manager reviewing this version.
This version contains exactly **two things**: (1) the removal of the
17-pass startup scan by migrating branch reading to `TTreeReader`, and
(2) fixes for Run3 branches that were silently read as garbage because of
integer-type mismatches. Nothing else: no output-format change, no new
user-facing features, no behavior knobs.

## Scope — files touched

| File | What |
|---|---|
| `Analyzers/src/SKNanoLoader.cc` | TTreeReader read path; legacy path kept for skimming; type-mismatch fixes in both paths |
| `Analyzers/include/SKNanoLoader.h` | Members for the above (`fReader`, filler vectors, `SkimmingMode`, raw `Buf_*` buffers) |
| `scripts/compare_outputs.C` | New: standalone macro that diffs two output ROOT files (histograms + trees) and prints IDENTICAL/DIFFER — used for the verification below |

**No analyzer code was touched.** The member interface the analyzers use
(`Muon_pt[i]`, `nJet`, `TriggerMap`, …) is unchanged; every existing
analyzer compiles and runs as-is.

## 1. The problem: startup did 17 full passes over the input

Old `SKNanoLoader::SetMaxLeafSize()` sized the `RVec` members by running
**17 sequential `RDataFrame` event loops over the whole chain at startup**
(one per collection counter) before a single event was analyzed. Measured
cost: **9.8 s for a 1-file job, 15.4 s for a 3-file chain** — pure waste,
paid by every condor job, and 17× redundant read traffic on `/gv0` at
job start.

## 2. Change 1: branch reading via TTreeReader

- `Init()` now dispatches to `InitTTreeReader()` (default) or
  `InitLegacy()` (the previous code, kept intact).
- `InitTTreeReader()` looks up each branch's **on-file leaf type**, creates
  a typed `TTreeReaderValue`/`TTreeReaderArray`, and registers a small
  "filler" lambda that copies the proxy into the existing class member.
  Runtime dispatch covers 11 leaf types, so the same binary handles
  NanoAOD v9/v12/v13 type differences without recompiling.
- `Loop()` per event: scalar fillers → existing Run2/Run3 counter-sync
  block (unchanged) → array fillers (they use the synced `nX` counters to
  reproduce the absent-branch zero-fill semantics of the old code).
- **No startup pre-scan is needed at all**: init drops from 9.8 s / 15.4 s
  to ~0.04 s. The event-loop rate itself is unchanged (same decompression
  work). TTreeReader also only ever reads branches that have a proxy, so
  the old `SetBranchStatus` bookkeeping is automatic.

### Why skimmers keep the legacy path (`SkimmingMode`)

`Skim_*` analyzers do `fChain->CloneTree(0)` + `newtree->Fill()`, which
requires the `SetBranchAddress` buffers that `fChain->GetEntry()` fills;
TTreeReader does not fill those. `Init()` therefore auto-sets
`SkimmingMode = true` when the class name contains "Skim" (it can also be
set manually before `Init()`). The legacy path's `SetMaxLeafSize()` was
also fixed to take array maxima from `TLeaf::GetMaximum()` metadata, so
even skim jobs no longer pay the 17-pass startup.

## 3. Change 2: Run3 branch-type bug fixes (the only output changes)

The old code bound every branch by the *member's* type. NanoAOD v12+
narrowed several branches (e.g. `Int_t` → `UChar_t`), and
`SetBranchAddress` byte-copies with no conversion or warning, so those
members were silently garbage. Concrete example: `Muon_nTrackerLayers` is
`UChar_t` on file; read into an `Int_t` array it yields values like
`3344` (= 16 + 13·256, two muons' bytes packed into one int) — **and this
member feeds the Rochester muon momentum correction**, so Run3 RoccoR
inputs were corrupted.

Fixed branches (Run3 values change, correctly):

| Branch | Old bug |
|---|---|
| `Muon_nTrackerLayers` | UChar→Int byte-pack garbage → fed RoccoR |
| `Jet_chMultiplicity`, `Jet_neMultiplicity` | same UChar→Int pattern |
| `GenVisTau_charge`, `GenVisTau_status`, `GenVisTau_genPartIdxMother` | Short/UChar→Int byte-pack |
| `FatJet_genJetAK8Idx` (and `_RunII`) | bound to a never-resized RVec → undefined behavior |

How the fix works in each path:

- **TTreeReader path**: reads the on-file type directly and `static_cast`s
  into the member — correct by construction.
- **Legacy (SkimmingMode) path**: Run3 binds these branches to raw
  `Buf_*` vectors of the on-file type and widens them into the members at
  the top of `Loop()` (bounds-capped). Run2 bindings are untouched.

Same-width mismatches (`run`/`event` UInt↔Int, `nLHEPdfWeight`) are
value-preserving and were deliberately left alone to keep the diff minimal.

## 4. Verification

All checks compare this version against the **unmodified original code**
compiled side-by-side into a standalone dump harness (no framework
changes involved), on real files:
2022 = Run3NanoAODv13p1 DYto2L-2Jets_MLL-50, 2016postVFP = Run2NanoAODv9p1
DYJetsToLL_M-10to50.

| Check | Result |
|---|---|
| All ~420 members + trigger map, original vs new, 2016postVFP (Run2), 2000 events | **bit-identical** |
| Same, 2022 (Run3), 2000 events | identical except the §3 bug-fixed branches (by design) |
| New code: TTreeReader path vs legacy (SkimmingMode) path, 2022, 2000 events | **identical** — the two reading modes agree, including the fixed branches |
| Analyzer-style histograms (lepton/jet kinematics, Σweights, HLT bit; 10k events), original vs new | **IDENTICAL** for 2022 and 2016postVFP (`scripts/compare_outputs.C`) |
| Skim output (`CloneTree(0)`+`Fill`, 5k events), original vs new | tree content **IDENTICAL** |
| `Muon_nTrackerLayers` after fix | physical values (11, 14, 12, …) instead of byte-packed garbage (3344, …) |

The only unrelated diff ever observed: `nTau_RunII`/`nJet_RunII`-type
scalars that are *never written in Run3* and hold uninitialized stack
values — pre-existing, harmless, unchanged by this work.

Suggested reviewer check after building: run one Run2 analyzer
before/after and compare outputs with
`root -l -b -q 'scripts/compare_outputs.C("old.root","new.root")'`
(expect IDENTICAL), and one Run3 analyzer using muons (expect differences
only through RoccoR ← `Muon_nTrackerLayers`).

## 5. Known behavior differences / risk notes

- A branch present in the first chain file but **missing in a later file**
  now stops the job with an error under TTreeReader; the old code
  silently kept stale values from the last file that had it. Failing
  loudly was chosen on purpose. (DATA trigger branches keep the old
  per-file existence check and read as `false` where missing.)
- `SkimmingMode` auto-detection relies on the `Skim*` class-name
  convention; a skimmer named otherwise must set `SkimmingMode = true`
  before `Init()`.
- Run3 outputs that consume the §3 branches (notably muon Rochester
  corrections) **will change** relative to previous productions — that is
  the bug fix, not a regression.
