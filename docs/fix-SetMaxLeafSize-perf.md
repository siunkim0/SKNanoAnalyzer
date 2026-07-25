# Task: Migrate `SKNanoLoader` reading to `TTreeReader` (performance refactor)

> Paste this whole file to Claude Code as the task prompt. It is self-contained.
> Reply in Korean (code/commands/paths stay English). Do NOT run builds or
> `SKNano.py` yourself — the user builds and runs jobs. For ROOT on the CLI, prefix
> with `rot &&` (e.g. `rot && root -l -b -q macro.C`); otherwise ROOT is not on PATH.

## Goal

This is a **pure performance refactor**, not a behavior change. Replace the manual
`SetBranchAddress` + `SetMaxLeafSize()` machinery with ROOT's `TTreeReader` /
`TTreeReaderValue` / `TTreeReaderArray`, which read variable-length arrays into a
self-growing internal buffer. This removes the need to precompute array sizes and
eliminates the startup bottleneck entirely (see below).

### ⚠️ Non-negotiable equivalence contract (read first)

**After this change the framework MUST behave EXACTLY as before. The analysis output
must be identical — same histograms, same trees, same event counts, same weights,
bit-for-bit where floating point allows.** Nothing about selection, corrections,
systematics, weights, or output format may change. This is *only* a change in HOW
branches are read from the input files, not WHAT is computed from them.

If any design choice risks changing output, do NOT make it — stop and flag it.
The task is done only when a before/after output comparison (see "Verification")
shows identical results.

## Why (the bottleneck being removed)

`SetMaxLeafSize()` in `Analyzers/src/SKNanoLoader.cc` (~line 76, called from `Init()`
~line 673) computes each array's max multiplicity with:

```cpp
auto maxValue = static_cast<int>(*(df.Max(branchName)));   // immediate dereference
```

`RDataFrame::Max()` is lazy, but dereferencing each result immediately triggers a
separate event loop. With 17 counter branches this does **17 full sequential passes
over the entire TChain before analysis starts**, each JIT-compiled, single-threaded.
`TTreeReader` needs none of this: no max, no fixed buffers, no pre-scan.

## Current design (what you are replacing)

- Members in `Analyzers/include/SKNanoLoader.h` are the public interface analyzers
  read, e.g. `RVec<Float_t> Jet_pt;`, `RVec<Int_t> GenPart_pdgId;`, scalars like
  `Int_t nJet;`, `Float_t genWeight;`.
- `SetMaxLeafSize()` resizes every `RVec` to the measured max.
- `Init()` binds each branch: arrays via `SetBranchAddress(name, vec.data())`,
  scalars via `SetBranchAddress(name, &scalar)`; helpers `SafeSetBranchAddress`
  (skips branch if absent), `SuperSafeSetBranchAddress` (skips if absent in ANY file,
  used for DATA triggers), `SetBranchWithRunCheck` (binds Run3 `Int_t` var or Run2
  `UInt_t` `_RunII` var depending on `Run`).
- `Loop()` does `fChain->GetEntry(jentry)`, then a Run2 sync block copying
  `nJet = nJet_RunII;` etc., then `executeEvent()`.
- Triggers: `TriggerMap[path].first` is a `Bool_t*` bound to the HLT branch; special
  cases: path `"Full"` is forced true; paths missing from the tree are dropped.

## Target design (keep the member interface; drive it from TTreeReader)

**Keep every existing data member exactly as-is** (same names, same types). Analyzers
must not need any change — that is what guarantees identical output. Only change how
those members get filled each event.

1. Add a reader and one reader-proxy per branch to the loader:
   ```cpp
   std::unique_ptr<TTreeReader> fReader;
   std::unique_ptr<TTreeReaderArray<Float_t>> rdr_Jet_pt;
   std::unique_ptr<TTreeReaderValue<Int_t>>   rdr_nJet;
   // ... one per branch, matching the branch's on-file type
   ```

2. In `Init()`, replace `SetMaxLeafSize()`, `SetBranchStatus("*",0)`, and all the
   `SetBranchAddress` calls with reader construction. Create a proxy ONLY if the
   branch exists, so absent branches keep today's behavior (member left default):
   ```cpp
   fReader = std::make_unique<TTreeReader>(fChain);
   auto makeArr = [&](auto& ptr, const char* name){
       using Proxy = typename std::remove_reference_t<decltype(ptr)>::element_type;
       if (fChain->GetBranch(name)) ptr = std::make_unique<Proxy>(*fReader, name);
   };
   makeArr(rdr_Jet_pt, "Jet_pt");
   // scalars: analogous with TTreeReaderValue
   ```
   (TTreeReader only reads requested branches, so the old `SetBranchStatus("*",0)`
   optimization is automatic — do not re-add it.)

3. In `Loop()`, replace `fChain->GetEntry(jentry)` with `fReader->SetEntry(jentry)`
   and copy each proxy into its existing member; keep the Run2 sync block and
   `executeEvent()` unchanged:
   ```cpp
   if (fReader->SetEntry(jentry) != TTreeReader::kEntryValid) { /* same error path as today */ }
   // arrays -> existing RVec members (empty if branch absent, matching today):
   if (rdr_Jet_pt) Jet_pt = RVec<Float_t>(rdr_Jet_pt->begin(), rdr_Jet_pt->end());
   else            Jet_pt.clear();
   // scalars:
   nJet = rdr_nJet ? **rdr_nJet : 0;
   ```
   Copying into the same members preserves the exact interface the analyzers use;
   the removed 17-pass scan is the real win, so the per-event copy is an acceptable
   and negligible cost. Keep `MaxEvent` / `NSkipEvent` logic identical.

### Cases that must be preserved 1:1

- **Run2/Run3 dual types**: where `SetBranchWithRunCheck` / the `if(Run==3)/else`
  blocks read a branch into a different member/type per run (e.g.
  `GenPart_genPartIdxMother` (Int_t) vs `GenPart_genPartIdxMother_RunII` (UInt_t),
  the `nX` vs `nX_RunII` counters), replicate the exact same per-run choice: create
  the proxy with the correct on-file type and fill the same member the old code did.
  Keep the existing Run2 sync block in `Loop()` as-is.
- **Missing branches**: `SafeSetBranchAddress` silently skips absent branches; a null
  proxy + default/empty member must reproduce that (no new warnings-as-errors, no
  crash).
- **DATA triggers / `SuperSafeSetBranchAddress`**: preserve "drop trigger if missing
  in any file" and the per-file check, the `"Full"` = true special case, and the
  dropped-path warning list. `TriggerMap[path].first` (`Bool_t*`) can be filled from
  a `TTreeReaderValue<Bool_t>` each event, or keep its current handling — whichever
  keeps output identical.
- The `if(fChain->GetEntries()==0) exit(0)` guard and all timing/log prints' intent.

### Things to delete

- `SetMaxLeafSize()` and its declaration in the header.
- All `.resize(kMax...)` sizing (the reader owns buffers now).
- `SetBranchAddress` calls and the now-unused Safe/SuperSafe/RunCheck address
  binding (keep the branch-existence and per-file trigger checks if you still need
  them for proxy creation).

## Do NOT change

- Any file under `Analyzers/src/*.cc` analyzers, `DataFormats/`, `AnalyzerTools/`,
  correction/weight/systematic logic, YAML/JSON configs, output file structure.
- Member names/types in `SKNanoLoader.h` that analyzers consume.
- Event selection, ordering, or which events are written.

## Verification (required — proves output is unchanged)

The user runs the jobs; you provide the recipe and the comparison macro.

1. **Before** the change, on a small fixed input (use `--reduction` so it is quick),
   run one representative analyzer for both a Run3 and a Run2 era, and keep the
   output ROOT files (e.g. `before_2022.root`, `before_2016postVFP.root`).
2. **After** the change, rebuild and run the identical commands → `after_*.root`.
3. Compare. Provide a ROOT macro that walks every histogram and TTree and asserts
   equality, e.g. compares histogram bin contents/entries and tree branch values;
   any nonzero diff = fail. `rot && root -l -b -q compare.C` should print `IDENTICAL`.
4. Also confirm the log no longer shows the 17-pass `SetMaxLeafSize` timing and that
   startup is visibly faster.

State explicitly in your summary that output was verified identical (or, if you could
not run it, that the user must run this comparison before trusting the change).

## Suggested approach

- This touches hundreds of branches — do it mechanically and completely; a partial
  migration that mixes `TTreeReader` and `SetBranchAddress` on the same `fChain` will
  not work. Migrate all branches in one consistent pass.
- If the full migration is judged too risky in one shot, a smaller drop-in that
  removes the SAME 17-pass bottleneck without touching the read path is to replace
  only the measurement in `SetMaxLeafSize()` with leaf-metadata lookups
  (`TLeaf::GetMaximum()` over `fChain->GetListOfFiles()`, no event loop) — verified
  to return the correct max from metadata alone. Mention this as the fallback but
  prefer the `TTreeReader` migration the user asked for.

## Constraints

- Do not run `./scripts/build.sh`, `./scripts/rebuild.sh`, or `SKNano.py` — the user
  builds and submits. Make the code change and give exact build + verification steps.
- Tree name is `"Events"`; prefer `fChain->GetName()` over a literal.
- ROOT 6.30.02. Reply in Korean.
