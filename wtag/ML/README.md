# Wtag v2 — ML training package (runs on the GPU server)

Pipeline (see `../PLAN_v2.md`):

1. **On tamsa2** (needs ROOT/uproot — LCG el9 view):
   `python3 make_dataset.py` reads the merged v2 ntuples in `../2022_v2/*.root`
   (tree `jets`) and writes `data/dataset.parquet` + `data/meta.json`
   (features, splits, flat-pT and (pT,m)-flat weights precomputed).
2. **Copy `ML/` (with `data/`) to the GPU server**, `pip install -r requirements.txt`.
3. `python3 train_bdt.py` → `models/bdt_nom.json`, `models/bdt_md.json`
4. `python3 train_partlite.py` → `models/part_nom.pt`, `models/part_md.pt`
5. `python3 evaluate.py` → `results/` (ROC + rejection table, eff-vs-pT at fixed
   mistag, mass-sculpting JSD, feature importance). Baselines computed from the
   same test set: τ21, τ21+mass window, N2^DDT (quantile-map scan), stored
   ParticleNet WvsQCD (nominal) and XqqVsQCD (mass-decorrelated).

Classes: signal = `label==1` (clean W) from TTLJ/TTJJ/WW/WZ; background =
`label==0` from the QCD_PT samples only (TT/VV unmatched jets are spectators).
`label==2` (t→bqq̄ capture) is kept in the dataset for a later top-vs-W study.

Selection: AK8 pT ∈ [200, 1200] GeV, |η| < 2.4, tight ID (in ntuple), 20 < m_SD < 250.
Training weights: per-class flat-pT (nominal) or flat-(pT, m_SD) (MD variants),
so neither the spectrum nor (for MD) the mass is learnable.
