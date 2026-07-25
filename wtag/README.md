# Boosted W Tagger on NanoAOD high-level features

A machine-learning study of hadronic **W-jet tagging** for CMS Run 3 (2022 +
2022EE), built entirely from the high-level features stored in NanoAOD
(Run3NanoAODv13p1) — no PF-candidate branches required. Two taggers (an XGBoost
BDT and a small attention model, *ParT-lite*) are trained in both a **nominal**
and a **mass-decorrelated** variant, and benchmarked against the CMS/ATLAS
classic W tags (τ₂₁, τ₂₁ + mass window, N₂ᴰᴰᵀ) and the constituent-based
ParticleNet ceiling.

> This repository contains **code and documentation only**. ROOT ntuples,
> datasets, trained model weights, and figures/plots are produced locally and
> are intentionally excluded via `.gitignore` (see below).

## Layout

| path | contents |
|---|---|
| `PLAN.md`, `PLAN_v2.md` | analysis plan, method, and results write-up |
| `Tools/` | plotting and post-processing (`wtag_plots.py`), sample registration, hadd/merge scripts |
| `ML/` | training package — dataset builder, BDT & ParT-lite models, evaluation |
| `ML/results/` | result tables (rejection, working points) as text |
| `paper/` | LaTeX source of the write-up (`wtag_v2.tex`) |

## ML pipeline (`ML/`)

See `ML/README.md` and `PLAN_v2.md` for the full recipe. In short:

1. `make_dataset.py` — read merged v2 ntuples → `data/dataset.parquet` (+ meta).
2. `train_bdt.py` — XGBoost nominal & mass-decorrelated (`models/bdt_*.json`).
3. `train_partlite.py` — attention model nominal & MD (`models/part_*.pt`).
4. `evaluate.py` — ROC, rejection tables, efficiency-vs-pT, mass sculpting (JSD),
   feature importance, plus all cut-based / ParticleNet baselines.

`pip install -r ML/requirements.txt` (a GPU is used for the ParT-lite training).

## Key result

Nominal HLF taggers beat the τ₂₁ + mass window by ≈2× in background rejection at
ε_sig = 0.5; the mass-decorrelated variant matches N₂ᴰᴰᵀ / ParticleNet-MD
sculpting while keeping higher rejection. Full numbers are in the "Results"
section of `PLAN_v2.md` and in `ML/results/`.

## What is not tracked

`*.root`, `*.parquet`, model weights (`ML/models/`, `*.pt`), and all figures
(`*.pdf`, `*.png`, `Tools/Plots/`) are git-ignored. Regenerate them from the
scripts above.
