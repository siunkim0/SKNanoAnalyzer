#!/usr/bin/env python3
"""Train ParT-lite (nominal + mass-decorrelated) on the GPU server.

  python3 train_partlite.py            # both variants
  python3 train_partlite.py nom|md     # one variant

Tokens: subjet1, subjet2, up to 3 in-cone SVs (dlenSig-ordered).
Globals: the jet-level HLF vector (variant-dependent, from meta.json,
minus the per-subjet columns which live in the tokens).
Outputs models/part_<variant>.pt (weights + normalization + config).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset

from models import ParTLite

BASE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(BASE, "models"), exist_ok=True)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# columns that live in tokens, not in the global vector
SUBJET_COLS = ["sj1_pt", "sj1_mass", "sj1_btag", "sj2_pt", "sj2_mass", "sj2_btag"]
EPOCHS, BATCH, LR, WD = 30, 1024, 1e-3, 1e-2
N_TOKEN_FEAT = 6  # [log pt, mass', extra1, extra2, is_sj, is_sv]


def build_tokens(df, md):
    """(N, 5, 6) float32 token array + (N, 5) bool mask."""
    n = len(df)
    tok = np.zeros((n, 5, N_TOKEN_FEAT), dtype=np.float32)
    mask = np.zeros((n, 5), dtype=bool)

    def logpt(x):
        return np.where(x > 0, np.log(np.maximum(x, 1e-3)), 0.0)

    for i, p in enumerate(("sj1", "sj2")):
        pt = df[f"{p}_pt"].values
        ok = pt > 0
        m = np.zeros(n) if md else df[f"{p}_mass"].values / 50.0
        tok[:, i, 0] = logpt(pt)
        tok[:, i, 1] = np.where(ok, m, 0.0)
        tok[:, i, 2] = np.where(ok, df[f"{p}_btag"].values, 0.0)
        tok[:, i, 3] = 0.0 if md else np.where(ok, df["sj_dr"].values, 0.0)
        tok[:, i, 4] = 1.0
        mask[:, i] = ok
    for k in (1, 2, 3):
        i = 1 + k
        pt = df[f"sv{k}_pt"].values
        ok = pt > 0
        tok[:, i, 0] = logpt(pt)
        tok[:, i, 1] = np.where(ok, df[f"sv{k}_mass"].values / 5.0, 0.0)
        tok[:, i, 2] = np.where(ok, df[f"sv{k}_dlensig"].values / 50.0, 0.0)
        tok[:, i, 3] = np.where(ok, df[f"sv{k}_ntracks"].values / 10.0, 0.0)
        tok[:, i, 5] = 1.0
        mask[:, i] = ok
    tok[:, :, :4] = np.clip(tok[:, :, :4], -5.0, 20.0)
    return np.nan_to_num(tok, nan=0.0, posinf=0.0, neginf=0.0), mask


def global_features(meta, variant):
    feats = meta["features"] if variant == "nom" else meta["features_md"]
    return [f for f in feats if f not in SUBJET_COLS]


def make_loader(df, gfeats, md, norm, wcol, shuffle):
    tok, mask = build_tokens(df, md)
    glob = ((df[gfeats].values - norm["mean"]) / norm["std"]).astype(np.float32)
    glob = np.nan_to_num(glob, nan=0.0, posinf=0.0, neginf=0.0)
    ds = TensorDataset(
        torch.from_numpy(tok), torch.from_numpy(mask),
        torch.from_numpy(glob),
        torch.from_numpy(df.y.values.astype(np.float32)),
        torch.from_numpy(df[wcol].values.astype(np.float32)),
    )
    return DataLoader(ds, batch_size=BATCH, shuffle=shuffle, num_workers=4,
                      pin_memory=(DEVICE == "cuda"), drop_last=shuffle)


def run_epoch(model, loader, opt=None):
    training = opt is not None
    model.train(training)
    tot, wsum = 0.0, 0.0
    with torch.set_grad_enabled(training):
        for tok, mask, glob, y, w in loader:
            tok, mask, glob = tok.to(DEVICE), mask.to(DEVICE), glob.to(DEVICE)
            y, w = y.to(DEVICE), w.to(DEVICE)
            logit = model(tok, mask, glob)
            loss = (torch.nn.functional.binary_cross_entropy_with_logits(
                logit, y, reduction="none") * w).sum() / w.sum()
            if training:
                opt.zero_grad()
                loss.backward()
                opt.step()
            tot += loss.item() * w.sum().item()
            wsum += w.sum().item()
    return tot / wsum


def train(variant):
    md = variant == "md"
    meta = json.load(open(os.path.join(BASE, "data", "meta.json")))
    gfeats = global_features(meta, variant)
    wcol = "w_flatpt" if variant == "nom" else "w_flatptm"

    df = pd.read_parquet(os.path.join(BASE, "data", "dataset.parquet"))
    df = df[df.y >= 0]          # drop t->bqq spectator rows
    tr, va = df[df.split == "train"], df[df.split == "val"]
    vals = tr[gfeats].values
    norm = {"mean": np.nanmean(vals, 0), "std": np.nanstd(vals, 0) + 1e-6}

    model = ParTLite(N_TOKEN_FEAT, len(gfeats)).to(DEVICE)
    print(f"[{variant}] params={sum(p.numel() for p in model.parameters()):,} "
          f"globals={len(gfeats)} device={DEVICE}")
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS)

    tr_loader = make_loader(tr, gfeats, md, norm, wcol, shuffle=True)
    va_loader = make_loader(va, gfeats, md, norm, wcol, shuffle=False)

    best, best_state = np.inf, None
    for ep in range(EPOCHS):
        tl = run_epoch(model, tr_loader, opt)
        vl = run_epoch(model, va_loader)
        sched.step()
        star = vl < best
        if star:
            best = vl
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        print(f"[{variant}] epoch {ep+1:02d}/{EPOCHS} train={tl:.4f} val={vl:.4f}{' *' if star else ''}")

    out = os.path.join(BASE, "models", f"part_{variant}.pt")
    torch.save({"state_dict": best_state, "gfeats": gfeats, "md": md,
                "norm_mean": norm["mean"], "norm_std": norm["std"],
                "n_token_feat": N_TOKEN_FEAT, "val_loss": best}, out)
    print(f"[{variant}] best val loss {best:.4f} -> {out}")


if __name__ == "__main__":
    for v in (sys.argv[1:] or ["nom", "md"]):
        train(v)
