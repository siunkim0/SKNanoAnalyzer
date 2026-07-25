"""ParT-lite: a small attention tagger over the token set available in
NanoAOD — 2 softdrop subjets + up to N_SV in-cone secondary vertices — with
the jet-level high-level features injected through the class token
(class-attention pooling as in the Particle Transformer).

Token features (per token, padded & masked):
  subjet:  [log pt, mass/50, btag, dr_to_other_sj, is_subjet=1, is_sv=0]
  sv:      [log pt*, mass/5, dlensig/50, ntracks/10, is_subjet=0, is_sv=1]
(* SVs enter with summed kinematics per vertex; per-SV branches are provided
   by build_tokens() in train_partlite.py from the aggregate columns when
   per-SV info is unavailable — see notes there.)
"""
import torch
import torch.nn as nn


class Block(nn.Module):
    def __init__(self, dim, heads, mlp_ratio=2.0, drop=0.0):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, heads, dropout=drop, batch_first=True)
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, int(dim * mlp_ratio)), nn.GELU(),
            nn.Linear(int(dim * mlp_ratio), dim),
        )

    def forward(self, x, pad_mask):
        h = self.norm1(x)
        a, _ = self.attn(h, h, h, key_padding_mask=pad_mask)
        x = x + a
        return x + self.mlp(self.norm2(x))


class ClassAttention(nn.Module):
    """CLS token (built from global jet features) attends to the tokens."""

    def __init__(self, dim, heads):
        super().__init__()
        self.norm = nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, heads, batch_first=True)

    def forward(self, cls, x, pad_mask):
        h = self.norm(torch.cat([cls, x], dim=1))
        mask = torch.cat([torch.zeros_like(pad_mask[:, :1]), pad_mask], dim=1)
        out, _ = self.attn(h[:, :1], h, h, key_padding_mask=mask)
        return cls + out


class ParTLite(nn.Module):
    def __init__(self, n_token_feat, n_global, dim=64, heads=8, depth=3, drop=0.1):
        super().__init__()
        self.embed = nn.Sequential(
            nn.Linear(n_token_feat, dim), nn.GELU(), nn.LayerNorm(dim),
            nn.Linear(dim, dim), nn.GELU(), nn.LayerNorm(dim),
        )
        self.global_embed = nn.Sequential(
            nn.Linear(n_global, dim), nn.GELU(), nn.LayerNorm(dim),
            nn.Linear(dim, dim), nn.GELU(), nn.LayerNorm(dim),
        )
        self.blocks = nn.ModuleList([Block(dim, heads, drop=drop) for _ in range(depth)])
        self.cls_attn = ClassAttention(dim, heads)
        self.head = nn.Sequential(
            nn.LayerNorm(dim), nn.Linear(dim, dim), nn.GELU(),
            nn.Dropout(drop), nn.Linear(dim, 1),
        )

    def forward(self, tokens, token_mask, glob):
        """tokens: (B, T, F); token_mask: (B, T) True = real token; glob: (B, G)"""
        pad = ~token_mask
        # fully-padded rows would NaN in attention; give them one live slot
        all_pad = pad.all(dim=1)
        pad = pad.clone()
        pad[all_pad, 0] = False
        x = self.embed(tokens)
        for blk in self.blocks:
            x = blk(x, pad)
        cls = self.global_embed(glob).unsqueeze(1)
        cls = self.cls_attn(cls, x, pad)
        return self.head(cls.squeeze(1)).squeeze(-1)
