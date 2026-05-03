#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ----------------------------
# Global plotting style
# ----------------------------
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 17,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "figure.titlesize": 18,
    "lines.linewidth": 2.2,
    "lines.markersize": 6.5,
})

# ----------------------------
# Helpers
# ----------------------------
def newest(paths):
    paths = [Path(p) for p in paths]
    if not paths:
        return None
    return max(paths, key=lambda p: p.stat().st_mtime)

def choose_col(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Could not find any of {candidates} in columns {list(df.columns)}")

# ----------------------------
# 1) Adaptivity figure
# ----------------------------
def build_adaptivity_gap():
    per_m_files = sorted(glob.glob(str(ROOT / "results" / "processed" / "adapt_final" / "*adapt_dense_floor81*_perM.csv")))
    if not per_m_files:
        raise FileNotFoundError("Could not find adaptivity per-M file under results/processed/adapt_final/")
    f = newest(per_m_files)
    df = pd.read_csv(f)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.9), sharex=True, sharey=True)

    panels = [
        (axes[0], ["GG1", "GG2"], "GG testbeds"),
        (axes[1], ["MM1", "MM2"], "MM testbeds"),
    ]

    for ax, labels, title in panels:
        handles = []
        names = []

        for label in labels:
            g = df[df["label"] == label].sort_values("M")
            x = g["M"].to_numpy(dtype=float)
            y = g["mean_gap"].to_numpy(dtype=float)
            q25 = g["q25_gap"].to_numpy(dtype=float)
            q75 = g["q75_gap"].to_numpy(dtype=float)

            lx = np.log10(x)
            lx_dense = np.linspace(lx.min(), lx.max(), 300)
            x_dense = 10 ** lx_dense

            y_dense = PchipInterpolator(lx, y)(lx_dense)
            q25_dense = PchipInterpolator(lx, q25)(lx_dense)
            q75_dense = PchipInterpolator(lx, q75)(lx_dense)

            line, = ax.plot(x_dense, y_dense, label=label)
            ax.fill_between(x_dense, q25_dense, q75_dense, alpha=0.18, color=line.get_color())
            ax.plot(x, y, marker="o", linestyle="none", color=line.get_color())

            handles.append(line)
            names.append(label)

        ref = ax.axhline(4.0, linestyle="--", linewidth=1.6, color="black")
        handles.append(ref)
        names.append("GL ref. 4")

        ax.set_xscale("log")
        ax.set_title(title, pad=10)
        ax.set_xlabel("Sample size $M$")
        ax.set_ylim(1.0, 4.2)
        ax.set_yticks([1, 2, 3, 4])
        ax.grid(True, alpha=0.25)
        ax.tick_params(axis="both", which="major", pad=4)

        ax.legend(
            handles, names,
            loc="upper left",
            frameon=True,
            facecolor="white",
            framealpha=0.95,
            edgecolor="0.8",
        )

    axes[0].set_ylabel(r"Empirical oracle ratio $\bar{\Gamma}_\nu(M)$")

    plt.tight_layout()
    out = FIG_DIR / "adaptivity_gap.pdf"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out}")

# ----------------------------
# 2) QQ plot figure
# ----------------------------
def find_clt_tag():
    env_tag = os.environ.get("CLT_TAG", "").strip()
    if env_tag:
        gg = ROOT / "results" / "raw" / "clt_runs" / env_tag / "gg_1d" / "pointwise_clt.csv"
        mm = ROOT / "results" / "raw" / "clt_runs" / env_tag / "mm_1d" / "pointwise_clt.csv"
        if gg.exists() and mm.exists():
            return env_tag
        raise FileNotFoundError(f"CLT_TAG={env_tag} was set, but gg_1d/mm_1d pointwise_clt.csv were not both found.")

    candidates = []
    for p in (ROOT / "results" / "raw" / "clt_runs").glob("*"):
        if not p.is_dir():
            continue
        name = p.name.lower()
        if any(bad in name for bad in ["timing_probe", "screen", "smoke", "alpha_screen"]):
            continue
        gg = p / "gg_1d" / "pointwise_clt.csv"
        mm = p / "mm_1d" / "pointwise_clt.csv"
        if gg.exists() and mm.exists():
            candidates.append(p)

    if not candidates:
        raise FileNotFoundError(
            "Could not auto-detect final CLT run. "
            "Set it explicitly, e.g. export CLT_TAG=<your_final_clt_tag> and rerun."
        )
    return max(candidates, key=lambda p: p.stat().st_mtime).name

def build_qqplot():
    tag = find_clt_tag()
    print(f"Using CLT tag: {tag}")

    gg_path = ROOT / "results" / "raw" / "clt_runs" / tag / "gg_1d" / "pointwise_clt.csv"
    mm_path = ROOT / "results" / "raw" / "clt_runs" / tag / "mm_1d" / "pointwise_clt.csv"

    gg = pd.read_csv(gg_path)
    mm = pd.read_csv(mm_path)

    zcol_gg = choose_col(gg, ["z", "z_stat", "zscore", "z_score", "Z"])
    zcol_mm = choose_col(mm, ["z", "z_stat", "zscore", "z_score", "Z"])

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.2))
    ms = [4000, 8000]

    cases = [
        ("GG1", gg, zcol_gg, 0),
        ("MM1", mm, zcol_mm, 1),
    ]

    for row_label, df, zcol, i in cases:
        for j, M in enumerate(ms):
            ax = axes[i, j]
            z = df.loc[df["M"] == M, zcol].dropna().to_numpy(dtype=float)

            (theo_q, samp_q), (slope, intercept, _) = stats.probplot(z, dist="norm")
            ax.scatter(theo_q, samp_q, s=20, alpha=0.85)
            ax.plot(theo_q, slope * theo_q + intercept, linewidth=2.0)

            ax.set_title(f"{row_label}, $M={M}$", pad=8)
            ax.grid(True, alpha=0.25)

            if i == 1:
                ax.set_xlabel("Theoretical quantiles")
            if j == 0:
                ax.set_ylabel("Sample quantiles")

    plt.tight_layout()
    out = FIG_DIR / "qqplot.pdf"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out}")

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    build_adaptivity_gap()
    build_qqplot()
