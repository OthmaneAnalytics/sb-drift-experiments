#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
from dataclasses import asdict, dataclass
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.estimator import DriftEstimator
from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import ensure_dir, load_yaml, save_json


@dataclass
class CLTConfig:
    config_path: str
    model_id: str
    name: str
    dim: int
    sample_sizes: list[int]
    reps: int
    t0: float
    x0: list[float]
    xi0: list[float]
    alpha: float
    c: float
    seed: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the pointwise CLT experiment.")
    p.add_argument("--config", type=str, required=True)
    p.add_argument("--sample-sizes", type=str, default="1000,2000,4000,8000")
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--alpha", type=float, default=0.24)
    p.add_argument("--c", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=12345)
    return p.parse_args()


def parse_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def bandwidth(M: int, alpha: float, c: float) -> float:
    return float(c * (M ** (-alpha)))


def qq_plot(z: np.ndarray, title: str, output_path: Path) -> None:
    z = np.sort(np.asarray(z, dtype=float))
    n = len(z)
    probs = (np.arange(1, n + 1) - 0.5) / n
    q_theory = stats.norm.ppf(probs)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(q_theory, z, s=12)
    lo = min(q_theory.min(), z.min())
    hi = max(q_theory.max(), z.max())
    ax.plot([lo, hi], [lo, hi], linestyle="--")
    ax.set_xlabel("Theoretical quantiles")
    ax.set_ylabel("Empirical quantiles")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def conditional_expectation(values: np.ndarray, weights: np.ndarray, fhat: float) -> float:
    return float(np.mean(values * weights) / fhat)


def run(args: argparse.Namespace) -> None:
    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = ROOT / cfg_path
    cfg = load_yaml(cfg_path)
    model = load_model_from_config(cfg)
    if model.dim != 1:
        raise ValueError("The current CLT runner is implemented for dim=1 only.")

    clt_cfg = cfg.get("clt_point")
    if clt_cfg is None:
        raise ValueError("Config must contain a clt_point block with t0, x0, and xi0.")

    t0 = float(clt_cfg["t0"])
    x0 = np.asarray(clt_cfg["x0"], dtype=float).reshape(model.dim)
    xi0 = np.asarray(clt_cfg["xi0"], dtype=float).reshape(model.dim)
    engine = TruthEngine(model)
    a_star = float(engine.a_star(t0, x=x0, xi=xi0)[0])
    delta_t = float(model.u - t0)
    sample_sizes = parse_int_list(args.sample_sizes)

    ccfg = CLTConfig(
        config_path=str(cfg_path),
        model_id=str(cfg.get("model_id", cfg_path.stem)),
        name=str(cfg.get("name", cfg_path.stem)),
        dim=model.dim,
        sample_sizes=sample_sizes,
        reps=args.reps,
        t0=t0,
        x0=x0.tolist(),
        xi0=xi0.tolist(),
        alpha=float(args.alpha),
        c=float(args.c),
        seed=args.seed,
    )

    raw_out = ensure_dir(ROOT / "results" / "raw" / "clt" / ccfg.model_id)
    proc_out = ensure_dir(ROOT / "results" / "processed" / "clt" / ccfg.model_id)
    fig_out = ensure_dir(ROOT / "results" / "figures")

    rows: list[dict] = []
    rng_master = np.random.default_rng(args.seed)

    for M in sample_sizes:
        h = bandwidth(M, args.alpha, args.c)
        seeds = rng_master.integers(0, 2**32 - 1, size=args.reps)
        for rep, rep_seed in enumerate(seeds):
            rng = np.random.default_rng(int(rep_seed))
            xs, xu = model.sample(M, rng)
            est = DriftEstimator(model=model, xs=xs, xu=xu)
            details = est.point_details(t=t0, x=x0, xi=xi0, h=h)
            ahat = float(np.asarray(details["a_hat"])[0])
            fhat = float(details["f_hat"])
            Dhat = float(np.asarray(details["D_hat"])[0])
            Fvals = np.asarray(details["F_matrix"])[0]
            w = np.asarray(details["kernel_weights"])
            xu_vec = xu.reshape(-1)
            hatpsi = (xu_vec - x0[0] - delta_t * ahat) * Fvals
            mean1 = conditional_expectation(hatpsi, w, fhat)
            mean2 = conditional_expectation(hatpsi**2, w, fhat)
            var_hat = max(mean2 - mean1**2, 1e-12)
            Rk = float(0.6)  # int k^2 for 1D Epanechnikov
            sigma_hat = (Rk / (fhat * (delta_t**2) * (Dhat**2))) * var_hat
            sigma_hat = max(float(sigma_hat), 1e-12)
            z = math.sqrt(M * h) * (ahat - a_star) / math.sqrt(sigma_hat)
            half_width = 1.96 * math.sqrt(sigma_hat / (M * h))
            covered = int((a_star >= ahat - half_width) and (a_star <= ahat + half_width))
            rows.append({
                "model_id": ccfg.model_id,
                "M": M,
                "rep": rep,
                "seed": int(rep_seed),
                "h": h,
                "a_hat": ahat,
                "a_star": a_star,
                "f_hat": fhat,
                "D_hat": Dhat,
                "sigma_hat": sigma_hat,
                "Z": float(z),
                "covered_95": covered,
            })
            print(f"[{ccfg.model_id}] M={M} rep={rep+1}/{args.reps} done")

    df = pd.DataFrame(rows)
    df.to_csv(raw_out / "pointwise_clt.csv", index=False)
    summary = df.groupby("M", as_index=False).agg(
        mean_Z=("Z", "mean"),
        var_Z=("Z", "var"),
        coverage_95=("covered_95", "mean"),
        mean_sigma_hat=("sigma_hat", "mean"),
        mean_h=("h", "mean"),
    )
    summary["coverage_95"] *= 100.0
    summary.to_csv(proc_out / "summary.csv", index=False)
    save_json(asdict(ccfg), proc_out / "run_config.json")

    for M in sample_sizes:
        z = df.loc[df["M"] == M, "Z"].to_numpy(dtype=float)
        qq_plot(z, f"QQ plot: {ccfg.name}, M={M}", fig_out / f"qq_{ccfg.model_id}_{M}.pdf")

    lines = [
        f"# CLT experiment summary: {ccfg.name}",
        "",
        f"- model_id: `{ccfg.model_id}`",
        f"- dim: `{ccfg.dim}`",
        f"- sample_sizes: `{ccfg.sample_sizes}`",
        f"- reps: `{ccfg.reps}`",
        f"- t0: `{ccfg.t0}`",
        f"- x0: `{ccfg.x0}`",
        f"- xi0: `{ccfg.xi0}`",
        f"- alpha: `{ccfg.alpha}`",
        f"- c: `{ccfg.c}`",
        "",
        "## Diagnostics",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- M={int(row['M'])}: mean Z = `{row['mean_Z']:.6f}`, var Z = `{row['var_Z']:.6f}`, coverage = `{row['coverage_95']:.2f}`%"
        )
    (proc_out / "summary.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote raw outputs to {raw_out}")
    print(f"Wrote processed outputs to {proc_out}")
    print(f"Wrote QQ plots to {fig_out}")


if __name__ == "__main__":
    run(parse_args())
