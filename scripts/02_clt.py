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
    out_tag: str
    qq_ms: list[int]
    no_qq: bool
    progress_every: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the pointwise CLT experiment.")
    p.add_argument("--config", type=str, required=True)
    p.add_argument("--sample-sizes", type=str, default="1000,2000,4000,8000")
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--alpha", type=float, default=0.24)
    p.add_argument("--c", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--out-tag", type=str, default="")
    p.add_argument("--qq-ms", type=str, default="")
    p.add_argument("--no-qq", action="store_true")
    p.add_argument("--progress-every", type=int, default=10)
    return p.parse_args()


def parse_int_list(s: str) -> list[int]:
    if not s.strip():
        return []
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


def get_output_dirs(model_id: str, out_tag: str) -> tuple[Path, Path, Path]:
    if out_tag:
        raw_out = ensure_dir(ROOT / "results" / "raw" / "clt_runs" / out_tag / model_id)
        proc_out = ensure_dir(ROOT / "results" / "processed" / "clt_runs" / out_tag / model_id)
        fig_out = ensure_dir(ROOT / "results" / "figures" / "clt_runs" / out_tag)
    else:
        raw_out = ensure_dir(ROOT / "results" / "raw" / "clt" / model_id)
        proc_out = ensure_dir(ROOT / "results" / "processed" / "clt" / model_id)
        fig_out = ensure_dir(ROOT / "results" / "figures")
    return raw_out, proc_out, fig_out


def shapiro_safe(z: np.ndarray) -> tuple[float, float]:
    z = np.asarray(z, dtype=float)
    if len(z) < 3:
        return float("nan"), float("nan")
    if len(z) > 5000:
        z = z[:5000]
    res = stats.shapiro(z)
    return float(res.statistic), float(res.pvalue)


def anderson_safe(z: np.ndarray) -> tuple[float, float, int]:
    z = np.asarray(z, dtype=float)
    if len(z) < 8:
        return float("nan"), float("nan"), -1
    res = stats.anderson(z, dist="norm")
    sig = np.asarray(res.significance_level, dtype=float)
    crit = np.asarray(res.critical_values, dtype=float)
    idx = int(np.argmin(np.abs(sig - 5.0)))
    stat = float(res.statistic)
    crit_5 = float(crit[idx])
    reject_5 = int(stat > crit_5)
    return stat, crit_5, reject_5


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for M, g in df.groupby("M"):
        z = g["Z"].to_numpy(dtype=float)
        shapiro_W, shapiro_p = shapiro_safe(z)
        ad_stat, ad_crit_5, ad_reject_5 = anderson_safe(z)
        rows.append(
            {
                "M": int(M),
                "mean_Z": float(np.mean(z)),
                "var_Z": float(np.var(z, ddof=1)) if len(z) > 1 else float("nan"),
                "coverage_95": float(100.0 * g["covered_95"].mean()),
                "mean_sigma_hat": float(g["sigma_hat"].mean()),
                "mean_h": float(g["h"].mean()),
                "shapiro_W": shapiro_W,
                "shapiro_p": shapiro_p,
                "anderson_stat": ad_stat,
                "anderson_crit_5pct": ad_crit_5,
                "anderson_reject_5pct": ad_reject_5,
            }
        )
    return pd.DataFrame(rows).sort_values("M").reset_index(drop=True)


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
    sample_sizes = parse_int_list(args.sample_sizes)
    qq_ms = parse_int_list(args.qq_ms)

    engine = TruthEngine(model)
    a_star = float(engine.a_star(t0, x=x0, xi=xi0)[0])
    delta_t = float(model.u - t0)

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
        out_tag=str(args.out_tag),
        qq_ms=qq_ms,
        no_qq=bool(args.no_qq),
        progress_every=max(1, int(args.progress_every)),
    )

    raw_out, proc_out, fig_out = get_output_dirs(ccfg.model_id, ccfg.out_tag)

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

            Rk = float(0.6)  # ∫ K^2 for 1D Epanechnikov
            sigma_hat = (Rk / (fhat * (delta_t**2) * (Dhat**2))) * var_hat
            sigma_hat = max(float(sigma_hat), 1e-12)

            z = math.sqrt(M * h) * (ahat - a_star) / math.sqrt(sigma_hat)
            half_width = 1.96 * math.sqrt(sigma_hat / (M * h))
            covered = int((a_star >= ahat - half_width) and (a_star <= ahat + half_width))

            rows.append(
                {
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
                }
            )

            if (rep + 1) % ccfg.progress_every == 0 or (rep + 1) == args.reps:
                print(f"[{ccfg.model_id}] M={M} rep={rep+1}/{args.reps} done")

    df = pd.DataFrame(rows)
    df.to_csv(raw_out / "pointwise_clt.csv", index=False)

    summary = summarize(df)
    summary.to_csv(proc_out / "summary.csv", index=False)
    save_json(asdict(ccfg), proc_out / "run_config.json")

    qq_set = set(qq_ms)
    if not args.no_qq:
        for M in sample_sizes:
            if qq_set and M not in qq_set:
                continue
            z = df.loc[df["M"] == M, "Z"].to_numpy(dtype=float)
            suffix = f"{ccfg.out_tag}_" if ccfg.out_tag else ""
            qq_plot(z, f"QQ plot: {ccfg.name}, M={M}", fig_out / f"qq_{suffix}{ccfg.model_id}_{M}.pdf")

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
        f"- out_tag: `{ccfg.out_tag}`",
        f"- qq_ms: `{ccfg.qq_ms}`",
        "",
        "## Diagnostics",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- M={int(row['M'])}: mean Z = `{row['mean_Z']:.6f}`, "
            f"var Z = `{row['var_Z']:.6f}`, "
            f"coverage = `{row['coverage_95']:.2f}`%, "
            f"Shapiro p = `{row['shapiro_p']:.6f}`, "
            f"Anderson stat = `{row['anderson_stat']:.6f}`, "
            f"5% crit = `{row['anderson_crit_5pct']:.6f}`, "
            f"reject@5% = `{int(row['anderson_reject_5pct']) if row['anderson_reject_5pct'] >= 0 else -1}`"
        )
    (proc_out / "summary.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote raw outputs to {raw_out}")
    print(f"Wrote processed outputs to {proc_out}")
    print(f"Wrote QQ plots to {fig_out}")


if __name__ == "__main__":
    run(parse_args())
