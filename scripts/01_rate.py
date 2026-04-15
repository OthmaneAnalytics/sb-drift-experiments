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

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.estimator import DriftEstimator
from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import ensure_dir, load_yaml, save_json


@dataclass
class RateConfig:
    config_path: str
    model_id: str
    name: str
    dim: int
    beta_effective: float
    sample_sizes: list[int]
    reps: int
    t0: float
    xi0: list[float]
    x_grid_per_dim: int
    truth_grid_2d: int
    h0: float
    q: float
    min_h_factor: float
    kappa: float
    seed: int


@dataclass
class SlopeSummary:
    method: str
    metric: str
    slope: float
    intercept: float
    theory_secant: float
    theory_asymptotic: float


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the finite-sample rate experiment.")
    p.add_argument("--config", type=str, required=True)
    p.add_argument("--sample-sizes", type=str, default="1000,2000,4000,8000")
    p.add_argument("--reps", type=int, default=50)
    p.add_argument("--x-grid-1d", type=int, default=200)
    p.add_argument("--x-grid-2d", type=int, default=21)
    p.add_argument("--truth-grid-2d", type=int, default=121)
    p.add_argument("--h0", type=float, default=1.2)
    p.add_argument("--q", type=float, default=2 ** (-0.5))
    p.add_argument("--min-h-factor", type=float, default=1.0)
    p.add_argument("--kappa", type=float, default=2.0)
    p.add_argument("--beta-effective", type=float, default=2.0)
    p.add_argument("--seed", type=int, default=12345)
    return p.parse_args()


def parse_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def grid_from_box(box: np.ndarray, n_per_dim: int) -> tuple[list[np.ndarray], np.ndarray]:
    axes = [np.linspace(box[k, 0], box[k, 1], n_per_dim) for k in range(box.shape[0])]
    mesh = np.meshgrid(*axes, indexing="ij")
    pts = np.stack([m.reshape(-1) for m in mesh], axis=1)
    return axes, pts


def vector_field_ise(est: np.ndarray, truth: np.ndarray, axes: list[np.ndarray]) -> float:
    diff = est - truth
    sq = np.sum(diff**2, axis=-1)
    if len(axes) == 1:
        return float(np.sqrt(np.trapezoid(sq, axes[0], axis=0)))
    shaped = sq.reshape(*[len(a) for a in axes])
    out = shaped
    for axis_vals in reversed(axes):
        out = np.trapezoid(out, axis_vals, axis=-1)
    return float(np.sqrt(out))


def sup_grid_error(est: np.ndarray, truth: np.ndarray) -> float:
    return float(np.max(np.linalg.norm(est - truth, axis=1)))


def sup_grid_distance(est_a: np.ndarray, est_b: np.ndarray) -> float:
    return float(np.max(np.linalg.norm(est_a - est_b, axis=1)))


def make_bandwidth_grid(M: int, dim: int, h0: float, q: float, min_h_factor: float) -> np.ndarray:
    min_h = min_h_factor * (M ** (-1.0 / dim))
    hs = []
    h = h0
    while h >= min_h:
        hs.append(float(h))
        h *= q
    if not hs or hs[-1] > min_h:
        hs.append(float(min_h))
    return np.array(sorted(set(round(v, 12) for v in hs)), dtype=float)


def finite_range_theory_slope(sample_sizes: list[int], beta: float, dim: int) -> float:
    M1, M2 = min(sample_sizes), max(sample_sizes)
    p = beta / (2.0 * beta + dim)
    return float(-p + p * math.log(math.log(M2) / math.log(M1)) / math.log(M2 / M1))


def asymptotic_theory_slope(beta: float, dim: int) -> float:
    p = beta / (2.0 * beta + dim)
    return float(-p)


def fit_loglog_slope(df: pd.DataFrame, metric: str, theory_secant: float, theory_asymptotic: float) -> SlopeSummary:
    x = np.log(df["M"].to_numpy(dtype=float))
    y = np.log(df[f"mean_{metric}"] .to_numpy(dtype=float))
    slope, intercept = np.polyfit(x, y, deg=1)
    return SlopeSummary(
        method=str(df["method"].iloc[0]),
        metric=metric,
        slope=float(slope),
        intercept=float(intercept),
        theory_secant=theory_secant,
        theory_asymptotic=theory_asymptotic,
    )


def make_truth_cache_path(model_id: str, t0: float, xi0: np.ndarray, n_grid: int) -> Path:
    xi_tag = "_".join(f"{v:+.3f}" for v in xi0)
    return ROOT / "results" / "processed" / "rate" / "truth_cache" / f"{model_id}_t{t0:.3f}_xi{xi_tag}_n{n_grid}.npz"


def get_truth(engine: TruthEngine, model_id: str, t0: float, xi0: np.ndarray, x_grid: np.ndarray, n_grid: int, truth_grid_2d: int) -> np.ndarray:
    cache_path = make_truth_cache_path(model_id, t0, xi0, n_grid)
    ensure_dir(cache_path.parent)
    if cache_path.exists():
        return np.asarray(np.load(cache_path)["truth"], dtype=float)
    vals = [engine.a_star(t0, x=x, xi=xi0, grid_points_2d=truth_grid_2d) for x in x_grid]
    truth = np.asarray(vals, dtype=float)
    np.savez(cache_path, truth=truth)
    return truth


def select_lepski(est_by_h: dict[float, np.ndarray], hs: np.ndarray, M: int, dim: int, kappa: float) -> float:
    penalties = {float(h): kappa * math.sqrt(math.log(M) / (M * (h ** dim))) for h in hs}
    scores: dict[float, float] = {}
    for h in hs:
        cand = []
        for hp in hs:
            if hp <= h:
                dist = sup_grid_distance(est_by_h[float(hp)], est_by_h[float(h)])
                cand.append(max(0.0, dist - penalties[float(hp)]))
        bias_proxy = max(cand) if cand else 0.0
        scores[float(h)] = bias_proxy + penalties[float(h)]
    return min(scores, key=scores.get)


def write_markdown_summary(path: Path, rcfg: RateConfig, slope_df: pd.DataFrame, mean_sup_df: pd.DataFrame, mean_ise_df: pd.DataFrame, diag_df: pd.DataFrame) -> None:
    theory_secant = slope_df["theory_secant"].iloc[0]
    theory_asymptotic = slope_df["theory_asymptotic"].iloc[0]
    lines = [
        f"# Rate experiment summary: {rcfg.name}",
        "",
        f"- model_id: `{rcfg.model_id}`",
        f"- dim: `{rcfg.dim}`",
        f"- beta_effective: `{rcfg.beta_effective}`",
        f"- sample_sizes: `{rcfg.sample_sizes}`",
        f"- reps: `{rcfg.reps}`",
        f"- t0: `{rcfg.t0}`",
        f"- xi0: `{rcfg.xi0}`",
        f"- x_grid_per_dim: `{rcfg.x_grid_per_dim}`",
        f"- truth_grid_2d: `{rcfg.truth_grid_2d}`",
        f"- h0: `{rcfg.h0}`",
        f"- q: `{rcfg.q}`",
        f"- min_h_factor: `{rcfg.min_h_factor}`",
        f"- kappa: `{rcfg.kappa}`",
        "",
        "## Theoretical benchmarks",
        "",
        f"- finite-range secant slope: `{theory_secant:.6f}`",
        f"- asymptotic slope: `{theory_asymptotic:.6f}`",
        "",
        "## Fitted slopes: sup-grid error",
        "",
    ]
    for _, row in slope_df[slope_df["metric"] == "sup_err"].iterrows():
        lines.append(f"- {row['method']}: slope = `{row['slope']:.6f}`")
    lines.extend(["", "## Mean sup-grid error by sample size", ""])
    for _, row in mean_sup_df.iterrows():
        lines.append(f"- {row['method']} at M={int(row['M'])}: mean sup-grid error = `{row['mean_sup_err']:.6e}`")
    lines.extend(["", "## Mean ISE by sample size", ""])
    for _, row in mean_ise_df.iterrows():
        lines.append(f"- {row['method']} at M={int(row['M'])}: mean ISE = `{row['mean_ise']:.6e}`")
    lines.extend(["", "## Adaptive/oracle gap (sup-grid error)", ""])
    gap_df = diag_df[diag_df["method"] == "lepski"].copy()
    for _, row in gap_df.iterrows():
        lines.append(f"- M={int(row['M'])}: mean gap = `{row['mean_gap_to_oracle']:.6f}`")
    path.write_text("\n".join(lines), encoding="utf-8")


def plot_rate_curve(mean_df: pd.DataFrame, ycol: str, ylabel: str, title: str, output_path: Path, theory_slope: float) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    for method, sub in mean_df.groupby("method"):
        ax.plot(sub["M"], sub[ycol], marker="o", label=method)
    ref = mean_df.sort_values("M")
    x = ref["M"].to_numpy(dtype=float)
    y0 = float(ref[ycol].iloc[0])
    x0 = float(ref["M"].iloc[0])
    yref = y0 * (x / x0) ** theory_slope
    ax.plot(x, yref, linestyle="--", label="finite-range theory")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = ROOT / cfg_path
    cfg = load_yaml(cfg_path)
    model = load_model_from_config(cfg)
    engine = TruthEngine(model)

    sample_sizes = parse_int_list(args.sample_sizes)
    x_grid_per_dim = args.x_grid_1d if model.dim == 1 else args.x_grid_2d
    truth_grid_2d = args.truth_grid_2d if model.dim == 2 else int(cfg.get("truth_grid_2d_fine", 121))
    eval_box = np.asarray(cfg["eval_box"], dtype=float)
    axes, x_grid = grid_from_box(eval_box, x_grid_per_dim)
    t0 = float(cfg["rate_time"])
    xi0 = np.asarray(cfg["xi0"], dtype=float)
    truth = get_truth(engine, str(cfg.get("model_id", cfg_path.stem)), t0, xi0, x_grid, x_grid_per_dim, truth_grid_2d)

    rcfg = RateConfig(
        config_path=str(cfg_path),
        model_id=str(cfg.get("model_id", cfg_path.stem)),
        name=str(cfg.get("name", cfg_path.stem)),
        dim=model.dim,
        beta_effective=float(args.beta_effective),
        sample_sizes=sample_sizes,
        reps=args.reps,
        t0=t0,
        xi0=xi0.tolist(),
        x_grid_per_dim=x_grid_per_dim,
        truth_grid_2d=truth_grid_2d,
        h0=args.h0,
        q=args.q,
        min_h_factor=args.min_h_factor,
        kappa=args.kappa,
        seed=args.seed,
    )

    raw_out = ensure_dir(ROOT / "results" / "raw" / "rate" / rcfg.model_id)
    proc_out = ensure_dir(ROOT / "results" / "processed" / "rate" / rcfg.model_id)
    fig_out = ensure_dir(ROOT / "results" / "figures")

    per_h_rows: list[dict] = []
    selected_rows: list[dict] = []
    rng_master = np.random.default_rng(args.seed)

    for M in sample_sizes:
        hs = make_bandwidth_grid(M, model.dim, args.h0, args.q, args.min_h_factor)
        seeds = rng_master.integers(0, 2**32 - 1, size=args.reps)
        for rep, rep_seed in enumerate(seeds):
            rng = np.random.default_rng(int(rep_seed))
            xs, xu = model.sample(M, rng)
            est = DriftEstimator(model=model, xs=xs, xu=xu)
            est_by_h: dict[float, np.ndarray] = {}
            sup_by_h: dict[float, float] = {}
            ise_by_h: dict[float, float] = {}
            f_by_h: dict[float, float] = {}
            minD_by_h: dict[float, float] = {}

            for h in hs:
                ah, details = est.a_hat_grid(t=t0, x_grid=x_grid, xi=xi0, h=float(h), return_details=True)
                ah = np.asarray(ah, dtype=float)
                est_by_h[float(h)] = ah
                sup_err = sup_grid_error(ah, truth)
                ise = vector_field_ise(ah, truth, axes)
                fhat = float(details["f_hat"])
                minD = float(np.min(np.asarray(details["D_hat"], dtype=float)))
                sup_by_h[float(h)] = sup_err
                ise_by_h[float(h)] = ise
                f_by_h[float(h)] = fhat
                minD_by_h[float(h)] = minD
                per_h_rows.append({
                    "model_id": rcfg.model_id,
                    "dim": rcfg.dim,
                    "M": M,
                    "rep": rep,
                    "seed": int(rep_seed),
                    "h": float(h),
                    "sup_err": float(sup_err),
                    "ise": float(ise),
                    "f_hat": fhat,
                    "min_D_hat": minD,
                })

            h_oracle = min(sup_by_h, key=sup_by_h.get)
            h_lepski = select_lepski(est_by_h, hs, M, model.dim, args.kappa)
            boundary_hs = {float(hs[0]), float(hs[-1])}
            oracle_sup = sup_by_h[h_oracle]
            for method, hsel in [("oracle", h_oracle), ("lepski", h_lepski)]:
                selected_rows.append({
                    "model_id": rcfg.model_id,
                    "dim": rcfg.dim,
                    "M": M,
                    "rep": rep,
                    "seed": int(rep_seed),
                    "method": method,
                    "selected_h": float(hsel),
                    "sup_err": float(sup_by_h[hsel]),
                    "ise": float(ise_by_h[hsel]),
                    "f_hat": float(f_by_h[hsel]),
                    "min_D_hat": float(minD_by_h[hsel]),
                    "gap_to_oracle": float(sup_by_h[hsel] / oracle_sup),
                    "is_grid_boundary": int(float(hsel) in boundary_hs),
                })
            print(f"[{rcfg.model_id}] M={M} rep={rep+1}/{args.reps} done")

    per_h_df = pd.DataFrame(per_h_rows)
    sel_df = pd.DataFrame(selected_rows)
    per_h_df.to_csv(raw_out / "per_h_metrics.csv", index=False)
    sel_df.to_csv(raw_out / "selected_metrics.csv", index=False)

    mean_sup_df = sel_df.groupby(["method", "M"], as_index=False)["sup_err"].mean().rename(columns={"sup_err": "mean_sup_err"})
    mean_ise_df = sel_df.groupby(["method", "M"], as_index=False)["ise"].mean().rename(columns={"ise": "mean_ise"})
    diag_df = sel_df.groupby(["method", "M"], as_index=False).agg(
        mean_selected_h=("selected_h", "mean"),
        mean_gap_to_oracle=("gap_to_oracle", "mean"),
        boundary_frequency=("is_grid_boundary", "mean"),
        mean_f_hat=("f_hat", "mean"),
        mean_min_D_hat=("min_D_hat", "mean"),
    )

    theory_secant = finite_range_theory_slope(sample_sizes, args.beta_effective, model.dim)
    theory_asymptotic = asymptotic_theory_slope(args.beta_effective, model.dim)
    slope_rows = []
    for method, sub in mean_sup_df.groupby("method"):
        slope_rows.append(asdict(fit_loglog_slope(sub.sort_values("M"), "sup_err", theory_secant, theory_asymptotic)))
    for method, sub in mean_ise_df.groupby("method"):
        tmp = sub.copy()
        slope_rows.append(asdict(fit_loglog_slope(tmp.sort_values("M").rename(columns={"mean_ise":"mean_ise"}), "ise", theory_secant, theory_asymptotic)))
    slope_df = pd.DataFrame(slope_rows)

    mean_sup_df.to_csv(proc_out / "mean_sup_err.csv", index=False)
    mean_ise_df.to_csv(proc_out / "mean_ise.csv", index=False)
    diag_df.to_csv(proc_out / "bandwidth_diagnostics.csv", index=False)
    slope_df.to_csv(proc_out / "slopes.csv", index=False)
    save_json(asdict(rcfg), proc_out / "run_config.json")
    write_markdown_summary(proc_out / "summary.md", rcfg, slope_df, mean_sup_df, mean_ise_df, diag_df)
    plot_rate_curve(mean_sup_df, "mean_sup_err", "mean sup-grid error", f"Rate experiment: {rcfg.name}", fig_out / f"rate_{rcfg.model_id}.pdf", theory_secant)
    plot_rate_curve(mean_ise_df, "mean_ise", "mean ISE", f"Rate experiment (ISE): {rcfg.name}", fig_out / f"rate_{rcfg.model_id}_ise.pdf", theory_secant)

    print(f"Wrote raw outputs to {raw_out}")
    print(f"Wrote processed outputs to {proc_out}")
    print(f"Wrote figures to {fig_out}")


if __name__ == "__main__":
    run(parse_args())
