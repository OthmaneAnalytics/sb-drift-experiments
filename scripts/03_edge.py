#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import math
import yaml
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.estimator import DriftEstimator


def load_yaml(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def make_bandwidth_grid(M: int, dim: int, h0: float, q: float, min_h_factor: float) -> np.ndarray:
    min_h = min_h_factor * (M ** (-1.0 / dim))
    hs = []
    h = h0
    while h >= min_h:
        hs.append(float(h))
        h *= q
    if not hs or hs[-1] > min_h:
        hs.append(float(min_h))
    hs = sorted(set(round(v, 12) for v in hs))
    return np.array(hs, dtype=float)


def x_grid_from_eval_box(eval_box, n_per_dim: int) -> np.ndarray:
    box = np.asarray(eval_box, dtype=float)
    dim = box.shape[0]
    if dim == 1:
        return np.linspace(box[0, 0], box[0, 1], n_per_dim).reshape(-1, 1)
    raise ValueError("This script is intended for the 1D edge figure runs.")


def compute_truth_grid(engine: TruthEngine, times, x_grid, xi0, truth_grid_2d: int):
    truth = {}
    for t in times:
        vals = [engine.a_star(t, x=x, xi=xi0, grid_points_2d=truth_grid_2d) for x in x_grid]
        truth[float(t)] = np.asarray(vals, dtype=float)
    return truth


def sup_error(a_hat: np.ndarray, a_star: np.ndarray) -> float:
    return float(np.max(np.linalg.norm(a_hat - a_star, axis=1)))


def pick_h_oracle(est: DriftEstimator, hs, t0, x_grid, xi0, truth_t0):
    errs = []
    for h in hs:
        ah = est.a_hat_grid(t=t0, x_grid=x_grid, xi=xi0, h=float(h))
        errs.append((float(h), sup_error(ah, truth_t0)))
    min_err = min(e for _, e in errs)
    # tie-break toward the larger bandwidth
    cands = [h for h, e in errs if abs(e - min_err) <= 1e-14]
    return max(cands)


def pick_h_lepski(est: DriftEstimator, hs, t0, x_grid, xi0, M: int, dim: int, kappa: float):
    penalties = {float(h): kappa * math.sqrt(math.log(M) / (M * (float(h) ** dim))) for h in hs}
    est_by_h = {}
    for h in hs:
        est_by_h[float(h)] = est.a_hat_grid(t=t0, x_grid=x_grid, xi=xi0, h=float(h))

    scores = {}
    for h in hs:
        h = float(h)
        vals = []
        for hp in hs:
            hp = float(hp)
            if hp <= h:
                dist = float(np.max(np.linalg.norm(est_by_h[hp] - est_by_h[h], axis=1)))
                vals.append(max(0.0, dist - penalties[hp]))
        bias_proxy = max(vals) if vals else 0.0
        scores[h] = bias_proxy + penalties[h]

    min_score = min(scores.values())
    cands = [h for h, sc in scores.items() if abs(sc - min_score) <= 1e-14]
    return max(cands)


def run_family(cfg_path: Path, M: int, reps: int, seed: int, h0: float, q: float, min_h_factor: float,
               x_grid_1d: int, truth_grid_2d: int, method: str, kappa: float):
    cfg = load_yaml(cfg_path)
    model = load_model_from_config(cfg)
    if model.dim != 1:
        raise ValueError(f"{cfg_path} is not a 1D config.")
    engine = TruthEngine(model)

    xi0 = np.asarray(cfg["xi0"], dtype=float)
    t0 = float(cfg.get("rate_time", 0.6))
    u = float(cfg["interval"][1])
    times = np.array([u - 0.40, u - 0.25, u - 0.15, u - 0.10, u - 0.05], dtype=float)

    x_grid = x_grid_from_eval_box(cfg["eval_box"], x_grid_1d)
    truth = compute_truth_grid(engine, times, x_grid, xi0, truth_grid_2d=truth_grid_2d)
    hs = make_bandwidth_grid(M=M, dim=model.dim, h0=h0, q=q, min_h_factor=min_h_factor)

    rng = np.random.default_rng(seed)

    err_mat = []
    h_sel = []

    for _ in range(reps):
        xs, xu = model.sample(M, rng)
        est = DriftEstimator(model=model, xs=xs, xu=xu)

        if method == "oracle":
            h = pick_h_oracle(est, hs, t0, x_grid, xi0, truth[float(t0)])
        elif method == "lepski":
            h = pick_h_lepski(est, hs, t0, x_grid, xi0, M=M, dim=model.dim, kappa=kappa)
        else:
            raise ValueError("method must be 'oracle' or 'lepski'")

        h_sel.append(h)
        errs = []
        for t in times:
            ah = est.a_hat_grid(t=float(t), x_grid=x_grid, xi=xi0, h=float(h))
            errs.append(sup_error(ah, truth[float(t)]))
        err_mat.append(errs)

    err_mat = np.asarray(err_mat, dtype=float)
    mean_err = err_mat.mean(axis=0)
    std_err = err_mat.std(axis=0, ddof=1) if reps > 1 else np.zeros_like(mean_err)
    mean_rescaled = (u - times) * mean_err

    return {
        "model_name": cfg.get("name", cfg_path.stem),
        "times": times,
        "mean_err": mean_err,
        "std_err": std_err,
        "mean_rescaled": mean_rescaled,
        "mean_h": float(np.mean(h_sel)),
        "median_h": float(np.median(h_sel)),
        "method": method,
        "M": M,
        "reps": reps,
    }


def plot_curves(outdir: Path, gg, mm):
    outdir.mkdir(parents=True, exist_ok=True)

    # Raw edge error
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(gg["times"], gg["mean_err"], marker="o", label=f'{gg["model_name"]}')
    ax.plot(mm["times"], mm["mean_err"], marker="o", label=f'{mm["model_name"]}')
    ax.set_xlabel("t")
    ax.set_ylabel("Sup-grid error")
    ax.set_title("Terminal-edge sup-grid error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "edge_error.pdf")
    fig.savefig(outdir / "edge_error_sup.pdf")
    plt.close(fig)

    # Rescaled edge error
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(gg["times"], gg["mean_rescaled"], marker="o", label=f'{gg["model_name"]}')
    ax.plot(mm["times"], mm["mean_rescaled"], marker="o", label=f'{mm["model_name"]}')
    ax.set_xlabel("t")
    ax.set_ylabel("Rescaled sup-grid error")
    ax.set_title("Terminal-edge rescaled sup-grid error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "edge_rescaled.pdf")
    fig.savefig(outdir / "edge_rescaled_sup.pdf")
    plt.close(fig)


def write_summary(outdir: Path, gg, mm):
    p = outdir / "edge_summary.md"
    lines = [
        "# Edge experiment summary",
        "",
        f"- method: `{gg['method']}`",
        f"- sample size: `{gg['M']}`",
        f"- reps: `{gg['reps']}`",
        "",
        "## GG1",
        f"- mean selected h: `{gg['mean_h']:.6f}`",
        f"- median selected h: `{gg['median_h']:.6f}`",
        f"- times: `{list(map(float, gg['times']))}`",
        f"- mean error: `{list(map(float, gg['mean_err']))}`",
        f"- mean rescaled error: `{list(map(float, gg['mean_rescaled']))}`",
        "",
        "## MM1",
        f"- mean selected h: `{mm['mean_h']:.6f}`",
        f"- median selected h: `{mm['median_h']:.6f}`",
        f"- times: `{list(map(float, mm['times']))}`",
        f"- mean error: `{list(map(float, mm['mean_err']))}`",
        f"- mean rescaled error: `{list(map(float, mm['mean_rescaled']))}`",
    ]
    p.write_text("\n".join(lines), encoding="utf-8")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gg-config", default="configs/gg_1d.yaml")
    ap.add_argument("--mm-config", default="configs/mm_1d.yaml")
    ap.add_argument("--M", type=int, default=8000)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20260416)
    ap.add_argument("--method", choices=["oracle", "lepski"], default="oracle")
    ap.add_argument("--kappa", type=float, default=2.0)
    ap.add_argument("--h0", type=float, default=1.2)
    ap.add_argument("--q", type=float, default=2**(-0.5))
    ap.add_argument("--min-h-factor", type=float, default=20.0)
    ap.add_argument("--x-grid-1d", type=int, default=200)
    ap.add_argument("--truth-grid-2d", type=int, default=121)
    ap.add_argument("--outdir", default="figures")
    return ap.parse_args()


def main():
    args = parse_args()
    gg = run_family(
        cfg_path=ROOT / args.gg_config,
        M=args.M,
        reps=args.reps,
        seed=args.seed,
        h0=args.h0,
        q=args.q,
        min_h_factor=args.min_h_factor,
        x_grid_1d=args.x_grid_1d,
        truth_grid_2d=args.truth_grid_2d,
        method=args.method,
        kappa=args.kappa,
    )
    mm = run_family(
        cfg_path=ROOT / args.mm_config,
        M=args.M,
        reps=args.reps,
        seed=args.seed + 1,
        h0=args.h0,
        q=args.q,
        min_h_factor=args.min_h_factor,
        x_grid_1d=args.x_grid_1d,
        truth_grid_2d=args.truth_grid_2d,
        method=args.method,
        kappa=args.kappa,
    )

    outdir = ROOT / args.outdir
    plot_curves(outdir, gg, mm)
    write_summary(outdir, gg, mm)
    print(f"Wrote: {outdir / 'edge_error.pdf'}")
    print(f"Wrote: {outdir / 'edge_rescaled.pdf'}")
    print(f"Wrote: {outdir / 'edge_error_sup.pdf'}")
    print(f"Wrote: {outdir / 'edge_rescaled_sup.pdf'}")
    print(f"Wrote: {outdir / 'edge_summary.md'}")


if __name__ == "__main__":
    main()
