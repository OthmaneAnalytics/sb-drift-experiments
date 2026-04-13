#!/usr/bin/env python
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import sys
from typing import Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import ensure_dir, load_yaml, save_json


@dataclass
class PreflightResult:
    model_id: str
    name: str
    dim: int
    xi0: list[float]
    rate_time: float
    xs_density_at_xi0: float
    xs_density_min_on_eval_xi_grid: float
    D_star_min_on_eval_grid: float
    D_star_argmin: list[float]
    truth_convergence_abs_err: float
    truth_coarse: list[float]
    truth_fine: list[float]
    sample_mean_xs: list[float]
    sample_mean_xu: list[float]
    support_accepts_sampling: bool


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pre-flight checks for SB synthetic families and truth engine.")
    p.add_argument("--config", type=str, required=True, help="Path to a YAML config file.")
    p.add_argument("--all", action="store_true", help="Treat --config as a common config and run all listed configs.")
    p.add_argument("--grid-1d", type=int, default=121)
    p.add_argument("--grid-2d", type=int, default=9)
    p.add_argument("--sample-size", type=int, default=500)
    return p.parse_args()


def grid_from_box(box: np.ndarray, n_per_dim: int) -> np.ndarray:
    dim = box.shape[0]
    axes = [np.linspace(box[k, 0], box[k, 1], n_per_dim) for k in range(dim)]
    mesh = np.meshgrid(*axes, indexing="ij")
    pts = np.stack([m.reshape(-1) for m in mesh], axis=1)
    return pts


def run_one(config_path: Path, args: argparse.Namespace, output_dir: Path) -> PreflightResult:
    cfg = load_yaml(config_path)
    model = load_model_from_config(cfg)
    engine = TruthEngine(model)
    rng = np.random.default_rng(12345)

    xi0 = np.asarray(cfg["xi0"], dtype=float)
    t0 = float(cfg["rate_time"])
    eval_box = np.asarray(cfg["eval_box"], dtype=float)

    grid_n = args.grid_1d if model.dim == 1 else args.grid_2d
    eval_pts = grid_from_box(eval_box, grid_n)

    # Marginal density checks.
    xs_density_at_xi0 = model.xs_density(xi0)
    xs_density_min = min(model.xs_density(pt) for pt in eval_pts)

    # Denominator floor on actual experiment domain.
    denom_truth_grid = int(cfg.get("denom_truth_grid_2d", 61)) if model.dim == 2 else int(cfg.get("truth_grid_2d_coarse", 121))
    D_vals = np.array([engine.D_star(t0, x=pt, xi=xi0, grid_points_2d=denom_truth_grid) for pt in eval_pts])
    argmin = eval_pts[int(np.argmin(D_vals))]

    # Truth convergence check at one representative point.
    x_probe = np.zeros(model.dim)
    if model.dim == 1:
        truth_coarse = engine.a_star(t0, x_probe, xi0)
        truth_fine = truth_coarse.copy()
    else:
        coarse_n = int(cfg.get("truth_grid_2d_coarse", 121))
        fine_n = int(cfg.get("truth_grid_2d_fine", 241))
        truth_coarse = engine.a_star(t0, x_probe, xi0, grid_points_2d=coarse_n)
        truth_fine = engine.a_star(t0, x_probe, xi0, grid_points_2d=fine_n)
    truth_err = float(np.linalg.norm(truth_fine - truth_coarse))

    # Small sampling smoke test.
    xs, xu = model.sample(args.sample_size, rng)
    support_ok = bool(np.all(xs >= model.low) and np.all(xs <= model.high) and np.all(xu >= model.low) and np.all(xu <= model.high))

    result = PreflightResult(
        model_id=str(cfg.get("model_id", config_path.stem)),
        name=str(cfg.get("name", config_path.stem)),
        dim=model.dim,
        xi0=xi0.tolist(),
        rate_time=t0,
        xs_density_at_xi0=float(xs_density_at_xi0),
        xs_density_min_on_eval_xi_grid=float(xs_density_min),
        D_star_min_on_eval_grid=float(np.min(D_vals)),
        D_star_argmin=argmin.tolist(),
        truth_convergence_abs_err=truth_err,
        truth_coarse=np.asarray(truth_coarse, dtype=float).tolist(),
        truth_fine=np.asarray(truth_fine, dtype=float).tolist(),
        sample_mean_xs=np.mean(xs, axis=0).tolist(),
        sample_mean_xu=np.mean(xu, axis=0).tolist(),
        support_accepts_sampling=support_ok,
    )

    save_json(asdict(result), output_dir / f"{result.model_id}.json")
    write_markdown(result, output_dir / f"{result.model_id}.md")
    return result


def write_markdown(result: PreflightResult, path: Path) -> None:
    lines = [
        f"# Pre-flight report: {result.name}",
        "",
        f"- model_id: `{result.model_id}`",
        f"- dim: `{result.dim}`",
        f"- xi0: `{result.xi0}`",
        f"- rate_time: `{result.rate_time}`",
        f"- xs_density_at_xi0: `{result.xs_density_at_xi0:.6e}`",
        f"- xs_density_min_on_eval_xi_grid: `{result.xs_density_min_on_eval_xi_grid:.6e}`",
        f"- D_star_min_on_eval_grid: `{result.D_star_min_on_eval_grid:.6e}`",
        f"- D_star_argmin: `{result.D_star_argmin}`",
        f"- truth_convergence_abs_err: `{result.truth_convergence_abs_err:.6e}`",
        f"- truth_coarse: `{result.truth_coarse}`",
        f"- truth_fine: `{result.truth_fine}`",
        f"- sample_mean_xs: `{result.sample_mean_xs}`",
        f"- sample_mean_xu: `{result.sample_mean_xu}`",
        f"- support_accepts_sampling: `{result.support_accepts_sampling}`",
        "",
        "## Interpretation",
        "",
        "Healthy runs should have: positive and non-tiny marginal density at the chosen conditioning point,",
        "a denominator floor on the actual evaluation box that is comfortably away from zero,",
        "and a small coarse-vs-fine truth error at the probe point.",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def load_config_paths(common_cfg_path: Path) -> Iterable[Path]:
    cfg = load_yaml(common_cfg_path)
    base = common_cfg_path.parent.parent if common_cfg_path.name == "common.yaml" else common_cfg_path.parent
    for rel in cfg["project"]["all_configs"]:
        yield ROOT / rel


def main() -> None:
    args = parse_args()
    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = ROOT / cfg_path

    if args.all:
        common = load_yaml(cfg_path)
        output_dir = ensure_dir(ROOT / common["project"]["output_dir"])
        results = [run_one(path, args, output_dir) for path in load_config_paths(cfg_path)]
        save_json([asdict(r) for r in results], output_dir / "summary.json")
        summary_md = ["# Pre-flight summary", ""]
        for r in results:
            summary_md.extend([
                f"## {r.name}",
                f"- xs_density_at_xi0: `{r.xs_density_at_xi0:.6e}`",
                f"- xs_density_min_on_eval_xi_grid: `{r.xs_density_min_on_eval_xi_grid:.6e}`",
                f"- D_star_min_on_eval_grid: `{r.D_star_min_on_eval_grid:.6e}`",
                f"- truth_convergence_abs_err: `{r.truth_convergence_abs_err:.6e}`",
                f"- support_accepts_sampling: `{r.support_accepts_sampling}`",
                "",
            ])
        (output_dir / "summary.md").write_text("\n".join(summary_md), encoding="utf-8")
        print(f"Wrote summary to {output_dir / 'summary.md'}")
    else:
        output_dir = ensure_dir(ROOT / "results/processed/preflight")
        result = run_one(cfg_path, args, output_dir)
        print(f"Wrote report to {output_dir / (result.model_id + '.md')}")


if __name__ == "__main__":
    main()
