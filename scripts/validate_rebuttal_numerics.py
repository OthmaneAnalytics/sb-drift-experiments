#!/usr/bin/env python
from __future__ import annotations

import math
import os
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from sbdrift.estimator import DriftEstimator
from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import load_yaml


RATE = runpy.run_path(str(ROOT / "scripts" / "01_rate.py"))
EDGE = runpy.run_path(str(ROOT / "scripts" / "03_edge.py"))

make_bandwidth_grid = RATE["make_bandwidth_grid"]
select_lepski = RATE["select_lepski"]

edge_run_family = EDGE["run_family"]
edge_make_bandwidth_grid = EDGE["make_bandwidth_grid"]
edge_x_grid = EDGE["x_grid_from_eval_box"]


CONFIGS = [
    "gg_1d.yaml",
    "gg_2d.yaml",
    "mm_1d.yaml",
    "mm_2d.yaml",
]

checks: list[dict[str, object]] = []


def record(
    name: str,
    error: float,
    tolerance: float,
    details: str = "",
) -> None:
    passed = bool(np.isfinite(error) and error <= tolerance)
    checks.append(
        {
            "check": name,
            "passed": passed,
            "error": float(error),
            "tolerance": float(tolerance),
            "details": details,
        }
    )
    status = "PASS" if passed else "FAIL"
    print(
        f"{status:4s} | {name:48s} "
        f"| error={error:.6e} | tol={tolerance:.6e}"
    )


def direct_truth_1d(model, t: float, x: np.ndarray, xi: np.ndarray) -> np.ndarray:
    x0 = float(x[0])
    xi0 = np.asarray(xi, dtype=float).reshape(1)
    dt = float(model.u - t)
    delta = float(model.u - model.s)
    qpdf = model.conditional_pdf_fn(xi0)

    def F(y: float) -> float:
        return math.exp(
            -((y - x0) ** 2) / (2.0 * dt)
            + ((y - xi0[0]) ** 2) / (2.0 * delta)
        )

    low = float(model.low[0])
    high = float(model.high[0])

    denominator, _ = quad(
        lambda y: F(y) * qpdf(np.array([y])),
        low,
        high,
        epsabs=1e-11,
        epsrel=1e-11,
        limit=300,
    )
    numerator, _ = quad(
        lambda y: y * F(y) * qpdf(np.array([y])),
        low,
        high,
        epsabs=1e-11,
        epsrel=1e-11,
        limit=300,
    )

    return np.array([(numerator / denominator - x0) / dt])


def check_truth() -> None:
    print("\n=== Truth calculations ===")

    for config_name in ["gg_1d.yaml", "mm_1d.yaml"]:
        cfg = load_yaml(ROOT / "configs" / config_name)
        model = load_model_from_config(cfg)
        engine = TruthEngine(model)

        t = float(cfg["rate_time"])
        xi = np.asarray(cfg["xi0"], dtype=float)
        box = np.asarray(cfg["eval_box"], dtype=float)

        test_xs = [
            np.array([0.25 * box[0, 0] + 0.75 * box[0, 1]]),
            np.array([0.50 * box[0, 0] + 0.50 * box[0, 1]]),
            np.array([0.75 * box[0, 0] + 0.25 * box[0, 1]]),
        ]

        errors = []
        for x in test_xs:
            expected = direct_truth_1d(model, t, x, xi)
            observed = engine.a_star(t, x=x, xi=xi)
            errors.append(float(np.max(np.abs(observed - expected))))

        record(
            f"{config_name}: truth vs independent quadrature",
            max(errors),
            2e-9,
        )

    for config_name in ["gg_2d.yaml", "mm_2d.yaml"]:
        cfg = load_yaml(ROOT / "configs" / config_name)
        model = load_model_from_config(cfg)
        engine = TruthEngine(model)

        t = float(cfg["rate_time"])
        xi = np.asarray(cfg["xi0"], dtype=float)
        box = np.asarray(cfg["eval_box"], dtype=float)
        x = np.mean(box, axis=1)

        coarse = engine.a_star(
            t,
            x=x,
            xi=xi,
            grid_points_2d=81,
        )
        fine = engine.a_star(
            t,
            x=x,
            xi=xi,
            grid_points_2d=121,
        )

        error = float(np.max(np.abs(coarse - fine)))
        tolerance = 5e-3 + 5e-3 * float(np.max(np.abs(fine)))

        record(
            f"{config_name}: 2D truth grid refinement",
            error,
            tolerance,
            details=f"truth81={coarse.tolist()}, truth121={fine.tolist()}",
        )


def evaluation_grid(cfg: dict, dim: int) -> tuple[list[np.ndarray], np.ndarray]:
    box = np.asarray(cfg["eval_box"], dtype=float)
    axes = [
        np.linspace(box[k, 0], box[k, 1], 3 if dim == 2 else 7)
        for k in range(dim)
    ]
    mesh = np.meshgrid(*axes, indexing="ij")
    points = np.stack([part.reshape(-1) for part in mesh], axis=1)
    return axes, points


def direct_estimator(
    model,
    xs: np.ndarray,
    xu: np.ndarray,
    t: float,
    x_grid: np.ndarray,
    xi: np.ndarray,
    h: float,
) -> np.ndarray:
    estimator = DriftEstimator(model=model, xs=xs, xu=xu)
    weights = estimator.kernel_weights(xi, h)

    dt = float(model.u - t)
    delta = float(model.u - model.s)

    xg = x_grid[:, None, :]
    yy = xu[None, :, :]
    xi3 = xi.reshape(1, 1, model.dim)

    F = np.exp(
        -np.sum((yy - xg) ** 2, axis=2) / (2.0 * dt)
        + np.sum((yy - xi3) ** 2, axis=2) / (2.0 * delta)
    )

    f_hat = max(float(np.mean(weights)), 1e-12)
    g1 = np.mean(F * weights[None, :], axis=1)
    g2 = np.mean(
        F[:, :, None]
        * weights[None, :, None]
        * xu[None, :, :],
        axis=1,
    )

    D_hat = np.maximum(g1 / f_hat, 1e-12)
    N_hat = g2 / f_hat

    return (N_hat / D_hat[:, None] - x_grid) / dt


def check_estimators() -> None:
    print("\n=== Drift estimator ===")

    for index, config_name in enumerate(CONFIGS):
        cfg = load_yaml(ROOT / "configs" / config_name)
        model = load_model_from_config(cfg)

        rng = np.random.default_rng(7000 + index)
        xs, xu = model.sample(128, rng)

        _, x_grid = evaluation_grid(cfg, model.dim)
        xi = np.asarray(cfg["xi0"], dtype=float)
        t = float(cfg["rate_time"])
        h = 0.7

        estimator = DriftEstimator(model=model, xs=xs, xu=xu)
        observed = estimator.a_hat_grid(
            t=t,
            x_grid=x_grid,
            xi=xi,
            h=h,
        )
        expected = direct_estimator(
            model=model,
            xs=xs,
            xu=xu,
            t=t,
            x_grid=x_grid,
            xi=xi,
            h=h,
        )

        record(
            f"{config_name}: estimator direct reconstruction",
            float(np.max(np.abs(observed - expected))),
            2e-12,
        )


def independent_lepski(
    est_by_h: dict[float, np.ndarray],
    hs: np.ndarray,
    M: int,
    dim: int,
    kappa_pair: float,
    kappa_final: float,
) -> tuple[float, dict[float, float]]:
    pair_penalty = {
        float(h): kappa_pair
        * math.sqrt(math.log(M) / (M * float(h) ** dim))
        for h in hs
    }
    final_penalty = {
        float(h): kappa_final
        * math.sqrt(math.log(M) / (M * float(h) ** dim))
        for h in hs
    }

    scores: dict[float, float] = {}

    for h_value in hs:
        h = float(h_value)
        comparisons = []

        for hp_value in hs:
            hp = float(hp_value)
            if hp <= h:
                discrepancy = float(
                    np.max(
                        np.linalg.norm(
                            est_by_h[hp] - est_by_h[h],
                            axis=1,
                        )
                    )
                )
                comparisons.append(
                    max(0.0, discrepancy - pair_penalty[hp])
                )

        bias_proxy = max(comparisons) if comparisons else 0.0
        scores[h] = bias_proxy + final_penalty[h]

    selected = min(scores, key=scores.get)
    return float(selected), scores


def check_selector() -> None:
    print("\n=== Adaptive selector ===")

    for index, config_name in enumerate(["gg_1d.yaml", "mm_1d.yaml"]):
        cfg = load_yaml(ROOT / "configs" / config_name)
        model = load_model_from_config(cfg)

        M = 256
        rng = np.random.default_rng(8100 + index)
        xs, xu = model.sample(M, rng)

        axes, x_grid = evaluation_grid(cfg, model.dim)
        xi = np.asarray(cfg["xi0"], dtype=float)
        t = float(cfg["rate_time"])

        hs = make_bandwidth_grid(
            M=M,
            dim=model.dim,
            h0=1.2,
            q=2 ** (-0.5),
            min_h_factor=5.0,
        )

        estimator = DriftEstimator(model=model, xs=xs, xu=xu)
        est_by_h = {
            float(h): estimator.a_hat_grid(
                t=t,
                x_grid=x_grid,
                xi=xi,
                h=float(h),
            )
            for h in hs
        }

        expected, scores = independent_lepski(
            est_by_h,
            hs,
            M,
            model.dim,
            kappa_pair=2.0,
            kappa_final=2.0,
        )

        observed = select_lepski(
            est_by_h=est_by_h,
            hs=hs,
            M=M,
            dim=model.dim,
            kappa_pair=2.0,
            kappa_final=2.0,
            selector_metric="raw_max",
            trim_frac=0.05,
            penalty_form="one_sided",
            axes=axes,
        )

        record(
            f"{config_name}: Lepski score reconstruction",
            abs(float(observed) - expected),
            1e-14,
            details=f"scores={scores}",
        )


def reconstruct_clt_row(
    cfg: dict,
    row: pd.Series,
) -> dict[str, float | int]:
    model = load_model_from_config(cfg)
    clt = cfg["clt_point"]

    t = float(clt["t0"])
    x = np.asarray(clt["x0"], dtype=float)
    xi = np.asarray(clt["xi0"], dtype=float)

    M = int(row["M"])
    h = float(row["h"])
    rep_seed = int(row["seed"])

    xs, xu = model.sample(
        M,
        np.random.default_rng(rep_seed),
    )
    estimator = DriftEstimator(model=model, xs=xs, xu=xu)
    details = estimator.point_details(t=t, x=x, xi=xi, h=h)

    a_hat = float(np.asarray(details["a_hat"])[0])
    f_hat = float(details["f_hat"])
    D_hat = float(np.asarray(details["D_hat"])[0])
    F_values = np.asarray(details["F_matrix"])[0]
    weights = np.asarray(details["kernel_weights"])
    xu_values = xu.reshape(-1)

    dt = float(model.u - t)

    psi = (
        xu_values
        - float(x[0])
        - dt * a_hat
    ) * F_values

    mean_psi = float(np.mean(psi * weights) / f_hat)
    mean_psi2 = float(np.mean(psi**2 * weights) / f_hat)
    conditional_variance = max(
        mean_psi2 - mean_psi**2,
        1e-12,
    )

    sigma_hat = (
        0.6
        / (f_hat * dt**2 * D_hat**2)
        * conditional_variance
    )
    sigma_hat = max(float(sigma_hat), 1e-12)

    truth = float(
        TruthEngine(model).a_star(
            t,
            x=x,
            xi=xi,
        )[0]
    )

    Z = (
        math.sqrt(M * h)
        * (a_hat - truth)
        / math.sqrt(sigma_hat)
    )

    half_width = 1.96 * math.sqrt(sigma_hat / (M * h))
    covered = int(
        truth >= a_hat - half_width
        and truth <= a_hat + half_width
    )

    return {
        "a_hat": a_hat,
        "a_star": truth,
        "f_hat": f_hat,
        "D_hat": D_hat,
        "sigma_hat": sigma_hat,
        "Z": Z,
        "covered_95": covered,
    }


def check_clt() -> None:
    print("\n=== CLT calculation ===")

    for index, config_name in enumerate(["gg_1d.yaml", "mm_1d.yaml"]):
        cfg = load_yaml(ROOT / "configs" / config_name)
        model_id = str(cfg.get("model_id", Path(config_name).stem))

        tag = f"numerical_validation_{os.getpid()}_{index}"

        command = [
            sys.executable,
            str(ROOT / "scripts" / "02_clt.py"),
            "--config",
            f"configs/{config_name}",
            "--sample-sizes",
            "256",
            "--reps",
            "1",
            "--alpha",
            "0.24",
            "--c",
            "1.0",
            "--seed",
            str(9100 + index),
            "--out-tag",
            tag,
            "--no-qq",
            "--progress-every",
            "1",
        ]

        subprocess.run(
            command,
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )

        csv_path = (
            ROOT
            / "results"
            / "raw"
            / "clt_runs"
            / tag
            / model_id
            / "pointwise_clt.csv"
        )

        row = pd.read_csv(csv_path).iloc[0]
        expected = reconstruct_clt_row(cfg, row)

        numeric_columns = [
            "a_hat",
            "a_star",
            "f_hat",
            "D_hat",
            "sigma_hat",
            "Z",
        ]

        max_error = max(
            abs(float(row[column]) - float(expected[column]))
            for column in numeric_columns
        )

        coverage_error = abs(
            int(row["covered_95"])
            - int(expected["covered_95"])
        )

        record(
            f"{config_name}: complete CLT row reconstruction",
            max(max_error, float(coverage_error)),
            2e-11,
        )

        for path in [
            ROOT / "results" / "raw" / "clt_runs" / tag,
            ROOT / "results" / "processed" / "clt_runs" / tag,
            ROOT / "results" / "figures" / "clt_runs" / tag,
        ]:
            shutil.rmtree(path, ignore_errors=True)


def independent_edge_result(
    config_path: Path,
    M: int,
    seed: int,
) -> dict[str, np.ndarray | float]:
    cfg = load_yaml(config_path)
    model = load_model_from_config(cfg)
    engine = TruthEngine(model)

    xi = np.asarray(cfg["xi0"], dtype=float)
    t0 = float(cfg["rate_time"])
    u = float(cfg["interval"][1])

    times = np.array(
        [
            u - 0.40,
            u - 0.25,
            u - 0.15,
            u - 0.10,
            u - 0.05,
        ]
    )

    x_grid = edge_x_grid(cfg["eval_box"], 9)
    truth = {
        float(t): np.asarray(
            [engine.a_star(float(t), x=x, xi=xi) for x in x_grid]
        )
        for t in times
    }

    hs = edge_make_bandwidth_grid(
        M=M,
        dim=1,
        h0=1.2,
        q=2 ** (-0.5),
        min_h_factor=5.0,
    )

    xs, xu = model.sample(M, np.random.default_rng(seed))
    estimator = DriftEstimator(model=model, xs=xs, xu=xu)

    est_by_h = {
        float(h): estimator.a_hat_grid(
            t=t0,
            x_grid=x_grid,
            xi=xi,
            h=float(h),
        )
        for h in hs
    }

    selected, _ = independent_lepski(
        est_by_h,
        hs,
        M,
        1,
        kappa_pair=2.0,
        kappa_final=2.0,
    )

    errors = []
    for t in times:
        estimate = estimator.a_hat_grid(
            t=float(t),
            x_grid=x_grid,
            xi=xi,
            h=selected,
        )
        errors.append(
            float(
                np.max(
                    np.linalg.norm(
                        estimate - truth[float(t)],
                        axis=1,
                    )
                )
            )
        )

    errors = np.asarray(errors)

    return {
        "selected_h": selected,
        "errors": errors,
        "rescaled": (u - times) * errors,
    }


def check_edge() -> None:
    print("\n=== Terminal-edge calculation ===")

    for index, config_name in enumerate(["gg_1d.yaml", "mm_1d.yaml"]):
        config_path = ROOT / "configs" / config_name
        seed = 10100 + index
        M = 256

        observed = edge_run_family(
            cfg_path=config_path,
            M=M,
            reps=1,
            seed=seed,
            h0=1.2,
            q=2 ** (-0.5),
            min_h_factor=5.0,
            x_grid_1d=9,
            truth_grid_2d=121,
            method="lepski",
            kappa=2.0,
        )

        expected = independent_edge_result(
            config_path=config_path,
            M=M,
            seed=seed,
        )

        error = max(
            abs(float(observed["mean_h"]) - float(expected["selected_h"])),
            float(
                np.max(
                    np.abs(
                        np.asarray(observed["mean_err"])
                        - np.asarray(expected["errors"])
                    )
                )
            ),
            float(
                np.max(
                    np.abs(
                        np.asarray(observed["mean_rescaled"])
                        - np.asarray(expected["rescaled"])
                    )
                )
            ),
        )

        record(
            f"{config_name}: complete edge reconstruction",
            error,
            2e-11,
        )


def write_report() -> Path:
    output = (
        ROOT
        / "results"
        / "rebuttal"
        / "readiness"
        / "numerical_validation.md"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    passed = all(bool(check["passed"]) for check in checks)

    lines = [
        "# Numerical validation",
        "",
        f"Overall result: **{'PASS' if passed else 'FAIL'}**",
        "",
        "All calculations use fresh simulated data. Archived experiment "
        "tables are not validation inputs.",
        "",
        "| Check | Result | Maximum error | Tolerance |",
        "|---|---:|---:|---:|",
    ]

    for check in checks:
        lines.append(
            f"| {check['check']} "
            f"| {'PASS' if check['passed'] else 'FAIL'} "
            f"| {check['error']:.6e} "
            f"| {check['tolerance']:.6e} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Truth values were checked independently in 1D and by grid "
            "refinement in 2D.",
            "- Drift estimates were reconstructed directly from kernel "
            "weights, numerator, and denominator.",
            "- Adaptive bandwidth scores were independently reconstructed.",
            "- CLT output rows were regenerated from their recorded sample "
            "seeds and compared quantity by quantity.",
            "- Terminal-edge bandwidths and errors were independently "
            "reconstructed.",
            "",
        ]
    )

    output.write_text("\n".join(lines), encoding="utf-8")
    return output


def main() -> None:
    check_truth()
    check_estimators()
    check_selector()
    check_clt()
    check_edge()

    report = write_report()
    passed = all(bool(check["passed"]) for check in checks)

    print(f"\nReport: {report}")
    print(f"Overall: {'PASS' if passed else 'FAIL'}")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
