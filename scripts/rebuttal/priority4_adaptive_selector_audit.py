#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from datetime import datetime
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.estimator import DriftEstimator
from sbdrift.models import load_model_from_config
from sbdrift.utils import ensure_dir, load_yaml


RATE_PATH = ROOT / "scripts" / "01_rate.py"
SPEC = importlib.util.spec_from_file_location(
    "priority4_rate_module",
    RATE_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot load {RATE_PATH}")

RATE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RATE
SPEC.loader.exec_module(RATE)


NUMERICAL_FLOOR = 1.0e-12
PRIMARY_VOLUME = 81.0


def csv_values(text: str, cast):
    values = [
        cast(value.strip())
        for value in text.split(",")
        if value.strip()
    ]

    if not values:
        raise ValueError("At least one value is required.")

    return values


def parse_kappa_pairs(text: str) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []

    for item in text.split(","):
        item = item.strip()

        if not item:
            continue

        pieces = item.split(":")

        if len(pieces) != 2:
            raise ValueError(
                "Kappa pairs must use kp:kf syntax, "
                "for example 1.5:2."
            )

        pair = (float(pieces[0]), float(pieces[1]))

        if pair not in pairs:
            pairs.append(pair)

    if not pairs:
        raise ValueError("At least one kappa pair is required.")

    return pairs


def estimator_diagnostics(
    details: dict,
    M: int,
    h: float,
    dim: int,
) -> dict[str, float | int]:
    f_hat = float(details["f_hat"])

    g1_hat = np.asarray(
        details["g1_hat"],
        dtype=float,
    ).reshape(-1)

    D_hat = np.asarray(
        details["D_hat"],
        dtype=float,
    ).reshape(-1)

    weights = np.asarray(
        details["kernel_weights"],
        dtype=float,
    ).reshape(-1)

    sum_weights = float(np.sum(weights))
    sum_squared_weights = float(np.sum(weights**2))

    kernel_ess = (
        sum_weights**2 / sum_squared_weights
        if sum_squared_weights > 0.0
        else 0.0
    )

    f_floor_hit = int(
        np.isclose(
            f_hat,
            NUMERICAL_FLOOR,
            rtol=0.0,
            atol=1.0e-15,
        )
    )

    D_floor_mask = np.isclose(
        D_hat,
        NUMERICAL_FLOOR,
        rtol=0.0,
        atol=1.0e-15,
    )

    return {
        "Mh_to_d": float(M * h**dim),
        "f_hat": f_hat,
        "f_floor_hit": f_floor_hit,
        "n_kernel_positive": int(
            np.count_nonzero(weights > 0.0)
        ),
        "kernel_weight_ess": float(kernel_ess),
        "min_g1_hat": float(np.min(g1_hat)),
        "median_g1_hat": float(np.median(g1_hat)),
        "min_D_hat": float(np.min(D_hat)),
        "q01_D_hat": float(
            np.quantile(D_hat, 0.01)
        ),
        "median_D_hat": float(np.median(D_hat)),
        "n_D_floor_grid": int(
            np.count_nonzero(D_floor_mask)
        ),
    }


def build_distance_matrix(
    estimates: dict[float, np.ndarray],
    hs: np.ndarray,
    metric: str,
    trim_frac: float,
    axes: list[np.ndarray],
) -> np.ndarray:
    n = len(hs)
    distances = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i, n):
            estimate_i = estimates[float(hs[i])]
            estimate_j = estimates[float(hs[j])]

            if metric == "raw_max":
                distance = RATE.discrepancy_raw_max(
                    estimate_i,
                    estimate_j,
                )
            elif metric == "trimmed_max":
                distance = RATE.discrepancy_trimmed_max(
                    estimate_i,
                    estimate_j,
                    trim_frac=trim_frac,
                )
            elif metric == "ise":
                distance = RATE.discrepancy_ise(
                    estimate_i,
                    estimate_j,
                    axes=axes,
                )
            else:
                raise ValueError(
                    f"Unknown selector metric: {metric}"
                )

            distances[i, j] = float(distance)
            distances[j, i] = float(distance)

    return distances


def select_from_matrix(
    distances: np.ndarray,
    hs: np.ndarray,
    M: int,
    dim: int,
    kappa_pair: float,
    kappa_final: float,
    penalty_form: str,
) -> tuple[float, np.ndarray]:
    hs = np.asarray(hs, dtype=float)

    base_penalty = np.sqrt(
        np.log(M) / (M * hs**dim)
    )

    pair_penalty = kappa_pair * base_penalty
    final_penalty = kappa_final * base_penalty

    scores = np.empty(len(hs), dtype=float)

    for j in range(len(hs)):
        candidates: list[float] = []

        for i in range(j + 1):
            if penalty_form == "one_sided":
                threshold = pair_penalty[i]
            elif penalty_form == "two_sided":
                threshold = (
                    pair_penalty[i]
                    + pair_penalty[j]
                )
            else:
                raise ValueError(
                    f"Unknown penalty form: {penalty_form}"
                )

            candidates.append(
                max(
                    0.0,
                    float(distances[i, j])
                    - float(threshold),
                )
            )

        bias_proxy = max(candidates) if candidates else 0.0

        scores[j] = (
            bias_proxy + final_penalty[j]
        )

    selected_index = int(np.argmin(scores))

    return float(hs[selected_index]), scores


def is_baseline_selector(
    metric: str,
    kappa_pair: float,
    kappa_final: float,
    penalty_form: str,
) -> bool:
    return (
        metric == "raw_max"
        and math.isclose(kappa_pair, 2.0)
        and math.isclose(kappa_final, 2.0)
        and penalty_form == "one_sided"
    )


def summarize_by_M(
    selected: pd.DataFrame,
) -> pd.DataFrame:
    group_columns = [
        "model_id",
        "M",
        "nominal_volume",
        "selector_metric",
        "kappa_pair",
        "kappa_final",
        "penalty_form",
    ]

    rows: list[dict[str, object]] = []

    for keys, group in selected.groupby(
        group_columns,
        sort=True,
    ):
        (
            model_id,
            M,
            volume,
            metric,
            kappa_pair,
            kappa_final,
            penalty_form,
        ) = keys

        rows.append(
            {
                "model_id": model_id,
                "M": int(M),
                "nominal_volume": float(volume),
                "selector_metric": metric,
                "kappa_pair": float(kappa_pair),
                "kappa_final": float(kappa_final),
                "penalty_form": penalty_form,
                "n": len(group),
                "mean_selected_h": float(
                    group["selected_h"].mean()
                ),
                "median_selected_h": float(
                    group["selected_h"].median()
                ),
                "mean_sup_err": float(
                    group["sup_err"].mean()
                ),
                "mean_ise": float(
                    group["ise"].mean()
                ),
                "mean_gap_to_oracle": float(
                    group["gap_to_oracle"].mean()
                ),
                "median_gap_to_oracle": float(
                    group["gap_to_oracle"].median()
                ),
                "q90_gap_to_oracle": float(
                    group["gap_to_oracle"].quantile(0.90)
                ),
                "max_gap_to_oracle": float(
                    group["gap_to_oracle"].max()
                ),
                "boundary_rate": float(
                    group["boundary"].mean()
                ),
                "lower_boundary_rate": float(
                    group["lower_boundary"].mean()
                ),
                "upper_boundary_rate": float(
                    group["upper_boundary"].mean()
                ),
                "same_as_floor_baseline_rate": float(
                    group["same_as_floor_baseline"].mean()
                ),
                "min_Mh_to_d": float(
                    group["Mh_to_d"].min()
                ),
                "min_f_hat": float(
                    group["f_hat"].min()
                ),
                "min_D_hat": float(
                    group["min_D_hat"].min()
                ),
                "min_kernel_ess": float(
                    group["kernel_weight_ess"].min()
                ),
                "f_floor_events": int(
                    group["f_floor_hit"].sum()
                ),
                "D_floor_grid_events": int(
                    group["n_D_floor_grid"].sum()
                ),
            }
        )

    return pd.DataFrame(rows).sort_values(
        group_columns
    ).reset_index(drop=True)


def aggregate_across_M(
    by_M: pd.DataFrame,
) -> pd.DataFrame:
    group_columns = [
        "model_id",
        "nominal_volume",
        "selector_metric",
        "kappa_pair",
        "kappa_final",
        "penalty_form",
    ]

    rows: list[dict[str, object]] = []

    for keys, group in by_M.groupby(
        group_columns,
        sort=True,
    ):
        (
            model_id,
            volume,
            metric,
            kappa_pair,
            kappa_final,
            penalty_form,
        ) = keys

        rows.append(
            {
                "model_id": model_id,
                "nominal_volume": float(volume),
                "selector_metric": metric,
                "kappa_pair": float(kappa_pair),
                "kappa_final": float(kappa_final),
                "penalty_form": penalty_form,
                "n_sample_sizes": len(group),
                "mean_sup_err_across_M": float(
                    group["mean_sup_err"].mean()
                ),
                "mean_ise_across_M": float(
                    group["mean_ise"].mean()
                ),
                "mean_gap_across_M": float(
                    group["mean_gap_to_oracle"].mean()
                ),
                "max_mean_gap_over_M": float(
                    group["mean_gap_to_oracle"].max()
                ),
                "max_q90_gap_over_M": float(
                    group["q90_gap_to_oracle"].max()
                ),
                "max_boundary_rate": float(
                    group["boundary_rate"].max()
                ),
                "mean_same_as_floor_baseline_rate": float(
                    group[
                        "same_as_floor_baseline_rate"
                    ].mean()
                ),
                "min_Mh_to_d": float(
                    group["min_Mh_to_d"].min()
                ),
                "min_f_hat": float(
                    group["min_f_hat"].min()
                ),
                "min_D_hat": float(
                    group["min_D_hat"].min()
                ),
                "min_kernel_ess": float(
                    group["min_kernel_ess"].min()
                ),
                "total_floor_events": int(
                    group[
                        [
                            "f_floor_events",
                            "D_floor_grid_events",
                        ]
                    ]
                    .to_numpy()
                    .sum()
                ),
            }
        )

    return pd.DataFrame(rows).sort_values(
        group_columns
    ).reset_index(drop=True)


def find_priority1_reference(
    model_id: str,
) -> Path | None:
    root = (
        ROOT
        / "results"
        / "raw"
        / "rate"
        / model_id
    )

    candidates = sorted(
        root.glob(
            f"{model_id}_priority1_floor81_main_*/"
            "selected_metrics.csv"
        )
    )

    return candidates[-1] if candidates else None


def compare_with_priority1(
    selected: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for model_id in sorted(
        selected["model_id"].unique()
    ):
        reference_path = find_priority1_reference(
            model_id
        )

        if reference_path is None:
            rows.append(
                {
                    "model_id": model_id,
                    "status": "MISSING_REFERENCE",
                }
            )
            continue

        reference = pd.read_csv(reference_path)

        baseline = selected[
            (selected["model_id"] == model_id)
            & np.isclose(
                selected["nominal_volume"],
                PRIMARY_VOLUME,
            )
            & (selected["selector_metric"] == "raw_max")
            & np.isclose(selected["kappa_pair"], 2.0)
            & np.isclose(selected["kappa_final"], 2.0)
            & (selected["penalty_form"] == "one_sided")
        ].copy()

        old_lepski = reference[
            reference["method"] == "lepski"
        ].copy()

        old_oracle = reference[
            reference["method"] == "oracle"
        ].copy()

        keys = ["M", "rep", "seed"]

        lepski = old_lepski.merge(
            baseline,
            on=keys,
            how="inner",
            suffixes=("_old", "_new"),
            validate="one_to_one",
        )

        oracle = old_oracle.merge(
            baseline,
            on=keys,
            how="inner",
            suffixes=("_old", "_new"),
            validate="one_to_one",
        )

        row: dict[str, object] = {
            "model_id": model_id,
            "reference_path": str(reference_path),
            "expected_rows": len(old_lepski),
            "observed_rows": len(baseline),
            "matched_lepski_rows": len(lepski),
            "matched_oracle_rows": len(oracle),
            "status": "PASS",
        }

        if (
            len(old_lepski) != len(baseline)
            or len(lepski) != len(old_lepski)
            or len(oracle) != len(old_oracle)
        ):
            row["status"] = "FAIL"

        if len(lepski):
            old_h = lepski[
                "selected_h_old"
            ].to_numpy(dtype=float)

            new_h = lepski[
                "selected_h_new"
            ].to_numpy(dtype=float)

            row["lepski_h_mismatches"] = int(
                np.count_nonzero(
                    ~np.isclose(
                        old_h,
                        new_h,
                        rtol=0.0,
                        atol=1.0e-12,
                    )
                )
            )

            row["max_abs_lepski_sup_diff"] = float(
                np.max(
                    np.abs(
                        lepski[
                            "sup_err_old"
                        ].to_numpy(dtype=float)
                        - lepski[
                            "sup_err_new"
                        ].to_numpy(dtype=float)
                    )
                )
            )

            row["max_abs_gap_diff"] = float(
                np.max(
                    np.abs(
                        lepski[
                            "gap_to_oracle_old"
                        ].to_numpy(dtype=float)
                        - lepski[
                            "gap_to_oracle_new"
                        ].to_numpy(dtype=float)
                    )
                )
            )

            row["boundary_mismatches"] = int(
                np.count_nonzero(
                    lepski[
                        "boundary_old"
                    ].to_numpy(dtype=int)
                    != lepski[
                        "boundary_new"
                    ].to_numpy(dtype=int)
                )
            )

            if (
                row["lepski_h_mismatches"]
                or row["boundary_mismatches"]
                or row["max_abs_lepski_sup_diff"] > 1.0e-8
                or row["max_abs_gap_diff"] > 1.0e-8
            ):
                row["status"] = "FAIL"

        if len(oracle):
            row["oracle_h_mismatches"] = int(
                np.count_nonzero(
                    ~np.isclose(
                        oracle[
                            "selected_h_old"
                        ].to_numpy(dtype=float),
                        oracle[
                            "oracle_h"
                        ].to_numpy(dtype=float),
                        rtol=0.0,
                        atol=1.0e-12,
                    )
                )
            )

            row["max_abs_oracle_sup_diff"] = float(
                np.max(
                    np.abs(
                        oracle[
                            "sup_err_old"
                        ].to_numpy(dtype=float)
                        - oracle[
                            "oracle_sup_err"
                        ].to_numpy(dtype=float)
                    )
                )
            )

            if (
                row["oracle_h_mismatches"]
                or row["max_abs_oracle_sup_diff"] > 1.0e-8
            ):
                row["status"] = "FAIL"

        rows.append(row)

    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Paired adaptive-selector robustness audit."
        )
    )

    parser.add_argument(
        "--configs",
        default=(
            "configs/gg_1d.yaml,"
            "configs/gg_2d.yaml,"
            "configs/mm_1d.yaml,"
            "configs/mm_2d.yaml"
        ),
    )

    parser.add_argument(
        "--sample-sizes",
        default="1000,2000,4000,8000",
    )

    parser.add_argument(
        "--reps-1d",
        type=int,
        default=50,
    )

    parser.add_argument(
        "--reps-2d",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--nominal-volumes",
        default="49,81,121",
    )

    parser.add_argument(
        "--metrics",
        default="raw_max,trimmed_max,ise",
    )

    parser.add_argument(
        "--kappa-pairs",
        default=(
            "1.5:1.5,"
            "1.5:2,"
            "2:1.5,"
            "2:2,"
            "2.5:2,"
            "2:2.5,"
            "2.5:2.5"
        ),
    )

    parser.add_argument(
        "--penalties",
        default="one_sided,two_sided",
    )

    parser.add_argument(
        "--trim-frac",
        type=float,
        default=0.05,
    )

    parser.add_argument(
        "--h0",
        type=float,
        default=1.2,
    )

    parser.add_argument(
        "--q",
        type=float,
        default=2 ** (-0.5),
    )

    parser.add_argument(
        "--x-grid-1d",
        type=int,
        default=200,
    )

    parser.add_argument(
        "--x-grid-2d",
        type=int,
        default=21,
    )

    parser.add_argument(
        "--truth-grid-2d",
        type=int,
        default=121,
    )

    parser.add_argument(
        "--truth-cache-mode",
        choices=["use", "refresh", "off"],
        default="use",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
    )

    parser.add_argument(
        "--tag",
        default="",
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--compare-priority1",
        action="store_true",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    config_paths = csv_values(args.configs, str)
    sample_sizes = csv_values(
        args.sample_sizes,
        int,
    )
    nominal_volumes = csv_values(
        args.nominal_volumes,
        float,
    )
    metrics = csv_values(args.metrics, str)
    penalties = csv_values(args.penalties, str)
    kappa_pairs = parse_kappa_pairs(
        args.kappa_pairs
    )

    allowed_metrics = {
        "raw_max",
        "trimmed_max",
        "ise",
    }

    invalid_metrics = set(metrics) - allowed_metrics

    if invalid_metrics:
        raise ValueError(
            f"Invalid selector metrics: "
            f"{sorted(invalid_metrics)}"
        )

    if "raw_max" not in metrics:
        raise ValueError(
            "The metrics list must contain raw_max "
            "for the baseline comparison."
        )

    if args.compare_priority1 and PRIMARY_VOLUME not in nominal_volumes:
        raise ValueError(
            "Priority-1 comparison requires nominal volume 81."
        )

    tag = (
        args.tag.strip()
        or datetime.now().strftime(
            "priority4_adaptive_%Y%m%d_%H%M%S"
        )
    )

    output_directory = (
        ROOT
        / "results"
        / "rebuttal"
        / "priority4"
        / tag
    )

    if output_directory.exists():
        raise FileExistsError(
            f"Output directory already exists: "
            f"{output_directory}"
        )

    ensure_dir(output_directory)

    candidate_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []

    selector_consistency_mismatches = 0

    for config_text in config_paths:
        config_path = Path(config_text)

        if not config_path.is_absolute():
            config_path = ROOT / config_path

        cfg = load_yaml(config_path)
        model = load_model_from_config(cfg)

        model_id = str(
            cfg.get("model_id", config_path.stem)
        )

        reps = (
            args.reps_1d
            if model.dim == 1
            else args.reps_2d
        )

        x_grid_per_dim = (
            args.x_grid_1d
            if model.dim == 1
            else args.x_grid_2d
        )

        t0 = float(cfg["rate_time"])
        xi0 = np.asarray(
            cfg["xi0"],
            dtype=float,
        )

        eval_box = np.asarray(
            cfg["eval_box"],
            dtype=float,
        )

        axes, x_grid = RATE.grid_from_box(
            eval_box,
            x_grid_per_dim,
        )

        engine = RATE.TruthEngine(model)

        truth = RATE.get_truth(
            engine=engine,
            cfg=cfg,
            model_id=model_id,
            t0=t0,
            xi0=xi0,
            x_grid=x_grid,
            truth_grid_2d=args.truth_grid_2d,
            cache_mode=args.truth_cache_mode,
        )

        rng_master = np.random.default_rng(
            args.seed
        )

        print()
        print("=" * 110)
        print(
            f"MODEL={model_id}, "
            f"dim={model.dim}, "
            f"reps={reps}"
        )
        print("=" * 110)

        for M in sample_sizes:
            grids_by_volume: dict[
                float,
                np.ndarray,
            ] = {}

            for volume in nominal_volumes:
                min_h_factor = volume ** (
                    1.0 / model.dim
                )

                grids_by_volume[volume] = (
                    RATE.make_bandwidth_grid(
                        M=M,
                        dim=model.dim,
                        h0=args.h0,
                        q=args.q,
                        min_h_factor=min_h_factor,
                    )
                )

            union_hs = np.array(
                sorted(
                    {
                        float(h)
                        for hs in grids_by_volume.values()
                        for h in hs
                    }
                ),
                dtype=float,
            )

            repetition_seeds = rng_master.integers(
                0,
                2**32 - 1,
                size=reps,
            )

            for rep, repetition_seed in enumerate(
                repetition_seeds
            ):
                rng = np.random.default_rng(
                    int(repetition_seed)
                )

                xs, xu = model.sample(M, rng)

                estimator = DriftEstimator(
                    model=model,
                    xs=xs,
                    xu=xu,
                )

                estimates: dict[
                    float,
                    np.ndarray,
                ] = {}

                sup_by_h: dict[float, float] = {}
                ise_by_h: dict[float, float] = {}

                diagnostics_by_h: dict[
                    float,
                    dict[str, float | int],
                ] = {}

                for h in union_hs:
                    h_value = float(h)

                    estimate_raw, details = (
                        estimator.a_hat_grid(
                            t=t0,
                            x_grid=x_grid,
                            xi=xi0,
                            h=h_value,
                            return_details=True,
                        )
                    )

                    estimate = np.asarray(
                        estimate_raw,
                        dtype=float,
                    )

                    diagnostics = estimator_diagnostics(
                        details=details,
                        M=M,
                        h=h_value,
                        dim=model.dim,
                    )

                    sup_error = RATE.sup_grid_error(
                        estimate,
                        truth,
                    )

                    ise = RATE.vector_field_ise(
                        estimate,
                        truth,
                        axes,
                    )

                    estimates[h_value] = estimate
                    sup_by_h[h_value] = float(
                        sup_error
                    )
                    ise_by_h[h_value] = float(ise)

                    diagnostics_by_h[
                        h_value
                    ] = diagnostics

                    included_volumes = [
                        volume
                        for volume, hs
                        in grids_by_volume.items()
                        if np.any(
                            np.isclose(
                                hs,
                                h_value,
                                rtol=0.0,
                                atol=1.0e-12,
                            )
                        )
                    ]

                    candidate_rows.append(
                        {
                            "model_id": model_id,
                            "dim": model.dim,
                            "M": M,
                            "rep": rep,
                            "seed": int(
                                repetition_seed
                            ),
                            "h": h_value,
                            "included_volumes": ",".join(
                                f"{volume:g}"
                                for volume
                                in included_volumes
                            ),
                            "sup_err": float(
                                sup_error
                            ),
                            "ise": float(ise),
                            **diagnostics,
                        }
                    )

                for volume in nominal_volumes:
                    hs = grids_by_volume[volume]

                    min_h = float(np.min(hs))
                    max_h = float(np.max(hs))

                    oracle_h = min(
                        (
                            float(h)
                            for h in hs
                        ),
                        key=lambda value: sup_by_h[
                            value
                        ],
                    )

                    oracle_sup = float(
                        sup_by_h[oracle_h]
                    )

                    distance_matrices = {
                        metric: build_distance_matrix(
                            estimates=estimates,
                            hs=hs,
                            metric=metric,
                            trim_frac=args.trim_frac,
                            axes=axes,
                        )
                        for metric in metrics
                    }

                    baseline_h, _ = (
                        select_from_matrix(
                            distances=
                                distance_matrices[
                                    "raw_max"
                                ],
                            hs=hs,
                            M=M,
                            dim=model.dim,
                            kappa_pair=2.0,
                            kappa_final=2.0,
                            penalty_form="one_sided",
                        )
                    )

                    reference_baseline = (
                        RATE.select_lepski(
                            est_by_h={
                                float(h):
                                    estimates[float(h)]
                                for h in hs
                            },
                            hs=hs,
                            M=M,
                            dim=model.dim,
                            kappa_pair=2.0,
                            kappa_final=2.0,
                            selector_metric="raw_max",
                            trim_frac=args.trim_frac,
                            penalty_form="one_sided",
                            axes=axes,
                        )
                    )

                    if not math.isclose(
                        baseline_h,
                        reference_baseline,
                        rel_tol=0.0,
                        abs_tol=1.0e-12,
                    ):
                        selector_consistency_mismatches += 1

                    for (
                        metric,
                        kappa_pair_tuple,
                        penalty_form,
                    ) in product(
                        metrics,
                        kappa_pairs,
                        penalties,
                    ):
                        (
                            kappa_pair,
                            kappa_final,
                        ) = kappa_pair_tuple

                        selected_h, scores = (
                            select_from_matrix(
                                distances=
                                    distance_matrices[
                                        metric
                                    ],
                                hs=hs,
                                M=M,
                                dim=model.dim,
                                kappa_pair=
                                    kappa_pair,
                                kappa_final=
                                    kappa_final,
                                penalty_form=
                                    penalty_form,
                            )
                        )

                        selected_sup = float(
                            sup_by_h[selected_h]
                        )

                        gap_to_oracle = (
                            selected_sup / oracle_sup
                        )

                        diagnostics = (
                            diagnostics_by_h[
                                selected_h
                            ]
                        )

                        matching_indices = np.where(
                            np.isclose(
                                hs,
                                selected_h,
                                rtol=0.0,
                                atol=1.0e-12,
                            )
                        )[0]

                        selected_index = int(
                            matching_indices[0]
                        )

                        selected_rows.append(
                            {
                                "model_id": model_id,
                                "dim": model.dim,
                                "M": M,
                                "rep": rep,
                                "seed": int(
                                    repetition_seed
                                ),
                                "nominal_volume":
                                    float(volume),
                                "min_h_factor": float(
                                    volume
                                    ** (
                                        1.0
                                        / model.dim
                                    )
                                ),
                                "n_bandwidths": len(hs),
                                "selector_metric":
                                    metric,
                                "trim_frac":
                                    args.trim_frac,
                                "kappa_pair":
                                    float(kappa_pair),
                                "kappa_final":
                                    float(kappa_final),
                                "penalty_form":
                                    penalty_form,
                                "selected_h":
                                    selected_h,
                                "selected_index":
                                    selected_index,
                                "selected_score":
                                    float(
                                        scores[
                                            selected_index
                                        ]
                                    ),
                                "sup_err":
                                    selected_sup,
                                "ise": float(
                                    ise_by_h[
                                        selected_h
                                    ]
                                ),
                                "oracle_h":
                                    oracle_h,
                                "oracle_sup_err":
                                    oracle_sup,
                                "oracle_ise":
                                    float(
                                        ise_by_h[
                                            oracle_h
                                        ]
                                    ),
                                "gap_to_oracle":
                                    gap_to_oracle,
                                "boundary": int(
                                    math.isclose(
                                        selected_h,
                                        min_h,
                                        abs_tol=1.0e-12,
                                    )
                                    or math.isclose(
                                        selected_h,
                                        max_h,
                                        abs_tol=1.0e-12,
                                    )
                                ),
                                "lower_boundary": int(
                                    math.isclose(
                                        selected_h,
                                        min_h,
                                        abs_tol=1.0e-12,
                                    )
                                ),
                                "upper_boundary": int(
                                    math.isclose(
                                        selected_h,
                                        max_h,
                                        abs_tol=1.0e-12,
                                    )
                                ),
                                "same_as_floor_baseline":
                                    int(
                                        math.isclose(
                                            selected_h,
                                            baseline_h,
                                            abs_tol=
                                                1.0e-12,
                                        )
                                    ),
                                "is_baseline_selector":
                                    int(
                                        is_baseline_selector(
                                            metric,
                                            kappa_pair,
                                            kappa_final,
                                            penalty_form,
                                        )
                                    ),
                                **diagnostics,
                            }
                        )

                if (
                    (rep + 1)
                    % max(
                        1,
                        args.progress_every,
                    )
                    == 0
                    or rep + 1 == reps
                ):
                    print(
                        f"[{model_id}] M={M} "
                        f"rep={rep + 1}/{reps}",
                        flush=True,
                    )

    candidates = pd.DataFrame(candidate_rows)
    selected = pd.DataFrame(selected_rows)

    by_M = summarize_by_M(selected)
    aggregate = aggregate_across_M(by_M)

    candidates.to_csv(
        output_directory
        / "paired_candidate_metrics.csv",
        index=False,
    )

    selected.to_csv(
        output_directory
        / "paired_selector_rows.csv",
        index=False,
    )

    by_M.to_csv(
        output_directory
        / "selector_summary_by_M.csv",
        index=False,
    )

    aggregate.to_csv(
        output_directory
        / "selector_aggregate_summary.csv",
        index=False,
    )

    comparison = pd.DataFrame()

    if args.compare_priority1:
        comparison = compare_with_priority1(
            selected
        )

        comparison.to_csv(
            output_directory
            / "priority1_reproduction.csv",
            index=False,
        )

    primary = aggregate[
        np.isclose(
            aggregate["nominal_volume"],
            PRIMARY_VOLUME,
        )
    ].copy()

    worst_primary = (
        primary.sort_values(
            [
                "model_id",
                "mean_gap_across_M",
            ],
            ascending=[True, False],
        )
        .groupby(
            "model_id",
            group_keys=False,
        )
        .head(10)
    )

    baseline_floor_sensitivity = aggregate[
        (aggregate["selector_metric"] == "raw_max")
        & np.isclose(
            aggregate["kappa_pair"],
            2.0,
        )
        & np.isclose(
            aggregate["kappa_final"],
            2.0,
        )
        & (
            aggregate["penalty_form"]
            == "one_sided"
        )
    ].copy()

    worst_primary.to_csv(
        output_directory
        / "floor81_worst_variants.csv",
        index=False,
    )

    baseline_floor_sensitivity.to_csv(
        output_directory
        / "baseline_floor_sensitivity.csv",
        index=False,
    )

    metadata = {
        "tag": tag,
        "configs": config_paths,
        "sample_sizes": sample_sizes,
        "reps_1d": args.reps_1d,
        "reps_2d": args.reps_2d,
        "nominal_volumes": nominal_volumes,
        "metrics": metrics,
        "kappa_pairs": kappa_pairs,
        "penalties": penalties,
        "trim_frac": args.trim_frac,
        "h0": args.h0,
        "q": args.q,
        "x_grid_1d": args.x_grid_1d,
        "x_grid_2d": args.x_grid_2d,
        "truth_grid_2d":
            args.truth_grid_2d,
        "truth_cache_mode":
            args.truth_cache_mode,
        "seed": args.seed,
        "shared_samples": True,
        "shared_candidate_estimates": True,
    }

    (
        output_directory / "run_config.json"
    ).write_text(
        json.dumps(
            metadata,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print()
    print("=" * 140)
    print("FLOOR-81 WORST SELECTOR VARIANTS")
    print("=" * 140)

    display_columns = [
        "model_id",
        "selector_metric",
        "kappa_pair",
        "kappa_final",
        "penalty_form",
        "mean_gap_across_M",
        "max_mean_gap_over_M",
        "max_q90_gap_over_M",
        "max_boundary_rate",
        "mean_same_as_floor_baseline_rate",
        "min_kernel_ess",
        "total_floor_events",
    ]

    with pd.option_context(
        "display.max_rows",
        100,
        "display.max_columns",
        30,
        "display.width",
        300,
    ):
        print(
            worst_primary[
                display_columns
            ].to_string(index=False)
        )

    print()
    print("=" * 140)
    print("BASELINE SELECTOR FLOOR SENSITIVITY")
    print("=" * 140)

    floor_columns = [
        "model_id",
        "nominal_volume",
        "mean_sup_err_across_M",
        "mean_gap_across_M",
        "max_mean_gap_over_M",
        "max_boundary_rate",
        "min_Mh_to_d",
        "min_f_hat",
        "min_D_hat",
        "min_kernel_ess",
        "total_floor_events",
    ]

    print(
        baseline_floor_sensitivity[
            floor_columns
        ].to_string(index=False)
    )

    if args.compare_priority1:
        print()
        print("=" * 140)
        print("PRIORITY-1 BASELINE REPRODUCTION")
        print("=" * 140)
        print(comparison.to_string(index=False))

    numerical_floor_events = int(
        candidates[
            [
                "f_floor_hit",
                "n_D_floor_grid",
            ]
        ]
        .to_numpy()
        .sum()
    )

    invalid_gaps = int(
        np.count_nonzero(
            selected[
                "gap_to_oracle"
            ].to_numpy(dtype=float)
            < 1.0 - 1.0e-10
        )
    )

    status_errors: list[str] = []

    if selector_consistency_mismatches:
        status_errors.append(
            f"{selector_consistency_mismatches} "
            "selector implementation mismatches"
        )

    if numerical_floor_events:
        status_errors.append(
            f"{numerical_floor_events} "
            "candidate numerical-floor events"
        )

    if invalid_gaps:
        status_errors.append(
            f"{invalid_gaps} oracle gaps below one"
        )

    if args.compare_priority1:
        if comparison.empty:
            status_errors.append(
                "Priority-1 comparison is empty"
            )
        elif np.any(
            comparison["status"] != "PASS"
        ):
            status_errors.append(
                "Priority-1 baseline reproduction failed"
            )

    print()
    print("=" * 140)
    print("FINAL STATUS")
    print("=" * 140)

    if status_errors:
        print("FAIL")

        for message in status_errors:
            print(f"- {message}")

        raise SystemExit(1)

    print("PASS")
    print(f"Outputs: {output_directory}")


if __name__ == "__main__":
    main()
