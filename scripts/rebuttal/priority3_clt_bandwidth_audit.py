#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.estimator import DriftEstimator
from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import ensure_dir, load_yaml


RK_EPANECHNIKOV_1D = 0.6
FLOOR = 1.0e-12

FINAL_ALPHA = {
    "gg_1d": 0.22,
    "mm_1d": 0.28,
}

FINAL_RUN_PATH = (
    ROOT
    / "results"
    / "raw"
    / "clt_runs"
    / "clt_final_20260417_101207"
)


def csv_values(text: str, cast):
    values = [
        cast(value.strip())
        for value in text.split(",")
        if value.strip()
    ]

    if not values:
        raise ValueError("At least one value is required.")

    return values


def bandwidth(M: int, alpha: float, constant: float) -> float:
    return float(constant * M ** (-alpha))


def conditional_expectation(
    values: np.ndarray,
    weights: np.ndarray,
    f_hat: float,
) -> float:
    return float(np.mean(values * weights) / f_hat)


def safe_skewness(values: np.ndarray) -> float:
    if len(values) < 3:
        return float("nan")

    return float(stats.skew(values, bias=False))


def safe_excess_kurtosis(values: np.ndarray) -> float:
    if len(values) < 4:
        return float("nan")

    return float(
        stats.kurtosis(
            values,
            fisher=True,
            bias=False,
        )
    )


def summarize_rows(rows: pd.DataFrame) -> pd.DataFrame:
    summaries: list[dict[str, object]] = []

    group_columns = [
        "model_id",
        "M",
        "alpha",
        "c",
    ]

    for keys, group in rows.groupby(
        group_columns,
        sort=True,
    ):
        model_id, M, alpha, constant = keys

        z = group["Z"].to_numpy(dtype=float)
        error = group["error"].to_numpy(dtype=float)
        a_hat = group["a_hat"].to_numpy(dtype=float)
        se_hat = group["se_hat"].to_numpy(dtype=float)
        covered = group["covered_95"].to_numpy(dtype=float)

        n = len(group)
        mean_z = float(np.mean(z))
        sd_z = (
            float(np.std(z, ddof=1))
            if n > 1
            else float("nan")
        )

        coverage = float(np.mean(covered))
        mean_error = float(np.mean(error))
        mean_se = float(np.mean(se_hat))

        empirical_sd = (
            float(np.std(a_hat, ddof=1))
            if n > 1
            else float("nan")
        )

        if n >= 2:
            jb = stats.jarque_bera(z)
            jb_stat = float(jb.statistic)
            jb_p = float(jb.pvalue)
        else:
            jb_stat = float("nan")
            jb_p = float("nan")

        summaries.append(
            {
                "model_id": str(model_id),
                "M": int(M),
                "alpha": float(alpha),
                "c": float(constant),
                "n": n,
                "h": float(group["h"].iloc[0]),
                "Mh": float(group["Mh"].iloc[0]),
                "mean_Z": mean_z,
                "mcse_mean_Z": (
                    sd_z / math.sqrt(n)
                    if n > 1
                    else float("nan")
                ),
                "sd_Z": sd_z,
                "var_Z": sd_z**2 if np.isfinite(sd_z) else np.nan,
                "q025_Z": float(np.quantile(z, 0.025)),
                "median_Z": float(np.quantile(z, 0.50)),
                "q975_Z": float(np.quantile(z, 0.975)),
                "coverage_95": coverage,
                "coverage_95_percent": 100.0 * coverage,
                "mcse_coverage": math.sqrt(
                    coverage * (1.0 - coverage) / n
                ),
                "skewness": safe_skewness(z),
                "excess_kurtosis":
                    safe_excess_kurtosis(z),
                "jarque_bera_stat": jb_stat,
                "jarque_bera_p": jb_p,
                "mean_error": mean_error,
                "rmse": float(
                    np.sqrt(np.mean(error**2))
                ),
                "empirical_sd_a_hat": empirical_sd,
                "mean_se_hat": mean_se,
                "bias_over_mean_se": (
                    mean_error / mean_se
                    if mean_se > 0.0
                    else np.nan
                ),
                "empirical_sd_over_mean_se": (
                    empirical_sd / mean_se
                    if mean_se > 0.0
                    else np.nan
                ),
                "mean_sigma_hat": float(
                    group["sigma_hat"].mean()
                ),
                "min_f_hat": float(
                    group["f_hat"].min()
                ),
                "min_D_hat": float(
                    group["D_hat"].min()
                ),
                "min_kernel_ess": float(
                    group["kernel_weight_ess"].min()
                ),
                "min_kernel_positive": int(
                    group["n_kernel_positive"].min()
                ),
                "max_abs_conditional_mean": float(
                    group["conditional_mean_hat"]
                    .abs()
                    .max()
                ),
                "f_floor_events": int(
                    group["f_floor_hit"].sum()
                ),
                "D_floor_events": int(
                    group["D_floor_hit"].sum()
                ),
                "variance_floor_events": int(
                    group["variance_floor_hit"].sum()
                ),
                "sigma_floor_events": int(
                    group["sigma_floor_hit"].sum()
                ),
            }
        )

    return pd.DataFrame(summaries).sort_values(
        group_columns
    ).reset_index(drop=True)


def aggregate_summary(
    summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for keys, group in summary.groupby(
        ["model_id", "alpha", "c"],
        sort=True,
    ):
        model_id, alpha, constant = keys

        rows.append(
            {
                "model_id": model_id,
                "alpha": float(alpha),
                "c": float(constant),
                "n_sample_sizes": len(group),
                "min_Mh": float(group["Mh"].min()),
                "max_abs_mean_Z": float(
                    group["mean_Z"].abs().max()
                ),
                "mean_abs_mean_Z": float(
                    group["mean_Z"].abs().mean()
                ),
                "max_abs_sd_Z_minus_1": float(
                    (group["sd_Z"] - 1.0).abs().max()
                ),
                "mean_sd_Z": float(
                    group["sd_Z"].mean()
                ),
                "min_coverage_95_percent": float(
                    group["coverage_95_percent"].min()
                ),
                "max_coverage_95_percent": float(
                    group["coverage_95_percent"].max()
                ),
                "mean_coverage_95_percent": float(
                    group["coverage_95_percent"].mean()
                ),
                "max_abs_bias_over_mean_se": float(
                    group["bias_over_mean_se"].abs().max()
                ),
                "max_abs_empirical_sd_ratio_minus_1":
                    float(
                        (
                            group[
                                "empirical_sd_over_mean_se"
                            ]
                            - 1.0
                        )
                        .abs()
                        .max()
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
                            "D_floor_events",
                            "variance_floor_events",
                            "sigma_floor_events",
                        ]
                    ]
                    .to_numpy()
                    .sum()
                ),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["model_id", "alpha", "c"]
    ).reset_index(drop=True)


def compare_with_final(
    rows: pd.DataFrame,
) -> pd.DataFrame:
    comparison_rows: list[dict[str, object]] = []

    numeric_columns = [
        "h",
        "a_hat",
        "a_star",
        "f_hat",
        "D_hat",
        "sigma_hat",
        "Z",
    ]

    for model_id, alpha in FINAL_ALPHA.items():
        final_path = (
            FINAL_RUN_PATH
            / model_id
            / "pointwise_clt.csv"
        )

        if not final_path.exists():
            comparison_rows.append(
                {
                    "model_id": model_id,
                    "status": "missing_final_file",
                }
            )
            continue

        expected = pd.read_csv(final_path)

        observed = rows[
            (rows["model_id"] == model_id)
            & np.isclose(
                rows["alpha"],
                alpha,
                rtol=0.0,
                atol=1e-15,
            )
            & np.isclose(
                rows["c"],
                1.0,
                rtol=0.0,
                atol=1e-15,
            )
        ].copy()

        merge_columns = ["M", "rep", "seed"]

        merged = expected.merge(
            observed,
            on=merge_columns,
            how="inner",
            suffixes=("_old", "_new"),
            validate="one_to_one",
        )

        row: dict[str, object] = {
            "model_id": model_id,
            "alpha": alpha,
            "expected_rows": len(expected),
            "observed_rows": len(observed),
            "matched_rows": len(merged),
            "status": "PASS",
        }

        if (
            len(expected) != len(observed)
            or len(merged) != len(expected)
        ):
            row["status"] = "FAIL"

        for column in numeric_columns:
            difference = np.abs(
                merged[f"{column}_new"].to_numpy(dtype=float)
                - merged[f"{column}_old"].to_numpy(dtype=float)
            )

            maximum = (
                float(np.max(difference))
                if len(difference)
                else np.nan
            )

            row[f"max_abs_diff_{column}"] = maximum

            if np.isfinite(maximum) and maximum > 1e-12:
                row["status"] = "FAIL"

        coverage_difference = np.abs(
            merged["covered_95_new"].to_numpy(dtype=int)
            - merged["covered_95_old"].to_numpy(dtype=int)
        )

        row["coverage_mismatches"] = int(
            np.count_nonzero(coverage_difference)
        )

        if row["coverage_mismatches"]:
            row["status"] = "FAIL"

        comparison_rows.append(row)

    return pd.DataFrame(comparison_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Paired CLT bandwidth sensitivity audit. "
            "Each simulated sample is reused for every "
            "(alpha, c) combination."
        )
    )

    parser.add_argument(
        "--configs",
        default=(
            "configs/gg_1d.yaml,"
            "configs/mm_1d.yaml"
        ),
    )
    parser.add_argument(
        "--sample-sizes",
        default="1000,2000,4000,8000",
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=300,
    )
    parser.add_argument(
        "--alphas",
        default="0.22,0.24,0.26,0.28",
    )
    parser.add_argument(
        "--constants",
        default="0.75,1.0,1.25",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260417,
    )
    parser.add_argument(
        "--tag",
        default="",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
    )
    parser.add_argument(
        "--compare-final",
        action="store_true",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    config_paths = csv_values(args.configs, str)
    sample_sizes = csv_values(args.sample_sizes, int)
    alphas = csv_values(args.alphas, float)
    constants = csv_values(args.constants, float)

    combinations = list(product(alphas, constants))

    tag = (
        args.tag.strip()
        or datetime.now().strftime(
            "priority3_clt_%Y%m%d_%H%M%S"
        )
    )

    output_directory = ensure_dir(
        ROOT
        / "results"
        / "rebuttal"
        / "priority3"
        / tag
    )

    rows: list[dict[str, object]] = []

    for config_text in config_paths:
        config_path = Path(config_text)

        if not config_path.is_absolute():
            config_path = ROOT / config_path

        cfg = load_yaml(config_path)
        model = load_model_from_config(cfg)

        if model.dim != 1:
            raise ValueError(
                "Priority 3 currently supports 1D models only."
            )

        model_id = str(
            cfg.get("model_id", config_path.stem)
        )

        clt_point = cfg.get("clt_point")

        if clt_point is None:
            raise ValueError(
                f"{config_path} has no clt_point block."
            )

        t0 = float(clt_point["t0"])
        x0 = np.asarray(
            clt_point["x0"],
            dtype=float,
        ).reshape(1)
        xi0 = np.asarray(
            clt_point["xi0"],
            dtype=float,
        ).reshape(1)

        engine = TruthEngine(model)
        a_star = float(
            engine.a_star(
                t0,
                x=x0,
                xi=xi0,
            )[0]
        )
        delta_t = float(model.u - t0)

        rng_master = np.random.default_rng(args.seed)

        print()
        print("=" * 100)
        print(
            f"MODEL={model_id}, "
            f"combinations={len(combinations)}, "
            f"reps={args.reps}"
        )
        print("=" * 100)

        for M in sample_sizes:
            repetition_seeds = rng_master.integers(
                0,
                2**32 - 1,
                size=args.reps,
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

                for alpha, constant in combinations:
                    h = bandwidth(
                        M,
                        alpha,
                        constant,
                    )

                    details = estimator.point_details(
                        t=t0,
                        x=x0,
                        xi=xi0,
                        h=h,
                    )

                    a_hat = float(
                        np.asarray(
                            details["a_hat"],
                            dtype=float,
                        ).reshape(-1)[0]
                    )
                    f_hat = float(details["f_hat"])
                    D_hat = float(
                        np.asarray(
                            details["D_hat"],
                            dtype=float,
                        ).reshape(-1)[0]
                    )

                    F_values = np.asarray(
                        details["F_matrix"],
                        dtype=float,
                    )[0]
                    weights = np.asarray(
                        details["kernel_weights"],
                        dtype=float,
                    ).reshape(-1)
                    xu_values = np.asarray(
                        xu,
                        dtype=float,
                    ).reshape(-1)

                    psi = (
                        xu_values
                        - x0[0]
                        - delta_t * a_hat
                    ) * F_values

                    conditional_mean = (
                        conditional_expectation(
                            psi,
                            weights,
                            f_hat,
                        )
                    )
                    conditional_second_moment = (
                        conditional_expectation(
                            psi**2,
                            weights,
                            f_hat,
                        )
                    )

                    variance_raw = (
                        conditional_second_moment
                        - conditional_mean**2
                    )
                    variance_hat = max(
                        float(variance_raw),
                        FLOOR,
                    )

                    sigma_raw = (
                        RK_EPANECHNIKOV_1D
                        / (
                            f_hat
                            * delta_t**2
                            * D_hat**2
                        )
                        * variance_hat
                    )
                    sigma_hat = max(
                        float(sigma_raw),
                        FLOOR,
                    )

                    Mh = float(M * h)
                    se_hat = math.sqrt(
                        sigma_hat / Mh
                    )
                    error = a_hat - a_star
                    Z = error / se_hat

                    half_width = 1.96 * se_hat
                    covered = int(
                        a_star >= a_hat - half_width
                        and a_star <= a_hat + half_width
                    )

                    sum_weights = float(
                        np.sum(weights)
                    )
                    sum_squared_weights = float(
                        np.sum(weights**2)
                    )
                    kernel_ess = (
                        sum_weights**2
                        / sum_squared_weights
                        if sum_squared_weights > 0.0
                        else 0.0
                    )

                    rows.append(
                        {
                            "model_id": model_id,
                            "config_path": str(
                                config_path
                            ),
                            "M": int(M),
                            "rep": int(rep),
                            "seed": int(
                                repetition_seed
                            ),
                            "alpha": float(alpha),
                            "c": float(constant),
                            "h": h,
                            "Mh": Mh,
                            "a_hat": a_hat,
                            "a_star": a_star,
                            "error": error,
                            "f_hat": f_hat,
                            "D_hat": D_hat,
                            "conditional_mean_hat":
                                conditional_mean,
                            "conditional_second_moment_hat":
                                conditional_second_moment,
                            "variance_raw": variance_raw,
                            "variance_hat": variance_hat,
                            "sigma_raw": sigma_raw,
                            "sigma_hat": sigma_hat,
                            "se_hat": se_hat,
                            "Z": Z,
                            "covered_95": covered,
                            "n_kernel_positive": int(
                                np.count_nonzero(
                                    weights > 0.0
                                )
                            ),
                            "kernel_weight_ess":
                                kernel_ess,
                            "f_floor_hit": int(
                                np.isclose(
                                    f_hat,
                                    FLOOR,
                                    rtol=0.0,
                                    atol=1e-15,
                                )
                            ),
                            "D_floor_hit": int(
                                np.isclose(
                                    D_hat,
                                    FLOOR,
                                    rtol=0.0,
                                    atol=1e-15,
                                )
                            ),
                            "variance_floor_hit": int(
                                variance_raw <= FLOOR
                            ),
                            "sigma_floor_hit": int(
                                sigma_raw <= FLOOR
                            ),
                        }
                    )

                if (
                    (rep + 1)
                    % max(1, args.progress_every)
                    == 0
                    or rep + 1 == args.reps
                ):
                    print(
                        f"[{model_id}] M={M} "
                        f"rep={rep + 1}/{args.reps}",
                        flush=True,
                    )

    raw = pd.DataFrame(rows)
    summary = summarize_rows(raw)
    aggregate = aggregate_summary(summary)

    raw.to_csv(
        output_directory / "paired_clt_rows.csv",
        index=False,
    )
    summary.to_csv(
        output_directory
        / "paired_clt_summary_by_M.csv",
        index=False,
    )
    aggregate.to_csv(
        output_directory
        / "paired_clt_aggregate_summary.csv",
        index=False,
    )

    comparison = pd.DataFrame()

    if args.compare_final:
        comparison = compare_with_final(raw)
        comparison.to_csv(
            output_directory
            / "final_run_reproduction.csv",
            index=False,
        )

    metadata = {
        "tag": tag,
        "configs": config_paths,
        "sample_sizes": sample_sizes,
        "reps": args.reps,
        "alphas": alphas,
        "constants": constants,
        "seed": args.seed,
        "shared_samples_across_bandwidths": True,
        "R_K": RK_EPANECHNIKOV_1D,
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
    print("=" * 120)
    print("AGGREGATE BANDWIDTH SENSITIVITY")
    print("=" * 120)

    with pd.option_context(
        "display.max_rows",
        200,
        "display.max_columns",
        40,
        "display.width",
        300,
    ):
        print(aggregate.to_string(index=False))

    if args.compare_final:
        print()
        print("=" * 120)
        print("FINAL-RUN REPRODUCTION")
        print("=" * 120)
        print(comparison.to_string(index=False))

    floor_events = int(
        raw[
            [
                "f_floor_hit",
                "D_floor_hit",
                "variance_floor_hit",
                "sigma_floor_hit",
            ]
        ]
        .to_numpy()
        .sum()
    )

    status_errors: list[str] = []

    if floor_events:
        status_errors.append(
            f"{floor_events} numerical floor events"
        )

    if args.compare_final and not comparison.empty:
        failed = comparison[
            comparison["status"] != "PASS"
        ]

        if not failed.empty:
            status_errors.append(
                "stored final run was not reproduced"
            )

    print()
    print("=" * 120)
    print("FINAL STATUS")
    print("=" * 120)

    if status_errors:
        print("FAIL")

        for error_message in status_errors:
            print(f"- {error_message}")

        raise SystemExit(1)

    print("PASS")
    print(f"Outputs: {output_directory}")


if __name__ == "__main__":
    main()
