#!/usr/bin/env python3
"""Priority 5: finite-sample rate and slope robustness audit.

This audit reuses saved per-repetition selected_metrics.csv files. It does not
run new simulations.

Checks:
  * exact reproduction of saved slopes.csv values;
  * paper-facing bandwidth grids versus floor-81 grids;
  * all-four, early-three, and late-three sample-size windows;
  * mean versus median aggregation;
  * sup-grid error versus ISE;
  * nonparametric bootstrap uncertainty over Monte Carlo repetitions;
  * selected-bandwidth boundary frequencies;
  * paper-number provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

MODEL_ORDER = ["gg_1d", "gg_2d", "mm_1d", "mm_2d"]

PAPER_EXPECTED_ORACLE_SUP_SLOPES = {
    "gg_1d": -0.3488994864487099,
    "gg_2d": -0.24343523241768025,
    "mm_1d": -0.3902290749462957,
    "mm_2d": -0.38582814230399737,
}


@dataclass(frozen=True)
class RunSpec:
    model_id: str
    protocol: str
    raw_dir: str
    processed_dir: str


RUN_SPECS = [
    RunSpec(
        model_id="gg_1d",
        protocol="paper_grid",
        raw_dir="results/raw/rate/gg_1d/gg1_final_rawmax_k2",
        processed_dir="results/processed/rate/gg_1d/gg1_final_rawmax_k2",
    ),
    RunSpec(
        model_id="gg_2d",
        protocol="paper_grid",
        raw_dir="results/raw/rate/gg_2d/gg2_final_rawmax_k2",
        processed_dir="results/processed/rate/gg_2d/gg2_final_rawmax_k2",
    ),
    RunSpec(
        model_id="mm_1d",
        protocol="paper_grid",
        raw_dir="results/raw/rate/mm_1d/mm1_final_rawmax_k2",
        processed_dir="results/processed/rate/mm_1d/mm1_final_rawmax_k2",
    ),
    RunSpec(
        model_id="mm_2d",
        protocol="paper_grid",
        raw_dir="results/raw/rate/mm_2d/mm2_final_rawmax_k2",
        processed_dir="results/processed/rate/mm_2d/mm2_final_rawmax_k2",
    ),
    RunSpec(
        model_id="gg_1d",
        protocol="floor81",
        raw_dir=(
            "results/raw/rate/gg_1d/"
            "gg_1d_priority1_floor81_main_20260711_180209"
        ),
        processed_dir=(
            "results/processed/rate/gg_1d/"
            "gg_1d_priority1_floor81_main_20260711_180209"
        ),
    ),
    RunSpec(
        model_id="gg_2d",
        protocol="floor81",
        raw_dir=(
            "results/raw/rate/gg_2d/"
            "gg_2d_priority1_floor81_main_20260711_180209"
        ),
        processed_dir=(
            "results/processed/rate/gg_2d/"
            "gg_2d_priority1_floor81_main_20260711_180209"
        ),
    ),
    RunSpec(
        model_id="mm_1d",
        protocol="floor81",
        raw_dir=(
            "results/raw/rate/mm_1d/"
            "mm_1d_priority1_floor81_main_20260711_180209"
        ),
        processed_dir=(
            "results/processed/rate/mm_1d/"
            "mm_1d_priority1_floor81_main_20260711_180209"
        ),
    ),
    RunSpec(
        model_id="mm_2d",
        protocol="floor81",
        raw_dir=(
            "results/raw/rate/mm_2d/"
            "mm_2d_priority1_floor81_main_20260711_180209"
        ),
        processed_dir=(
            "results/processed/rate/mm_2d/"
            "mm_2d_priority1_floor81_main_20260711_180209"
        ),
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit finite-sample rate-slope robustness."
    )
    parser.add_argument(
        "--bootstrap-reps",
        type=int,
        default=5000,
        help="Number of within-M bootstrap resamples.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260712,
        help="Master seed used only for bootstrap resampling.",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default="priority5_rate_slope",
        help="Output directory tag.",
    )
    return parser.parse_args()


def deterministic_seed(master_seed: int, *parts: Any) -> int:
    text = "|".join([str(master_seed), *map(str, parts)])
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little") % (2**32)


def finite_range_theory_slope(
    sample_sizes: list[int],
    beta: float,
    dim: int,
) -> float:
    m1 = min(sample_sizes)
    m2 = max(sample_sizes)

    if m1 == m2:
        return float("nan")

    p = beta / (2.0 * beta + dim)

    return float(
        -p
        + p
        * math.log(math.log(m2) / math.log(m1))
        / math.log(m2 / m1)
    )


def asymptotic_theory_slope(beta: float, dim: int) -> float:
    return float(-beta / (2.0 * beta + dim))


def fit_slope(
    sample_sizes: np.ndarray,
    aggregated_errors: np.ndarray,
) -> tuple[float, float, float]:
    sample_sizes = np.asarray(sample_sizes, dtype=float)
    aggregated_errors = np.asarray(aggregated_errors, dtype=float)

    if len(sample_sizes) < 2:
        return float("nan"), float("nan"), float("nan")

    if np.any(~np.isfinite(aggregated_errors)):
        return float("nan"), float("nan"), float("nan")

    if np.any(aggregated_errors <= 0.0):
        return float("nan"), float("nan"), float("nan")

    x = np.log(sample_sizes)
    y = np.log(aggregated_errors)

    slope, intercept = np.polyfit(x, y, deg=1)
    fitted = slope * x + intercept

    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))

    if ss_tot == 0.0:
        r_squared = 1.0 if ss_res == 0.0 else float("nan")
    else:
        r_squared = 1.0 - ss_res / ss_tot

    return float(slope), float(intercept), float(r_squared)


def bootstrap_slopes(
    frame: pd.DataFrame,
    metric: str,
    sample_sizes: list[int],
    aggregation: str,
    bootstrap_reps: int,
    rng: np.random.Generator,
) -> np.ndarray:
    x = np.log(np.asarray(sample_sizes, dtype=float))
    x_centered = x - np.mean(x)
    denominator = float(np.sum(x_centered**2))

    boot_aggregates: list[np.ndarray] = []

    for sample_size in sample_sizes:
        values = frame.loc[
            frame["M"] == sample_size,
            metric,
        ].to_numpy(dtype=float)

        if len(values) == 0:
            raise ValueError(
                f"No values for M={sample_size}, metric={metric}"
            )

        indices = rng.integers(
            0,
            len(values),
            size=(bootstrap_reps, len(values)),
        )
        samples = values[indices]

        if aggregation == "mean":
            aggregate = np.mean(samples, axis=1)
        elif aggregation == "median":
            aggregate = np.median(samples, axis=1)
        else:
            raise ValueError(f"Unknown aggregation: {aggregation}")

        boot_aggregates.append(aggregate)

    aggregate_matrix = np.column_stack(boot_aggregates)

    if np.any(aggregate_matrix <= 0.0):
        raise ValueError("Bootstrap produced a nonpositive error aggregate.")

    y = np.log(aggregate_matrix)
    y_centered = y - np.mean(y, axis=1, keepdims=True)

    slopes = (y_centered @ x_centered) / denominator
    return np.asarray(slopes, dtype=float)


def aggregation_function(name: str):
    if name == "mean":
        return np.mean
    if name == "median":
        return np.median
    raise ValueError(name)


def get_windows(sample_sizes: list[int]) -> dict[str, list[int]]:
    sample_sizes = sorted(sample_sizes)

    windows = {
        "all": sample_sizes,
    }

    if len(sample_sizes) >= 3:
        windows["early3"] = sample_sizes[:3]
        windows["late3"] = sample_sizes[-3:]

    return windows


def read_run(
    spec: RunSpec,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    raw_dir = ROOT / spec.raw_dir
    processed_dir = ROOT / spec.processed_dir

    selected_path = raw_dir / "selected_metrics.csv"
    slopes_path = processed_dir / "slopes.csv"
    config_path = processed_dir / "run_config.json"

    for path in [selected_path, slopes_path, config_path]:
        if not path.is_file():
            raise FileNotFoundError(path)

    selected = pd.read_csv(selected_path)
    stored_slopes = pd.read_csv(slopes_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))

    return selected, stored_slopes, config


def validate_selected_metrics(
    frame: pd.DataFrame,
    config: dict[str, Any],
    spec: RunSpec,
) -> list[str]:
    errors: list[str] = []

    required_columns = {
        "model_id",
        "M",
        "rep",
        "seed",
        "method",
        "selected_h",
        "sup_err",
        "ise",
        "boundary",
    }

    missing = required_columns - set(frame.columns)

    if missing:
        errors.append(
            f"{spec.model_id}/{spec.protocol}: missing columns "
            f"{sorted(missing)}"
        )
        return errors

    if set(frame["model_id"].astype(str).unique()) != {spec.model_id}:
        errors.append(
            f"{spec.model_id}/{spec.protocol}: unexpected model_id values"
        )

    expected_sample_sizes = sorted(
        int(value) for value in config["sample_sizes"]
    )
    observed_sample_sizes = sorted(
        int(value) for value in frame["M"].unique()
    )

    if expected_sample_sizes != observed_sample_sizes:
        errors.append(
            f"{spec.model_id}/{spec.protocol}: sample-size mismatch; "
            f"expected={expected_sample_sizes}, "
            f"observed={observed_sample_sizes}"
        )

    methods = set(frame["method"].astype(str).unique())

    if methods != {"lepski", "oracle"}:
        errors.append(
            f"{spec.model_id}/{spec.protocol}: unexpected methods {methods}"
        )

    duplicate_count = int(
        frame.duplicated(["M", "rep", "method"]).sum()
    )

    if duplicate_count:
        errors.append(
            f"{spec.model_id}/{spec.protocol}: "
            f"{duplicate_count} duplicate M/rep/method rows"
        )

    expected_reps = int(config["reps"])

    counts = (
        frame.groupby(["M", "method"])
        .size()
        .rename("rows")
        .reset_index()
    )

    bad_counts = counts[counts["rows"] != expected_reps]

    if not bad_counts.empty:
        errors.append(
            f"{spec.model_id}/{spec.protocol}: incorrect replicate counts"
        )

    for metric in ["sup_err", "ise"]:
        values = frame[metric].to_numpy(dtype=float)

        if np.any(~np.isfinite(values)):
            errors.append(
                f"{spec.model_id}/{spec.protocol}: "
                f"nonfinite {metric} values"
            )

        if np.any(values <= 0.0):
            errors.append(
                f"{spec.model_id}/{spec.protocol}: "
                f"nonpositive {metric} values"
            )

    boundary_values = set(
        frame["boundary"].dropna().astype(float).unique().tolist()
    )

    if not boundary_values.issubset({0.0, 1.0}):
        errors.append(
            f"{spec.model_id}/{spec.protocol}: "
            f"unexpected boundary values {boundary_values}"
        )

    return errors


def main() -> int:
    args = parse_args()

    if args.bootstrap_reps < 100:
        raise ValueError("--bootstrap-reps must be at least 100")

    output_dir = (
        ROOT
        / "results"
        / "rebuttal"
        / "priority5"
        / args.tag
    )
    output_dir.mkdir(parents=True, exist_ok=False)

    all_errors: list[str] = []
    slope_rows: list[dict[str, Any]] = []
    reproduction_rows: list[dict[str, Any]] = []
    boundary_rows: list[dict[str, Any]] = []
    inventory_rows: list[dict[str, Any]] = []

    for spec in RUN_SPECS:
        selected, stored_slopes, config = read_run(spec)

        validation_errors = validate_selected_metrics(
            selected,
            config,
            spec,
        )
        all_errors.extend(validation_errors)

        dim = int(config.get("dim", selected["dim"].iloc[0]))
        beta = float(config.get("beta_effective", 2.0))
        sample_sizes = sorted(
            int(value) for value in selected["M"].unique()
        )
        windows = get_windows(sample_sizes)

        inventory_rows.append(
            {
                "model_id": spec.model_id,
                "protocol": spec.protocol,
                "raw_dir": spec.raw_dir,
                "processed_dir": spec.processed_dir,
                "rows": len(selected),
                "dim": dim,
                "beta_effective": beta,
                "reps": int(config["reps"]),
                "sample_sizes": ",".join(map(str, sample_sizes)),
                "min_h_factor": float(config["min_h_factor"]),
                "h0": float(config["h0"]),
                "q": float(config["q"]),
                "selector_metric": str(config["selector_metric"]),
                "penalty_form": str(config["penalty_form"]),
                "kappa_pair": float(config["kappa_pair"]),
                "kappa_final": float(config["kappa_final"]),
                "seed": int(config["seed"]),
            }
        )

        for (method, sample_size), group in selected.groupby(
            ["method", "M"],
            sort=True,
        ):
            boundary_rows.append(
                {
                    "model_id": spec.model_id,
                    "protocol": spec.protocol,
                    "method": method,
                    "M": int(sample_size),
                    "boundary_rate": float(
                        group["boundary"].astype(float).mean()
                    ),
                    "n_rows": len(group),
                    "mean_selected_h": float(
                        group["selected_h"].astype(float).mean()
                    ),
                    "median_selected_h": float(
                        group["selected_h"].astype(float).median()
                    ),
                }
            )

        for method in ["lepski", "oracle"]:
            method_frame = selected[
                selected["method"] == method
            ].copy()

            for metric in ["sup_err", "ise"]:
                for aggregation in ["mean", "median"]:
                    aggregate_fn = aggregation_function(aggregation)

                    for window_name, window_sizes in windows.items():
                        point_values = []

                        for sample_size in window_sizes:
                            values = method_frame.loc[
                                method_frame["M"] == sample_size,
                                metric,
                            ].to_numpy(dtype=float)

                            point_values.append(
                                float(aggregate_fn(values))
                            )

                        slope, intercept, r_squared = fit_slope(
                            np.asarray(window_sizes, dtype=float),
                            np.asarray(point_values, dtype=float),
                        )

                        theory_finite = finite_range_theory_slope(
                            window_sizes,
                            beta,
                            dim,
                        )
                        theory_asymptotic = (
                            asymptotic_theory_slope(beta, dim)
                        )

                        seed = deterministic_seed(
                            args.seed,
                            spec.model_id,
                            spec.protocol,
                            method,
                            metric,
                            aggregation,
                            window_name,
                        )
                        rng = np.random.default_rng(seed)

                        bootstrap = bootstrap_slopes(
                            frame=method_frame,
                            metric=metric,
                            sample_sizes=window_sizes,
                            aggregation=aggregation,
                            bootstrap_reps=args.bootstrap_reps,
                            rng=rng,
                        )

                        ci_low, ci_high = np.quantile(
                            bootstrap,
                            [0.025, 0.975],
                        )

                        monotonic_increases = int(
                            np.sum(np.diff(point_values) > 0.0)
                        )

                        slope_rows.append(
                            {
                                "model_id": spec.model_id,
                                "protocol": spec.protocol,
                                "method": method,
                                "metric": metric,
                                "aggregation": aggregation,
                                "window": window_name,
                                "sample_sizes": ",".join(
                                    map(str, window_sizes)
                                ),
                                "n_sample_sizes": len(window_sizes),
                                "slope": slope,
                                "intercept": intercept,
                                "r_squared": r_squared,
                                "bootstrap_mean": float(
                                    np.mean(bootstrap)
                                ),
                                "bootstrap_sd": float(
                                    np.std(bootstrap, ddof=1)
                                ),
                                "ci95_low": float(ci_low),
                                "ci95_high": float(ci_high),
                                "theory_finite": theory_finite,
                                "theory_asymptotic": theory_asymptotic,
                                "slope_minus_finite_theory": (
                                    slope - theory_finite
                                ),
                                "finite_theory_in_ci95": bool(
                                    ci_low
                                    <= theory_finite
                                    <= ci_high
                                ),
                                "asymptotic_theory_in_ci95": bool(
                                    ci_low
                                    <= theory_asymptotic
                                    <= ci_high
                                ),
                                "monotonic_increases": (
                                    monotonic_increases
                                ),
                                "bootstrap_reps": args.bootstrap_reps,
                                "bootstrap_seed": seed,
                            }
                        )

        computed = pd.DataFrame(slope_rows)

        current = computed[
            (computed["model_id"] == spec.model_id)
            & (computed["protocol"] == spec.protocol)
            & (computed["aggregation"] == "mean")
            & (computed["window"] == "all")
        ]

        for _, stored_row in stored_slopes.iterrows():
            method = str(stored_row["method"])
            metric = str(stored_row["metric"])

            match = current[
                (current["method"] == method)
                & (current["metric"] == metric)
            ]

            if len(match) != 1:
                all_errors.append(
                    f"{spec.model_id}/{spec.protocol}: "
                    f"could not uniquely match stored slope "
                    f"{method}/{metric}"
                )
                continue

            calculated_slope = float(match.iloc[0]["slope"])
            stored_slope = float(stored_row["slope"])
            absolute_difference = abs(
                calculated_slope - stored_slope
            )
            status = (
                "PASS"
                if absolute_difference <= 1e-12
                else "FAIL"
            )

            if status == "FAIL":
                all_errors.append(
                    f"{spec.model_id}/{spec.protocol}: "
                    f"slope reproduction failed for "
                    f"{method}/{metric}: "
                    f"{absolute_difference:.3e}"
                )

            reproduction_rows.append(
                {
                    "model_id": spec.model_id,
                    "protocol": spec.protocol,
                    "method": method,
                    "metric": metric,
                    "stored_slope": stored_slope,
                    "calculated_slope": calculated_slope,
                    "absolute_difference": absolute_difference,
                    "status": status,
                }
            )

    slope_df = pd.DataFrame(slope_rows)
    reproduction_df = pd.DataFrame(reproduction_rows)
    boundary_df = pd.DataFrame(boundary_rows)
    inventory_df = pd.DataFrame(inventory_rows)

    paper_oracle = slope_df[
        (slope_df["protocol"] == "paper_grid")
        & (slope_df["method"] == "oracle")
        & (slope_df["metric"] == "sup_err")
        & (slope_df["aggregation"] == "mean")
        & (slope_df["window"] == "all")
    ].copy()

    floor_oracle = slope_df[
        (slope_df["protocol"] == "floor81")
        & (slope_df["method"] == "oracle")
        & (slope_df["metric"] == "sup_err")
        & (slope_df["aggregation"] == "mean")
        & (slope_df["window"] == "all")
    ].copy()

    paper_columns = {
        "slope": "paper_slope",
        "ci95_low": "paper_ci95_low",
        "ci95_high": "paper_ci95_high",
        "theory_finite": "theory_finite",
        "r_squared": "paper_r_squared",
    }
    floor_columns = {
        "slope": "floor81_slope",
        "ci95_low": "floor81_ci95_low",
        "ci95_high": "floor81_ci95_high",
        "r_squared": "floor81_r_squared",
    }

    comparison = (
        paper_oracle[
            ["model_id", *paper_columns.keys()]
        ]
        .rename(columns=paper_columns)
        .merge(
            floor_oracle[
                ["model_id", *floor_columns.keys()]
            ].rename(columns=floor_columns),
            on="model_id",
            how="outer",
            validate="one_to_one",
        )
    )

    comparison["floor81_minus_paper"] = (
        comparison["floor81_slope"]
        - comparison["paper_slope"]
    )
    comparison["absolute_protocol_difference"] = (
        comparison["floor81_minus_paper"].abs()
    )
    comparison["paper_expected_slope"] = comparison[
        "model_id"
    ].map(PAPER_EXPECTED_ORACLE_SUP_SLOPES)
    comparison["paper_expected_abs_diff"] = (
        comparison["paper_slope"]
        - comparison["paper_expected_slope"]
    ).abs()
    comparison["paper_number_status"] = np.where(
        comparison["paper_expected_abs_diff"] <= 1e-12,
        "PASS",
        "FAIL",
    )

    failed_paper = comparison[
        comparison["paper_number_status"] != "PASS"
    ]

    if not failed_paper.empty:
        for _, row in failed_paper.iterrows():
            all_errors.append(
                f"{row['model_id']}: paper slope does not reproduce "
                f"expected paper number"
            )

    model_rank = {
        model_id: index
        for index, model_id in enumerate(MODEL_ORDER)
    }

    for frame in [
        slope_df,
        reproduction_df,
        boundary_df,
        inventory_df,
        comparison,
    ]:
        if "model_id" in frame.columns:
            frame["_model_rank"] = frame["model_id"].map(model_rank)
            frame.sort_values(
                [column for column in [
                    "_model_rank",
                    "protocol",
                    "method",
                    "metric",
                    "aggregation",
                    "window",
                    "M",
                ] if column in frame.columns],
                inplace=True,
            )
            frame.drop(columns=["_model_rank"], inplace=True)

    slope_df.to_csv(
        output_dir / "slope_robustness.csv",
        index=False,
    )
    reproduction_df.to_csv(
        output_dir / "stored_slope_reproduction.csv",
        index=False,
    )
    boundary_df.to_csv(
        output_dir / "boundary_rates.csv",
        index=False,
    )
    inventory_df.to_csv(
        output_dir / "run_inventory.csv",
        index=False,
    )
    comparison.to_csv(
        output_dir / "paper_vs_floor81_oracle_sup.csv",
        index=False,
    )

    run_config = {
        "bootstrap_reps": args.bootstrap_reps,
        "seed": args.seed,
        "tag": args.tag,
        "run_specs": [asdict(spec) for spec in RUN_SPECS],
        "paper_expected_oracle_sup_slopes": (
            PAPER_EXPECTED_ORACLE_SUP_SLOPES
        ),
    }
    (output_dir / "run_config.json").write_text(
        json.dumps(run_config, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print()
    print("=" * 136)
    print("PAPER GRID VS FLOOR-81: ORACLE SUP-ERROR SLOPES")
    print("=" * 136)

    comparison_columns = [
        "model_id",
        "paper_slope",
        "paper_ci95_low",
        "paper_ci95_high",
        "floor81_slope",
        "floor81_ci95_low",
        "floor81_ci95_high",
        "floor81_minus_paper",
        "theory_finite",
        "paper_number_status",
    ]
    print(
        comparison[comparison_columns].to_string(
            index=False,
            float_format=lambda value: f"{value:.6f}",
        )
    )

    print()
    print("=" * 136)
    print("ORACLE SUP-ERROR SENSITIVITY")
    print("=" * 136)

    sensitivity = slope_df[
        (slope_df["method"] == "oracle")
        & (slope_df["metric"] == "sup_err")
    ][
        [
            "model_id",
            "protocol",
            "aggregation",
            "window",
            "sample_sizes",
            "slope",
            "ci95_low",
            "ci95_high",
            "r_squared",
            "theory_finite",
            "slope_minus_finite_theory",
            "finite_theory_in_ci95",
            "monotonic_increases",
        ]
    ]

    print(
        sensitivity.to_string(
            index=False,
            float_format=lambda value: f"{value:.6f}",
        )
    )

    print()
    print("=" * 136)
    print("MAXIMUM BOUNDARY RATES")
    print("=" * 136)

    maximum_boundary = (
        boundary_df.groupby(
            ["model_id", "protocol", "method"],
            as_index=False,
        )
        .agg(
            max_boundary_rate=("boundary_rate", "max"),
            mean_boundary_rate=("boundary_rate", "mean"),
        )
    )

    print(
        maximum_boundary.to_string(
            index=False,
            float_format=lambda value: f"{value:.6f}",
        )
    )

    print()
    print("=" * 136)
    print("STORED-SLOPE REPRODUCTION")
    print("=" * 136)

    print(
        reproduction_df.to_string(
            index=False,
            float_format=lambda value: f"{value:.12e}",
        )
    )

    status = "PASS" if not all_errors else "FAIL"

    summary_lines = [
        "# Priority 5 rate-slope robustness audit",
        "",
        f"- status: `{status}`",
        f"- bootstrap repetitions: `{args.bootstrap_reps}`",
        f"- bootstrap seed: `{args.seed}`",
        "- simulations rerun: `no`",
        "- source: saved per-repetition selected_metrics.csv files",
        "",
        "## Paper-grid versus floor-81 oracle sup-error slopes",
        "",
        "```text",
        comparison[comparison_columns].to_string(index=False),
        "```",
        "",
        "## Validation errors",
        "",
    ]

    if all_errors:
        summary_lines.extend(
            f"- {error}" for error in all_errors
        )
    else:
        summary_lines.append("- none")

    (output_dir / "summary.md").write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    print()
    print("=" * 136)
    print("FINAL STATUS")
    print("=" * 136)
    print(status)

    if all_errors:
        print()
        print("Validation errors:")
        for error in all_errors:
            print(f"- {error}")

    print(f"Outputs: {output_dir}")

    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
