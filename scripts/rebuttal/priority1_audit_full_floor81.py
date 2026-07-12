#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


RUN_ID = os.environ.get("RUN_ID")
if not RUN_ID:
    raise SystemExit("RUN_ID is not set.")

ROOT = Path("results")
OUT = ROOT / "rebuttal" / "priority1" / RUN_ID
OUT.mkdir(parents=True, exist_ok=True)

MODELS = {
    "GG1": ("gg_1d", 1),
    "GG2": ("gg_2d", 2),
    "MM1": ("mm_1d", 1),
    "MM2": ("mm_2d", 2),
}

FLOOR = 1.0e-12
TOL = 1.0e-7

REQUIRED = {
    "M",
    "rep",
    "seed",
    "h",
    "Mh_to_d",
    "sup_err",
    "ise",
    "f_hat",
    "f_floor_hit",
    "n_kernel_positive",
    "p_kernel_positive",
    "sum_kernel_weights",
    "kernel_weight_ess",
    "min_g1_hat",
    "median_g1_hat",
    "max_g1_hat",
    "min_D_hat",
    "q01_D_hat",
    "median_D_hat",
    "max_D_hat",
    "n_D_floor_grid",
    "p_D_floor_grid",
}

SELECTED_REQUIRED = REQUIRED - {"h"} | {
    "method",
    "selected_h",
    "boundary",
    "gap_to_oracle",
}


def q(series: pd.Series, probability: float) -> float:
    values = pd.to_numeric(series, errors="coerce")
    values = values[np.isfinite(values)]

    if values.empty:
        return np.nan

    return float(values.quantile(probability))


all_per_h: list[pd.DataFrame] = []
all_selected: list[pd.DataFrame] = []
candidate_summary_rows: list[dict[str, object]] = []
selected_summary_rows: list[dict[str, object]] = []
table4_rows: list[dict[str, object]] = []
errors: list[str] = []


for label, (model_id, dim) in MODELS.items():
    run_name = f"{model_id}_{RUN_ID}"
    directory = (
        ROOT
        / "raw"
        / "rate"
        / model_id
        / run_name
    )

    per_h_path = directory / "per_h_metrics.csv"
    selected_path = directory / "selected_metrics.csv"

    print("=" * 110)
    print(f"{label}: {directory}")
    print("=" * 110)

    if not per_h_path.exists():
        errors.append(f"{label}: missing {per_h_path}")
        continue

    if not selected_path.exists():
        errors.append(f"{label}: missing {selected_path}")
        continue

    per_h = pd.read_csv(per_h_path)
    selected = pd.read_csv(selected_path)

    missing_per_h = sorted(REQUIRED - set(per_h.columns))
    missing_selected = sorted(
        SELECTED_REQUIRED - set(selected.columns)
    )

    if missing_per_h:
        errors.append(
            f"{label}: missing per-h columns {missing_per_h}"
        )

    if missing_selected:
        errors.append(
            f"{label}: missing selected columns "
            f"{missing_selected}"
        )

    per_h.insert(0, "label", label)
    selected.insert(0, "label", label)

    all_per_h.append(per_h)
    all_selected.append(selected)

    # ----------------------------------------------------------
    # Basic validity checks
    # ----------------------------------------------------------
    finite_columns = [
        "h",
        "Mh_to_d",
        "sup_err",
        "ise",
        "f_hat",
        "kernel_weight_ess",
        "min_g1_hat",
        "min_D_hat",
        "median_D_hat",
        "max_D_hat",
    ]

    for column in finite_columns:
        values = pd.to_numeric(
            per_h[column],
            errors="coerce",
        )

        bad = int((~np.isfinite(values)).sum())

        if bad:
            errors.append(
                f"{label}: {bad} nonfinite values in {column}"
            )

    below_81 = per_h[
        per_h["Mh_to_d"] < 81.0 - TOL
    ]

    if not below_81.empty:
        errors.append(
            f"{label}: {len(below_81)} candidate rows "
            "have Mh^d < 81"
        )

    f_floor_count = int(per_h["f_floor_hit"].sum())
    D_floor_grid_count = int(
        per_h["n_D_floor_grid"].sum()
    )

    if f_floor_count:
        errors.append(
            f"{label}: {f_floor_count} f-floor events"
        )

    if D_floor_grid_count:
        errors.append(
            f"{label}: {D_floor_grid_count} "
            "D-floor grid points"
        )

    if (per_h["f_hat"] <= FLOOR + 1e-15).any():
        errors.append(
            f"{label}: f_hat reaches the numerical floor"
        )

    if (per_h["min_D_hat"] <= FLOOR + 1e-15).any():
        errors.append(
            f"{label}: min_D_hat reaches the numerical floor"
        )

    # ----------------------------------------------------------
    # Selected/per-bandwidth consistency
    # ----------------------------------------------------------
    diagnostic_columns = [
        "Mh_to_d",
        "f_hat",
        "f_floor_hit",
        "n_kernel_positive",
        "p_kernel_positive",
        "sum_kernel_weights",
        "kernel_weight_ess",
        "min_g1_hat",
        "median_g1_hat",
        "max_g1_hat",
        "min_D_hat",
        "q01_D_hat",
        "median_D_hat",
        "max_D_hat",
        "n_D_floor_grid",
        "p_D_floor_grid",
    ]

    mismatch_count = 0

    for row in selected.itertuples(index=False):
        match = per_h[
            (per_h["M"] == row.M)
            & (per_h["rep"] == row.rep)
            & np.isclose(
                per_h["h"],
                row.selected_h,
                rtol=0.0,
                atol=1e-12,
            )
        ]

        if len(match) != 1:
            errors.append(
                f"{label}: selected/per-h match count "
                f"is {len(match)} for M={row.M}, "
                f"rep={row.rep}, method={row.method}"
            )
            continue

        match_row = match.iloc[0]

        for column in diagnostic_columns:
            left = float(getattr(row, column))
            right = float(match_row[column])

            if not np.isclose(
                left,
                right,
                rtol=1e-12,
                atol=1e-15,
                equal_nan=True,
            ):
                mismatch_count += 1

    if mismatch_count:
        errors.append(
            f"{label}: {mismatch_count} selected/per-h "
            "diagnostic mismatches"
        )

    # ----------------------------------------------------------
    # Candidate-grid summaries
    # ----------------------------------------------------------
    for M, group in per_h.groupby("M", sort=True):
        candidate_summary_rows.append(
            {
                "label": label,
                "model_id": model_id,
                "dim": dim,
                "M": int(M),
                "n_rows": len(group),
                "n_reps": int(group["rep"].nunique()),
                "n_bandwidths": int(group["h"].nunique()),
                "min_h": float(group["h"].min()),
                "max_h": float(group["h"].max()),
                "min_Mh_to_d": float(
                    group["Mh_to_d"].min()
                ),
                "min_f_hat": float(
                    group["f_hat"].min()
                ),
                "q01_f_hat": q(group["f_hat"], 0.01),
                "median_f_hat": q(group["f_hat"], 0.50),
                "min_D_hat": float(
                    group["min_D_hat"].min()
                ),
                "q01_D_hat": q(
                    group["min_D_hat"],
                    0.01,
                ),
                "q05_D_hat": q(
                    group["min_D_hat"],
                    0.05,
                ),
                "median_min_D_hat": q(
                    group["min_D_hat"],
                    0.50,
                ),
                "min_kernel_positive": int(
                    group["n_kernel_positive"].min()
                ),
                "q01_kernel_positive": q(
                    group["n_kernel_positive"],
                    0.01,
                ),
                "min_kernel_ess": float(
                    group["kernel_weight_ess"].min()
                ),
                "q01_kernel_ess": q(
                    group["kernel_weight_ess"],
                    0.01,
                ),
                "median_kernel_ess": q(
                    group["kernel_weight_ess"],
                    0.50,
                ),
                "n_f_floor": int(
                    group["f_floor_hit"].sum()
                ),
                "n_D_floor_grid": int(
                    group["n_D_floor_grid"].sum()
                ),
                "max_sup_err": float(
                    group["sup_err"].max()
                ),
            }
        )

    # ----------------------------------------------------------
    # Selected summaries
    # ----------------------------------------------------------
    for (method, M), group in selected.groupby(
        ["method", "M"],
        sort=True,
    ):
        selected_summary_rows.append(
            {
                "label": label,
                "model_id": model_id,
                "method": method,
                "M": int(M),
                "n": len(group),
                "mean_selected_h": float(
                    group["selected_h"].mean()
                ),
                "min_Mh_to_d": float(
                    group["Mh_to_d"].min()
                ),
                "boundary_frequency": float(
                    group["boundary"].mean()
                ),
                "mean_gap_to_oracle": float(
                    group["gap_to_oracle"].mean()
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
                "max_sup_err": float(
                    group["sup_err"].max()
                ),
                "n_f_floor": int(
                    group["f_floor_hit"].sum()
                ),
                "n_D_floor_grid": int(
                    group["n_D_floor_grid"].sum()
                ),
            }
        )

    adaptive = selected[
        selected["method"].astype(str).str.lower()
        == "lepski"
    ]
    oracle = selected[
        selected["method"].astype(str).str.lower()
        == "oracle"
    ]

    gaps_by_M = adaptive.groupby("M")[
        "gap_to_oracle"
    ].mean()

    table4_rows.append(
        {
            "testbed": label,
            "oracle_bandwidth": float(
                oracle["selected_h"].mean()
            ),
            "adaptive_bandwidth": float(
                adaptive["selected_h"].mean()
            ),
            "boundary_frequency": float(
                adaptive["boundary"].mean()
            ),
            "Cavg": float(gaps_by_M.mean()),
            "Cmax": float(gaps_by_M.max()),
            "min_candidate_f_hat": float(
                per_h["f_hat"].min()
            ),
            "min_candidate_D_hat": float(
                per_h["min_D_hat"].min()
            ),
            "min_candidate_kernel_ess": float(
                per_h["kernel_weight_ess"].min()
            ),
            "f_floor_events": int(
                per_h["f_floor_hit"].sum()
            ),
            "D_floor_grid_points": int(
                per_h["n_D_floor_grid"].sum()
            ),
            "source_run_tag": run_name,
        }
    )

    print(
        f"candidate rows={len(per_h):,}, "
        f"selected rows={len(selected):,}, "
        f"min Mh^d={per_h['Mh_to_d'].min():.9f}, "
        f"min f_hat={per_h['f_hat'].min():.6e}, "
        f"min D_hat={per_h['min_D_hat'].min():.6e}, "
        f"min kernel ESS="
        f"{per_h['kernel_weight_ess'].min():.3f}, "
        f"f floors={f_floor_count}, "
        f"D-floor grid points={D_floor_grid_count}"
    )


per_h_all = pd.concat(
    all_per_h,
    ignore_index=True,
    sort=False,
)
selected_all = pd.concat(
    all_selected,
    ignore_index=True,
    sort=False,
)

candidate_summary = pd.DataFrame(
    candidate_summary_rows
)
selected_summary = pd.DataFrame(
    selected_summary_rows
)
table4 = pd.DataFrame(table4_rows)

per_h_all.to_csv(
    OUT / "all_candidate_bandwidth_rows.csv",
    index=False,
)
selected_all.to_csv(
    OUT / "all_selected_rows.csv",
    index=False,
)
candidate_summary.to_csv(
    OUT / "candidate_denominator_summary.csv",
    index=False,
)
selected_summary.to_csv(
    OUT / "selected_denominator_summary.csv",
    index=False,
)
table4.to_csv(
    OUT / "table4_floor81_with_denominators.csv",
    index=False,
)


# --------------------------------------------------------------
# Compare against the prior consistent table.
# --------------------------------------------------------------
reference_path = (
    ROOT
    / "rebuttal"
    / "day4"
    / "table4_floor81_consistent.csv"
)

comparison_rows: list[dict[str, object]] = []

if reference_path.exists():
    reference = pd.read_csv(reference_path)

    for row in table4.itertuples(index=False):
        ref = reference[
            reference["testbed"] == row.testbed
        ]

        if len(ref) != 1:
            continue

        ref_row = ref.iloc[0]

        comparison_rows.append(
            {
                "testbed": row.testbed,
                "oracle_bandwidth_new":
                    row.oracle_bandwidth,
                "oracle_bandwidth_old":
                    float(ref_row["oracle_bandwidth"]),
                "adaptive_bandwidth_new":
                    row.adaptive_bandwidth,
                "adaptive_bandwidth_old":
                    float(ref_row["adaptive_bandwidth"]),
                "boundary_frequency_new":
                    row.boundary_frequency,
                "boundary_frequency_old":
                    float(ref_row["boundary_frequency"]),
                "Cavg_new": row.Cavg,
                "Cavg_old": float(ref_row["Cavg"]),
                "Cmax_new": row.Cmax,
                "Cmax_old": float(ref_row["Cmax"]),
            }
        )

comparison = pd.DataFrame(comparison_rows)

if not comparison.empty:
    for name in [
        "oracle_bandwidth",
        "adaptive_bandwidth",
        "boundary_frequency",
        "Cavg",
        "Cmax",
    ]:
        comparison[f"{name}_abs_diff"] = np.abs(
            comparison[f"{name}_new"]
            - comparison[f"{name}_old"]
        )

    comparison.to_csv(
        OUT / "table4_reproduction_comparison.csv",
        index=False,
    )


print()
print("=" * 110)
print("CANDIDATE-GRID DENOMINATOR SUMMARY")
print("=" * 110)

with pd.option_context(
    "display.max_rows",
    100,
    "display.max_columns",
    40,
    "display.width",
    280,
):
    print(candidate_summary.to_string(index=False))


print()
print("=" * 110)
print("REGENERATED FLOOR-81 TABLE")
print("=" * 110)
print(table4.to_string(index=False))


print()
print("=" * 110)
print("COMPARISON WITH PRIOR CONSISTENT TABLE")
print("=" * 110)

if comparison.empty:
    print("No reference comparison available.")
else:
    print(comparison.to_string(index=False))


print()
print("=" * 110)
print("FINAL STATUS")
print("=" * 110)

if errors:
    print("FAIL")
    for error in errors:
        print(f"- {error}")
    raise SystemExit(1)

print("PASS")
print(f"Outputs saved under: {OUT}")
