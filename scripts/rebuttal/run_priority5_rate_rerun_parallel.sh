#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT" || exit 1

RATE_SCRIPT="scripts/01_rate.py"
LOG_ROOT="results/rebuttal/priority5_rerun_logs"
BASE="priority5_rerun_$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${LOG_ROOT}/${BASE}"
MANIFEST="${LOG_DIR}/manifest.tsv"

mkdir -p "$LOG_DIR"
printf '%s\n' "$BASE" > "${LOG_ROOT}/LATEST"

if [[ ! -s "$RATE_SCRIPT" ]]; then
    echo "ERROR: missing rate script: $RATE_SCRIPT" >&2
    exit 1
fi

python -m py_compile "$RATE_SCRIPT" || exit 1

cat > "$MANIFEST" <<EOF
label	model_id	protocol	config	reps	min_h_factor	tag	reference_processed
gg_1d_paper	gg_1d	paper_grid	configs/gg_1d.yaml	50	20	${BASE}_gg_1d_paper	results/processed/rate/gg_1d/gg1_final_rawmax_k2
gg_1d_floor81	gg_1d	floor81	configs/gg_1d.yaml	50	81	${BASE}_gg_1d_floor81	results/processed/rate/gg_1d/gg_1d_priority1_floor81_main_20260711_180209
mm_1d_paper	mm_1d	paper_grid	configs/mm_1d.yaml	50	20	${BASE}_mm_1d_paper	results/processed/rate/mm_1d/mm1_final_rawmax_k2
mm_1d_floor81	mm_1d	floor81	configs/mm_1d.yaml	50	81	${BASE}_mm_1d_floor81	results/processed/rate/mm_1d/mm_1d_priority1_floor81_main_20260711_180209
gg_2d_paper	gg_2d	paper_grid	configs/gg_2d.yaml	20	5	${BASE}_gg_2d_paper	results/processed/rate/gg_2d/gg2_final_rawmax_k2
gg_2d_floor81	gg_2d	floor81	configs/gg_2d.yaml	20	9	${BASE}_gg_2d_floor81	results/processed/rate/gg_2d/gg_2d_priority1_floor81_main_20260711_180209
mm_2d_paper	mm_2d	paper_grid	configs/mm_2d.yaml	20	5	${BASE}_mm_2d_paper	results/processed/rate/mm_2d/mm2_final_rawmax_k2
mm_2d_floor81	mm_2d	floor81	configs/mm_2d.yaml	20	9	${BASE}_mm_2d_floor81	results/processed/rate/mm_2d/mm_2d_priority1_floor81_main_20260711_180209
EOF

run_one() {
    local label="$1"
    local model_id="$2"
    local config="$3"
    local reps="$4"
    local min_h_factor="$5"
    local tag="$6"

    local log="${LOG_DIR}/${label}.log"
    local status_file="${LOG_DIR}/${label}.status"

    date --iso-8601=seconds > "${LOG_DIR}/${label}.started"

    (
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        export VECLIB_MAXIMUM_THREADS=1
        export BLIS_NUM_THREADS=1

        python -u "$RATE_SCRIPT" \
            --config "$config" \
            --sample-sizes 1000,2000,4000,8000 \
            --reps "$reps" \
            --x-grid-1d 200 \
            --x-grid-2d 21 \
            --truth-grid-2d 121 \
            --truth-cache-mode off \
            --h0 1.2 \
            --q 0.7071067811865476 \
            --min-h-factor "$min_h_factor" \
            --kappa-pair 2 \
            --kappa-final 2 \
            --beta-effective 2 \
            --selector-metric raw_max \
            --trim-frac 0.05 \
            --penalty-form one_sided \
            --seed 12345 \
            --tag "$tag"

        code=$?
        printf '%s\n' "$code" > "$status_file"
        date --iso-8601=seconds > "${LOG_DIR}/${label}.finished"
        exit "$code"
    ) > "$log" 2>&1
}

declare -a LABELS=()
declare -a PIDS=()

start_job() {
    local label="$1"
    shift

    run_one "$label" "$@" &
    LABELS+=("$label")
    PIDS+=("$!")

    printf '%s\n' "${PIDS[-1]}" > "${LOG_DIR}/${label}.pid"
    echo "Started $label as PID ${PIDS[-1]}"
}

echo "======================================================================"
echo "Priority 5 same-seed end-to-end rerun"
echo "Base: $BASE"
echo "Workers: 8"
echo "Threads per worker: 1"
echo "Truth cache: off"
echo "Logs: $LOG_DIR"
echo "======================================================================"

start_job \
    gg_1d_paper \
    gg_1d configs/gg_1d.yaml 50 20 "${BASE}_gg_1d_paper"

start_job \
    gg_1d_floor81 \
    gg_1d configs/gg_1d.yaml 50 81 "${BASE}_gg_1d_floor81"

start_job \
    mm_1d_paper \
    mm_1d configs/mm_1d.yaml 50 20 "${BASE}_mm_1d_paper"

start_job \
    mm_1d_floor81 \
    mm_1d configs/mm_1d.yaml 50 81 "${BASE}_mm_1d_floor81"

start_job \
    gg_2d_paper \
    gg_2d configs/gg_2d.yaml 20 5 "${BASE}_gg_2d_paper"

start_job \
    gg_2d_floor81 \
    gg_2d configs/gg_2d.yaml 20 9 "${BASE}_gg_2d_floor81"

start_job \
    mm_2d_paper \
    mm_2d configs/mm_2d.yaml 20 5 "${BASE}_mm_2d_paper"

start_job \
    mm_2d_floor81 \
    mm_2d configs/mm_2d.yaml 20 9 "${BASE}_mm_2d_floor81"

overall_status=0

for index in "${!PIDS[@]}"; do
    label="${LABELS[$index]}"
    pid="${PIDS[$index]}"

    if wait "$pid"; then
        echo "$label finished successfully."
    else
        code=$?
        echo "$label failed with wait status $code."
        overall_status=1
    fi
done

echo
echo "================ JOB STATUSES ================"

for label in "${LABELS[@]}"; do
    printf '%-20s ' "$label"

    if [[ -f "${LOG_DIR}/${label}.status" ]]; then
        cat "${LOG_DIR}/${label}.status"
    else
        echo "missing"
        overall_status=1
    fi
done

if [[ "$overall_status" -ne 0 ]]; then
    printf '%s\n' "$overall_status" > "${LOG_DIR}/overall.status"
    echo "One or more experiment jobs failed; slope verification skipped."
    exit "$overall_status"
fi

echo
echo "================ SLOPE VERIFICATION ================"

python - "$MANIFEST" "$BASE" "$LOG_DIR" <<'PY'
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path.cwd()
manifest_path = Path(sys.argv[1])
base = sys.argv[2]
log_dir = Path(sys.argv[3])

tolerance = 1e-8
manifest = pd.read_csv(manifest_path, sep="\t")

comparison_rows: list[dict[str, object]] = []
failures: list[str] = []

for row in manifest.itertuples(index=False):
    new_dir = (
        root
        / "results"
        / "processed"
        / "rate"
        / row.model_id
        / row.tag
    )
    reference_dir = root / row.reference_processed

    new_path = new_dir / "slopes.csv"
    reference_path = reference_dir / "slopes.csv"

    if not new_path.is_file():
        failures.append(f"{row.label}: missing rerun slopes: {new_path}")
        continue

    if not reference_path.is_file():
        failures.append(
            f"{row.label}: missing reference slopes: {reference_path}"
        )
        continue

    new = pd.read_csv(new_path)
    reference = pd.read_csv(reference_path)

    required = {"method", "metric", "slope"}

    if not required.issubset(new.columns):
        failures.append(f"{row.label}: rerun slopes schema is incomplete")
        continue

    if not required.issubset(reference.columns):
        failures.append(f"{row.label}: reference slopes schema is incomplete")
        continue

    merged = reference[
        ["method", "metric", "slope"]
    ].rename(
        columns={"slope": "reference_slope"}
    ).merge(
        new[["method", "metric", "slope"]].rename(
            columns={"slope": "rerun_slope"}
        ),
        on=["method", "metric"],
        how="outer",
        indicator="merge_status",
        validate="one_to_one",
    )

    if not (merged["merge_status"] == "both").all():
        failures.append(
            f"{row.label}: method/metric rows do not match reference"
        )

    merged["absolute_difference"] = (
        merged["rerun_slope"] - merged["reference_slope"]
    ).abs()

    for slope_row in merged.itertuples(index=False):
        status = (
            "PASS"
            if (
                slope_row.merge_status == "both"
                and np.isfinite(slope_row.absolute_difference)
                and slope_row.absolute_difference <= tolerance
            )
            else "FAIL"
        )

        comparison_rows.append(
            {
                "label": row.label,
                "model_id": row.model_id,
                "protocol": row.protocol,
                "method": slope_row.method,
                "metric": slope_row.metric,
                "reference_slope": slope_row.reference_slope,
                "rerun_slope": slope_row.rerun_slope,
                "absolute_difference": slope_row.absolute_difference,
                "tolerance": tolerance,
                "status": status,
            }
        )

        if status != "PASS":
            failures.append(
                f"{row.label}/{slope_row.method}/{slope_row.metric}: "
                f"absolute slope difference "
                f"{slope_row.absolute_difference:.3e}"
            )

comparison = pd.DataFrame(comparison_rows)

if not comparison.empty:
    comparison.to_csv(
        log_dir / "slope_reproduction.csv",
        index=False,
    )

    print(
        comparison.to_string(
            index=False,
            float_format=lambda value: f"{value:.12e}",
        )
    )

print()
print("=" * 90)
print("ORACLE SUP-ERROR SUMMARY")
print("=" * 90)

oracle_sup = comparison[
    (comparison["method"] == "oracle")
    & (comparison["metric"] == "sup_err")
].copy()

if not oracle_sup.empty:
    print(
        oracle_sup[
            [
                "model_id",
                "protocol",
                "reference_slope",
                "rerun_slope",
                "absolute_difference",
                "status",
            ]
        ].to_string(
            index=False,
            float_format=lambda value: f"{value:.12e}",
        )
    )

print()
print("=" * 90)
print("FINAL VERIFICATION STATUS")
print("=" * 90)

if failures:
    print("FAIL")
    for failure in failures:
        print(f"- {failure}")

    (log_dir / "verification.status").write_text(
        "1\n",
        encoding="utf-8",
    )
    raise SystemExit(1)

print("PASS")
print(f"All slope differences are <= {tolerance:.1e}.")

(log_dir / "verification.status").write_text(
    "0\n",
    encoding="utf-8",
)
PY

verification_status=$?

if [[ "$verification_status" -ne 0 ]]; then
    overall_status=1
fi

printf '%s\n' "$overall_status" > "${LOG_DIR}/overall.status"

echo
echo "================ FINAL STATUS ================"
echo "Overall status: $overall_status"
echo "Base: $BASE"
echo "Logs: $LOG_DIR"
echo "Slope comparison: ${LOG_DIR}/slope_reproduction.csv"

exit "$overall_status"
