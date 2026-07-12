#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

AUDIT_SCRIPT="scripts/rebuttal/priority4_adaptive_selector_audit.py"

if [[ ! -s "$AUDIT_SCRIPT" ]]; then
    echo "ERROR: missing or empty audit script: $AUDIT_SCRIPT" >&2
    exit 1
fi

if ! python -m py_compile "$AUDIT_SCRIPT"; then
    echo "ERROR: Priority 4 audit script does not compile." >&2
    exit 1
fi

P4_BASE="${P4_BASE:-priority4_parallel_$(date +%Y%m%d_%H%M%S)}"
LOG_ROOT="results/rebuttal/priority4_parallel_logs"
P4_LOG_DIR="${LOG_ROOT}/${P4_BASE}"

mkdir -p "$P4_LOG_DIR"
printf '%s\n' "$P4_BASE" > "${LOG_ROOT}/LATEST"

PHYSICAL_CORES="$(
    lscpu -p=CORE 2>/dev/null \
        | grep -v '^#' \
        | sort -u \
        | wc -l
)"

if [[ -z "$PHYSICAL_CORES" || "$PHYSICAL_CORES" -lt 1 ]]; then
    PHYSICAL_CORES="$(nproc)"
fi

THREADS_PER_JOB=$(( PHYSICAL_CORES / 4 ))

if [[ "$THREADS_PER_JOB" -lt 1 ]]; then
    THREADS_PER_JOB=1
fi

# Avoid excessive nested parallelism.
if [[ "$THREADS_PER_JOB" -gt 2 ]]; then
    THREADS_PER_JOB=2
fi

echo "============================================================"
echo "Priority 4 parallel run"
echo "Base tag:          $P4_BASE"
echo "Physical cores:    $PHYSICAL_CORES"
echo "Threads per job:   $THREADS_PER_JOB"
echo "Logs:              $P4_LOG_DIR"
echo "============================================================"

run_model() {
    local model="$1"
    local log="${P4_LOG_DIR}/${model}.log"
    local status_file="${P4_LOG_DIR}/${model}.status"
    local start_file="${P4_LOG_DIR}/${model}.started"
    local finish_file="${P4_LOG_DIR}/${model}.finished"

    date --iso-8601=seconds > "$start_file"

    (
        export OMP_NUM_THREADS="$THREADS_PER_JOB"
        export OPENBLAS_NUM_THREADS="$THREADS_PER_JOB"
        export MKL_NUM_THREADS="$THREADS_PER_JOB"
        export NUMEXPR_NUM_THREADS="$THREADS_PER_JOB"
        export VECLIB_MAXIMUM_THREADS="$THREADS_PER_JOB"

        python "$AUDIT_SCRIPT" \
            --configs "configs/${model}.yaml" \
            --sample-sizes 1000,2000,4000,8000 \
            --reps-1d 50 \
            --reps-2d 20 \
            --nominal-volumes 49,81,121 \
            --metrics raw_max,trimmed_max,ise \
            --kappa-pairs 1.5:1.5,1.5:2,2:1.5,2:2,2.5:2,2:2.5,2.5:2.5 \
            --penalties one_sided,two_sided \
            --trim-frac 0.05 \
            --h0 1.2 \
            --q 0.7071067811865476 \
            --x-grid-1d 200 \
            --x-grid-2d 21 \
            --truth-grid-2d 121 \
            --truth-cache-mode use \
            --seed 12345 \
            --tag "${P4_BASE}_${model}" \
            --progress-every 10 \
            --compare-priority1
    ) > "$log" 2>&1

    local status=$?
    printf '%s\n' "$status" > "$status_file"
    date --iso-8601=seconds > "$finish_file"

    return "$status"
}

run_model gg_1d &
PID_GG1=$!

run_model gg_2d &
PID_GG2=$!

run_model mm_1d &
PID_MM1=$!

run_model mm_2d &
PID_MM2=$!

printf '%s\n' "$PID_GG1" > "${P4_LOG_DIR}/gg_1d.pid"
printf '%s\n' "$PID_GG2" > "${P4_LOG_DIR}/gg_2d.pid"
printf '%s\n' "$PID_MM1" > "${P4_LOG_DIR}/mm_1d.pid"
printf '%s\n' "$PID_MM2" > "${P4_LOG_DIR}/mm_2d.pid"

echo "Started:"
echo "  gg_1d PID $PID_GG1"
echo "  gg_2d PID $PID_GG2"
echo "  mm_1d PID $PID_MM1"
echo "  mm_2d PID $PID_MM2"

overall_status=0

wait "$PID_GG1" || overall_status=1
wait "$PID_GG2" || overall_status=1
wait "$PID_MM1" || overall_status=1
wait "$PID_MM2" || overall_status=1

echo
echo "================ FINAL JOB STATUSES ================"

for model in gg_1d gg_2d mm_1d mm_2d; do
    status_file="${P4_LOG_DIR}/${model}.status"

    if [[ -f "$status_file" ]]; then
        status="$(cat "$status_file")"
    else
        status="missing"
        overall_status=1
    fi

    echo "$model: $status"
done

printf '%s\n' "$overall_status" > "${P4_LOG_DIR}/overall.status"

echo
echo "Overall status: $overall_status"
echo "Base tag: $P4_BASE"
echo "Logs: $P4_LOG_DIR"

exit "$overall_status"
