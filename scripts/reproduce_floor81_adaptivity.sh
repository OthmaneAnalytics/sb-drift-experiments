#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-}"
TARGET="${2:-}"

if [[ "$MODE" != "main" && "$MODE" != "dense" ]]; then
    echo "Usage: $0 {main|dense} {gg_1d|gg_2d|mm_1d|mm_2d|all}" >&2
    exit 2
fi

case "$TARGET" in
    gg_1d|gg_2d|mm_1d|mm_2d|all)
        ;;
    *)
        echo "Usage: $0 {main|dense} {gg_1d|gg_2d|mm_1d|mm_2d|all}" >&2
        exit 2
        ;;
esac

PYTHON="${PYTHON:-.venv/bin/python}"
RUN_ID="${RUN_ID:-rebuttal_floor81_${MODE}_$(date +%Y%m%d_%H%M%S)}"

if [[ ! -x "$PYTHON" ]]; then
    echo "Python executable not found: $PYTHON" >&2
    exit 1
fi

if [[ "$MODE" == "main" ]]; then
    SAMPLE_SIZES="1000,2000,4000,8000"
else
    SAMPLE_SIZES="1000,1500,2000,3000,4000,6000,8000"
fi

mkdir -p results/rebuttal/day4/logs

run_model() {
    local model="$1"
    local reps
    local min_h_factor
    local tag
    local log
    local processed_dir
    local raw_dir

    case "$model" in
        gg_1d|mm_1d)
            reps=50
            min_h_factor=81
            ;;
        gg_2d|mm_2d)
            reps=20
            min_h_factor=9
            ;;
        *)
            echo "Unsupported model: $model" >&2
            exit 2
            ;;
    esac

    tag="${model}_${RUN_ID}"
    log="results/rebuttal/day4/logs/${tag}.log"
    processed_dir="results/processed/rate/${model}/${tag}"
    raw_dir="results/raw/rate/${model}/${tag}"

    if [[ -e "$processed_dir" || -e "$raw_dir" ]]; then
        echo "Refusing to overwrite existing run: $tag" >&2
        exit 1
    fi

    echo "============================================================"
    echo "Model:          $model"
    echo "Mode:           $MODE"
    echo "Tag:            $tag"
    echo "Sample sizes:   $SAMPLE_SIZES"
    echo "Repetitions:    $reps"
    echo "min_h_factor:   $min_h_factor"
    echo "Implied floor:  M*h^d >= 81"
    echo "Log:            $log"
    echo "============================================================"

    /usr/bin/time -v "$PYTHON" scripts/01_rate.py \
        --config "configs/${model}.yaml" \
        --sample-sizes "$SAMPLE_SIZES" \
        --reps "$reps" \
        --x-grid-1d 200 \
        --truth-grid-2d 121 \
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
        --tag "$tag" \
        2>&1 | tee "$log"
}

if [[ "$TARGET" == "all" ]]; then
    for model in gg_1d gg_2d mm_1d mm_2d; do
        run_model "$model"
    done
else
    run_model "$TARGET"
fi
