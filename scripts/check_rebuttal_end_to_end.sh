#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

STAMP="e2e_$(date +%Y%m%d_%H%M%S)_$$"
export STAMP

LOG_DIR="results/rebuttal/readiness/e2e_logs/$STAMP"
EDGE_DIR="results/rebuttal/readiness/e2e/$STAMP/edge"
mkdir -p "$LOG_DIR"

run_rate() {
    local config="$1"
    local label="$2"

    echo "=== Rate: $label ==="

    python scripts/01_rate.py \
        --config "$config" \
        --sample-sizes 64,128 \
        --reps 2 \
        --x-grid-1d 11 \
        --x-grid-2d 5 \
        --truth-grid-2d 31 \
        --h0 1.2 \
        --q 0.7071067811865476 \
        --min-h-factor 5 \
        --kappa-pair 2 \
        --kappa-final 2 \
        --selector-metric raw_max \
        --penalty-form one_sided \
        --seed 12345 \
        --tag "${STAMP}_${label}" \
        > "$LOG_DIR/rate_${label}.log" 2>&1
}

run_rate configs/gg_1d.yaml gg1
run_rate configs/gg_2d.yaml gg2
run_rate configs/mm_1d.yaml mm1
run_rate configs/mm_2d.yaml mm2
run_rate configs/mm_1d_stress_wide_strong.yaml mm1_stress

echo "=== CLT: GG1 ==="

python scripts/02_clt.py \
    --config configs/gg_1d.yaml \
    --sample-sizes 128 \
    --reps 8 \
    --alpha 0.22 \
    --c 1.0 \
    --seed 22345 \
    --out-tag "${STAMP}_gg1_clt" \
    --no-qq \
    --progress-every 8 \
    > "$LOG_DIR/clt_gg1.log" 2>&1

echo "=== CLT: MM1 ==="

python scripts/02_clt.py \
    --config configs/mm_1d.yaml \
    --sample-sizes 128 \
    --reps 8 \
    --alpha 0.28 \
    --c 1.0 \
    --seed 32345 \
    --out-tag "${STAMP}_mm1_clt" \
    --no-qq \
    --progress-every 8 \
    > "$LOG_DIR/clt_mm1.log" 2>&1

echo "=== Terminal edge ==="

python scripts/03_edge.py \
    --gg-config configs/gg_1d.yaml \
    --mm-config configs/mm_1d.yaml \
    --M 128 \
    --reps 2 \
    --seed 42345 \
    --method lepski \
    --kappa 2 \
    --h0 1.2 \
    --q 0.7071067811865476 \
    --min-h-factor 5 \
    --x-grid-1d 11 \
    --truth-grid-2d 31 \
    --outdir "$EDGE_DIR" \
    > "$LOG_DIR/edge.log" 2>&1

python - <<'PY'
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path.cwd()
STAMP = os.environ["STAMP"]

RATE_CASES = [
    ("configs/gg_1d.yaml", "gg1"),
    ("configs/gg_2d.yaml", "gg2"),
    ("configs/mm_1d.yaml", "mm1"),
    ("configs/mm_2d.yaml", "mm2"),
    ("configs/mm_1d_stress_wide_strong.yaml", "mm1_stress"),
]

CLT_CASES = [
    ("configs/gg_1d.yaml", f"{STAMP}_gg1_clt"),
    ("configs/mm_1d.yaml", f"{STAMP}_mm1_clt"),
]

checks: list[str] = []
cleanup_paths: list[Path] = []
cleanup_files: list[Path] = []


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def require_file(path: Path) -> None:
    require(path.is_file(), f"Missing file: {path}")
    require(path.stat().st_size > 0, f"Empty file: {path}")


def require_finite(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    values = frame[columns].to_numpy(dtype=float)
    require(
        np.all(np.isfinite(values)),
        f"{label} contains NaN or infinity in {columns}",
    )


def read_config(path: str) -> dict:
    return yaml.safe_load((ROOT / path).read_text(encoding="utf-8"))


# Rate and stress pipelines
for config_path, label in RATE_CASES:
    cfg = read_config(config_path)
    model_id = str(cfg.get("model_id", Path(config_path).stem))
    tag = f"{STAMP}_{label}"

    raw_dir = ROOT / "results/raw/rate" / model_id / tag
    proc_dir = ROOT / "results/processed/rate" / model_id / tag

    selected_path = raw_dir / "selected_metrics.csv"
    per_h_path = raw_dir / "per_h_metrics.csv"

    required_processed = [
        proc_dir / "mean_sup_err.csv",
        proc_dir / "mean_ise.csv",
        proc_dir / "slopes.csv",
        proc_dir / "bandwidth_diagnostics.csv",
        proc_dir / "run_config.json",
        proc_dir / "summary.md",
    ]

    require_file(selected_path)
    require_file(per_h_path)

    for path in required_processed:
        require_file(path)

    selected = pd.read_csv(selected_path)
    per_h = pd.read_csv(per_h_path)

    require(len(selected) == 8, f"{label}: expected 8 selected rows")
    require(set(selected["M"]) == {64, 128}, f"{label}: wrong M values")
    require(
        set(selected["method"]) == {"oracle", "lepski"},
        f"{label}: missing method",
    )

    require_finite(
        selected,
        [
            "selected_h",
            "sup_err",
            "ise",
            "gap_to_oracle",
            "f_hat",
            "min_D_hat",
        ],
        label,
    )
    require_finite(
        per_h,
        ["h", "sup_err", "ise"],
        f"{label} per-h",
    )

    require(
        np.all(selected["selected_h"] > 0),
        f"{label}: nonpositive bandwidth",
    )
    require(
        np.all(selected["f_hat"] > 0),
        f"{label}: nonpositive f_hat",
    )
    require(
        np.all(selected["min_D_hat"] > 0),
        f"{label}: nonpositive D_hat",
    )

    oracle = selected[selected["method"] == "oracle"]
    require(
        np.allclose(oracle["gap_to_oracle"], 1.0),
        f"{label}: oracle ratio is not one",
    )

    for row in selected.itertuples(index=False):
        candidates = per_h.loc[
            (per_h["M"] == row.M) & (per_h["rep"] == row.rep),
            "h",
        ].to_numpy(dtype=float)

        require(
            np.any(
                np.isclose(
                    candidates,
                    float(row.selected_h),
                    rtol=0,
                    atol=1e-12,
                )
            ),
            f"{label}: selected bandwidth is absent from candidate grid",
        )

    for metric in ["sup", "ise"]:
        figure = (
            ROOT
            / "results/figures"
            / f"rate_{metric}_{model_id}_{tag}.pdf"
        )
        require_file(figure)
        cleanup_files.append(figure)

    cleanup_paths.extend([raw_dir, proc_dir])
    checks.append(f"PASS: rate pipeline {label}")


# CLT pipelines
for config_path, tag in CLT_CASES:
    cfg = read_config(config_path)
    model_id = str(cfg.get("model_id", Path(config_path).stem))

    raw_dir = ROOT / "results/raw/clt_runs" / tag / model_id
    proc_dir = ROOT / "results/processed/clt_runs" / tag / model_id
    fig_dir = ROOT / "results/figures/clt_runs" / tag

    raw_path = raw_dir / "pointwise_clt.csv"
    summary_path = proc_dir / "summary.csv"

    require_file(raw_path)
    require_file(summary_path)
    require_file(proc_dir / "run_config.json")
    require_file(proc_dir / "summary.md")

    raw = pd.read_csv(raw_path)
    summary = pd.read_csv(summary_path)

    require(len(raw) == 8, f"{model_id}: expected 8 CLT rows")
    require(len(summary) == 1, f"{model_id}: expected one summary row")
    require(set(raw["M"]) == {128}, f"{model_id}: wrong CLT M")

    require_finite(
        raw,
        [
            "h",
            "a_hat",
            "a_star",
            "f_hat",
            "D_hat",
            "sigma_hat",
            "Z",
            "covered_95",
        ],
        f"{model_id} CLT",
    )

    require(
        np.all(raw["h"] > 0),
        f"{model_id}: nonpositive CLT bandwidth",
    )
    require(
        np.all(raw["f_hat"] > 0),
        f"{model_id}: nonpositive CLT f_hat",
    )
    require(
        np.all(raw["D_hat"] > 0),
        f"{model_id}: nonpositive CLT D_hat",
    )
    require(
        np.all(raw["sigma_hat"] > 0),
        f"{model_id}: nonpositive variance estimate",
    )
    require(
        set(raw["covered_95"]).issubset({0, 1}),
        f"{model_id}: invalid coverage indicator",
    )

    cleanup_paths.extend(
        [
            ROOT / "results/raw/clt_runs" / tag,
            ROOT / "results/processed/clt_runs" / tag,
            fig_dir,
        ]
    )
    checks.append(f"PASS: CLT pipeline {model_id}")


# Edge pipeline
edge_dir = ROOT / "results/rebuttal/readiness/e2e" / STAMP / "edge"
edge_json = edge_dir / "plot_data/edge_plot_data.json"

require_file(edge_json)
require_file(edge_dir / "edge_summary.md")
require_file(edge_dir / "edge_error.pdf")
require_file(edge_dir / "edge_error_sup.pdf")
require_file(edge_dir / "edge_rescaled.pdf")
require_file(edge_dir / "edge_rescaled_sup.pdf")

payload = json.loads(edge_json.read_text(encoding="utf-8"))

require(set(payload) == {"GG1", "MM1"}, "Edge output has wrong families")

for family in ["GG1", "MM1"]:
    result = payload[family]

    require(result["reps"] == 2, f"{family}: wrong edge repetition count")
    require(result["M"] == 128, f"{family}: wrong edge sample size")
    require(result["method"] == "lepski", f"{family}: wrong edge method")
    require(result["mean_h"] > 0, f"{family}: nonpositive edge bandwidth")

    for key in ["times", "mean_err", "std_err", "mean_rescaled"]:
        values = np.asarray(result[key], dtype=float)
        require(values.shape == (5,), f"{family}: malformed {key}")
        require(
            np.all(np.isfinite(values)),
            f"{family}: non-finite values in {key}",
        )

checks.append("PASS: terminal-edge pipeline")
cleanup_paths.append(ROOT / "results/rebuttal/readiness/e2e" / STAMP)


# Write the single readiness report.
report = ROOT / "results/rebuttal/readiness/end_to_end.md"
report.parent.mkdir(parents=True, exist_ok=True)

lines = [
    "# End-to-end execution validation",
    "",
    "Overall result: **PASS**",
    "",
    "Fresh small runs completed successfully through every experimental "
    "pipeline required for likely rebuttal questions.",
    "",
]

lines.extend(f"- {check}" for check in checks)

lines.extend(
    [
        "",
        "Validated properties:",
        "",
        "- expected files were generated;",
        "- expected repetition and sample-size counts were present;",
        "- central numerical outputs were finite;",
        "- density and denominator diagnostics were positive;",
        "- selected bandwidths belonged to the generated candidate grids;",
        "- CLT variance estimates and coverage indicators were valid;",
        "- terminal-edge raw and rescaled errors were finite.",
        "",
        "Temporary numerical run outputs were deleted after validation.",
        "",
        f"Validation stamp: `{STAMP}`",
    ]
)

report.write_text("\n".join(lines) + "\n", encoding="utf-8")

# Delete temporary numerical outputs only after every check has passed.
for path in cleanup_paths:
    shutil.rmtree(path, ignore_errors=True)

for path in cleanup_files:
    path.unlink(missing_ok=True)

log_dir = ROOT / "results/rebuttal/readiness/e2e_logs" / STAMP
shutil.rmtree(log_dir, ignore_errors=True)

print(report.read_text(encoding="utf-8"))
PY

echo "Overall: PASS"
