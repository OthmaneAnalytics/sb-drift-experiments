# Schrödinger-Bridge Drift Estimation Experiments

This repository accompanies the paper on direct nonparametric estimation of Schrödinger-bridge drifts from a closed-form conditional-ratio representation.

It contains:
- synthetic model definitions,
- deterministic truth computation,
- kernel drift estimators,
- rate / CLT / adaptivity / stress-test drivers,
- YAML configurations,
- saved artifacts used for the paper figures and tables.

## Repository layout

- `src/sbdrift/` — core implementation
- `configs/` — experiment configurations
- `scripts/` — experiment drivers and summary scripts
- `figures/` — paper-facing figures
- `results/` — saved experiment outputs and summaries

## Canonical paper-facing artifacts

The final paper primarily uses:
- `results/processed/adapt_final/`
- `results/processed/clt_runs/`
- `results/processed/stress/latest/`
- `figures/`

See `results/README.md` for a guide to saved outputs.

## Reproducing the main experiments

Typical entry points are:

```bash
python scripts/00_preflight.py --config configs/gg_1d.yaml
python scripts/01_rate.py --config configs/gg_1d.yaml ...
python scripts/02_clt.py --config configs/gg_1d.yaml ...
python scripts/03_stress_summary_raw_only.py

### 2) Add a results guide

```bash
mkdir -p results

cat > results/README.md <<'MD'
# Results directory guide

This directory contains both final paper-facing artifacts and some exploratory outputs.

## Canonical paper-facing outputs

- `processed/adapt_final/` — final adaptivity summaries and tables
- `processed/clt_runs/` — final CLT summaries and diagnostics
- `processed/stress/latest/` — canonical bounded-support stress summary
- `../figures/` — paper-facing figures

## Raw experiment outputs

- `raw/rate/`
- `raw/clt_runs/`

These contain repetition-level saved outputs used to reconstruct summaries.

## Exploratory outputs

Some subdirectories correspond to timing probes, pilot runs, calibration runs, or superseded experiments. They are kept for transparency but are not part of the paper’s quantitative claims.

## Recommendation

For manuscript verification, start with:
- `processed/adapt_final/`
- `processed/clt_runs/`
- `processed/stress/latest/`
