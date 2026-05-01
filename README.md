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

    python scripts/00_preflight.py --config configs/gg_1d.yaml
    python scripts/01_rate.py --config configs/gg_1d.yaml ...
    python scripts/02_clt.py --config configs/gg_1d.yaml ...
    python scripts/03_stress_summary_raw_only.py

Run-specific parameters, seeds, and saved outputs are documented in the processed summaries under `results/processed/`.

## Notes

- All reported experiments were CPU-only.
- The repository contains both final reported artifacts and some exploratory outputs.
- The final paper-facing summaries are documented in `results/README.md`.

## Citation

If you use this repository, please cite the associated paper.
