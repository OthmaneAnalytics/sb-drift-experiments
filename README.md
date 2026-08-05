# Direct Estimation of Schrödinger Bridge Time-Series Drifts

This is the official experimental repository for:

> Othmane Mazhar and Huyên Pham,  
> **Direct Estimation of Schrödinger Bridge Time-Series Drifts: Finite-Sample, Asymptotic, and Adaptive Guarantees**, 2026.

The paper develops direct nonparametric estimators of time-dependent Schrödinger-bridge drifts and establishes finite-sample, asymptotic, adaptive, and minimax guarantees.

- Paper: https://arxiv.org/abs/2605.05432
- Author ORCID: https://orcid.org/0000-0003-1521-9091

## Contents

The repository contains:

- synthetic Schrödinger-bridge model families;
- deterministic computation of the population drift;
- kernel estimators of the bridge drift;
- finite-sample rate experiments;
- pointwise central-limit-theorem experiments;
- adaptive bandwidth-selection experiments;
- stress tests;
- YAML experiment configurations;
- saved outputs, figures, and tables used in the paper.

## Requirements

- Python 3.10 or later
- NumPy
- SciPy
- pandas
- Matplotlib
- PyYAML
- pytest for development and testing

## Installation

Clone the repository:

```bash
git clone https://github.com/OthmaneAnalytics/sb-drift-experiments.git
cd sb-drift-experiments
```

Create and activate a virtual environment.

On Linux or macOS:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

Install the package and development dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

## Verification

Run the tests:

```bash
pytest -q
```

Run a preflight check for the one-dimensional Gaussian model:

```bash
python scripts/00_preflight.py --config configs/gg_1d.yaml
```

Run preflight checks for all configured model families:

```bash
python scripts/00_preflight.py --config configs/common.yaml --all
```

## Quick smoke tests

A small finite-sample rate experiment can be run with:

```bash
python scripts/01_rate.py \
  --config configs/gg_1d.yaml \
  --sample-sizes 200,400 \
  --reps 5 \
  --tag smoke
```

A small CLT experiment can be run with:

```bash
python scripts/02_clt.py \
  --config configs/gg_1d.yaml \
  --sample-sizes 200,400 \
  --reps 20 \
  --out-tag smoke \
  --no-qq
```

These commands are intended to verify that the installation and experiment pipeline work. They are not the full paper experiments.

## Main experiment entry points

Finite-sample rates and adaptive bandwidth selection:

```bash
python scripts/01_rate.py --config configs/gg_1d.yaml
```

Pointwise Gaussian approximation:

```bash
python scripts/02_clt.py --config configs/gg_1d.yaml
```

Stress-test summary:

```bash
python scripts/03_stress_summary_raw_only.py
```

The available configurations are:

- `configs/gg_1d.yaml`
- `configs/gg_2d.yaml`
- `configs/mm_1d.yaml`
- `configs/mm_1d_stress_wide_strong.yaml`
- `configs/mm_2d.yaml`

Run-specific parameters, random seeds, and output locations are controlled through the command-line arguments and YAML configuration files.

## Repository layout

```text
configs/             experiment configurations
figures/             paper-facing figures
paper_figure_upload/ final figure material
results/             raw and processed experimental outputs
scripts/             experiment drivers and summary scripts
src/sbdrift/         core implementation
tests/               automated tests
```

## Paper-facing artifacts

The principal processed outputs used for the paper are located in:

```text
results/processed/adapt_final/
results/processed/clt_runs/
results/processed/stress/latest/
figures/
```

See `results/README.md` and `FINAL_ARTIFACTS.md` for additional information about saved outputs and final artifacts.

## Computational notes

- All reported experiments were run on CPUs.
- The repository contains both final reported artifacts and exploratory outputs.
- Full experiments may require substantially more time than the smoke tests.
- The supplied random seeds and saved outputs support reproducibility of the paper figures and tables.

## Citation

Please cite the associated paper when using this repository:

```bibtex
@misc{mazhar2026direct,
  title        = {Direct Estimation of Schrödinger Bridge Time-Series Drifts: Finite-Sample, Asymptotic, and Adaptive Guarantees},
  author       = {Mazhar, Othmane and Pham, Huyên},
  year         = {2026},
  eprint       = {2605.05432},
  archivePrefix = {arXiv},
  primaryClass = {math.ST},
  url          = {https://arxiv.org/abs/2605.05432}
}
```

## Licence

This repository is released under the MIT Licence. See `LICENSE` for details.
