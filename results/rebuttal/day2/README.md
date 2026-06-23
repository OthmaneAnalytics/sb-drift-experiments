# Day 2 — GG1 rate reproduction

## Outcome

The submitted GG1 rate experiment was reproduced successfully on the new
Linux Mint environment.

- Baseline commit: d1906322c1650813207023acfacd59fee4cb92d9
- Model: GG 1D
- Sample sizes: 1000, 2000, 4000, 8000
- Repetitions: 50
- Seed: 12345
- Selector metric: raw_max
- Penalty form: one_sided
- kappa_pair: 2
- kappa_final: 2
- min_h_factor: 20
- Wall-clock runtime: 5 minutes 47.92 seconds
- Peak resident memory: 208296 KB
- Exit status: 0

## Verification

The normalized run configurations match after excluding the machine-specific
config_path and the deliberately different output tag.

The generated summary differs from the submitted summary only in run_name.

All numerical CSV outputs pass comparison with:

- relative tolerance: 1e-12
- absolute tolerance: 1e-14

Observed differences are limited to last-bit floating-point representation.

See:

- gg1_reproduction_check.txt
- logs/gg1_reproduction.log
- ../../processed/rate/gg_1d/day2_gg1_reproduction/

## Protocol note requiring follow-up

The submitted GG1 run uses:

    h_min = 20 * M^(-1/d)

For d=1, this corresponds to M*h >= 20. This should be checked against the
stable-regime threshold described in the manuscript before preparing the
rebuttal response. No configuration or manuscript text has been changed yet.
