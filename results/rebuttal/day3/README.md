# Day 3 — MM1 rate reproduction

## Outcome

The submitted MM1 rate experiment was reproduced successfully on the new
Linux Mint environment.

- Baseline commit: 7eb00b5f82cdb2b7760618d0da3b33e512f973bf
- Model: MM 1D
- Sample sizes: 1000, 2000, 4000, 8000
- Repetitions: 50
- Seed: 12345
- Selector metric: raw_max
- Penalty form: one_sided
- kappa_pair: 2
- kappa_final: 2
- min_h_factor: 20
- Wall-clock runtime: 4 minutes 53.07 seconds
- Peak resident memory: 207784 KB
- Exit status: 0

## Verification

The normalized run configurations match after excluding the machine-specific
config_path and the deliberately different output tag.

All four processed numerical CSV files match the submitted reference within:

- relative tolerance: 1e-12
- absolute tolerance: 1e-14

Observed differences are limited to last-bit floating-point representation.
The mean ISE output matched byte-for-byte.

See:

- mm1_reproduction_check.txt
- logs/mm1_reproduction.log
- ../../processed/rate/mm_1d/day3_mm1_reproduction/

## Observed MM1 behavior

The reproduced adaptive error is non-monotone across the tested finite sample
sizes, including the larger error at M=4000. This behavior exactly matches the
submitted seeded experiment and is therefore not a machine or environment
regression.

## Protocol note

As with GG1, the submitted MM1 run uses:

    h_min = 20 * M^(-1/d)

For d=1, this corresponds to M*h >= 20. The relationship between this stored
configuration and the manuscript's stated stable-regime threshold should be
clarified before the rebuttal is finalized.
