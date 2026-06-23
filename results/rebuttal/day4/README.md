# Day 4 — Adaptivity provenance audit

## Outcome

The submitted Table 4 was confirmed to combine values from two experiment
protocols.

- Audit baseline commit: dab5b2bdf5a660f94dadd2cf82f52886cffdc0c6
- Frozen protocol commit: 26d82f0
- Frozen adaptive floor: M h^d >= 81
- 1D min_h_factor: 81
- 2D min_h_factor: 9

## Confirmed assembly

The submitted oracle-bandwidth column matches the older saved baseline runs
with implied floors 20 in 1D and 25 in 2D.

The submitted adaptive-bandwidth, boundary-frequency, Cavg, and Cmax columns
match the floor-81 aggregate runs.

Therefore, the submitted table is a hybrid assembly rather than a single
internally consistent experiment table.

See:

- table4_provenance_audit.txt
- table4_provenance_audit.csv
- table4_floor81_consistent.csv
- table4_floor81_consistent.md

## Provenance limitation

Git history preserves the floor-81 aggregate CSV files and the frozen protocol
note, but it does not preserve the original raw run directories,
run_config.json files, or original April launch commands.

The script `scripts/reproduce_floor81_adaptivity.sh` is therefore an explicit
reconstruction from the frozen protocol. It must not be described as the
original launch script.

## Recommended resolution

For consistency with the stated M h^d >= 81 protocol, replace the oracle
bandwidth column with the floor-81 oracle values. The other adaptive diagnostic
columns already use the floor-81 results.

The reconstructed runs will be used to verify the preserved aggregates before
any manuscript or rebuttal wording is finalized.
