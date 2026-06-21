# Rebuttal readiness — Day 1

## Repository

- Working branch: `rebuttal-prep`
- Baseline main commit: `f3b2794`
- Submitted snapshot: tag `neurips2026-submission`, commit `85cbebd`
- Separate submitted-version worktree created successfully

## Environment

- Linux Mint
- Python 3.12.3
- Project virtual environment: `.venv`
- Project installed in editable mode

## Preflight status

All configured preflight runs completed successfully:

- GG 1D: passed
- GG 2D: passed
- MM 1D: passed
- MM 2D: passed

For every model:

- `support_accepts_sampling` is `True`
- the estimated denominator floor is positive
- coarse and fine truth calculations agree up to machine precision

## Observation requiring attention

For MM 2D:

- `xs_density_min_on_eval_xi_grid = 4.918337e-09`

This is an extremely low-density point on the evaluation grid. It is not a runtime failure, but it should be reviewed before using MM 2D for new rebuttal experiments involving uniform-density-floor claims.

## Generated-output comparison

Regenerating the preflight reports produced only last-bit floating-point differences and changes from approximately `1e-16` to exact zero. These generated changes were not committed.
