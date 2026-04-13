# Pre-flight report: GG 1D

- model_id: `gg_1d`
- dim: `1`
- xi0: `[0.0]`
- rate_time: `0.6`
- xs_density_at_xi0: `4.000223e-01`
- xs_density_min_on_eval_xi_grid: `5.413713e-02`
- D_star_min_on_eval_grid: `6.140077e-03`
- D_star_argmin: `[-2.0]`
- truth_convergence_abs_err: `0.000000e+00`
- truth_coarse: `[0.6504065040650405]`
- truth_fine: `[0.6504065040650405]`
- sample_mean_xs: `[-0.0025441947297519737]`
- sample_mean_xu: `[0.2876299080813779]`
- support_accepts_sampling: `True`

## Interpretation

Healthy runs should have: positive and non-tiny marginal density at the chosen conditioning point,
a denominator floor on the actual evaluation box that is comfortably away from zero,
and a small coarse-vs-fine truth error at the probe point.