# Pre-flight report: MM 2D

- model_id: `mm_2d`
- dim: `2`
- xi0: `[0.8, -0.8]`
- rate_time: `0.6`
- xs_density_at_xi0: `3.115879e-01`
- xs_density_min_on_eval_xi_grid: `4.394140e-19`
- D_star_min_on_eval_grid: `2.452504e-04`
- D_star_argmin: `[-2.0, -2.0]`
- truth_convergence_abs_err: `2.220446e-16`
- truth_coarse: `[-1.2065958856210606, 1.278216723513223]`
- truth_fine: `[-1.2065958856210606, 1.2782167235132227]`
- sample_mean_xs: `[-0.12354823876469599, 0.13280716753496055]`
- sample_mean_xu: `[0.49415532979114823, -0.575717007718261]`
- support_accepts_sampling: `True`

## Interpretation

Healthy runs should have: positive and non-tiny marginal density at the chosen conditioning point,
a denominator floor on the actual evaluation box that is comfortably away from zero,
and a small coarse-vs-fine truth error at the probe point.