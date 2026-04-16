# Pre-flight report: MM 2D

- model_id: `mm_2d`
- dim: `2`
- xi0: `[0.8, -0.8]`
- rate_time: `0.6`
- xs_density_at_xi0: `4.672258e-01`
- xs_density_min_on_eval_xi_grid: `4.918337e-09`
- D_star_min_on_eval_grid: `8.261432e-03`
- D_star_argmin: `[-1.5, -1.5]`
- truth_convergence_abs_err: `0.000000e+00`
- truth_coarse: `[-1.2065958856210606, 1.278216723513223]`
- truth_fine: `[-1.2065958856210606, 1.278216723513223]`
- sample_mean_xs: `[-0.09976941573108106, 0.1103510486113835]`
- sample_mean_xu: `[0.3575178236561281, -0.40689163792406446]`
- support_accepts_sampling: `True`

## Interpretation

Healthy runs should have: positive and non-tiny marginal density at the chosen conditioning point,
a denominator floor on the actual evaluation box that is comfortably away from zero,
and a small coarse-vs-fine truth error at the probe point.