# Pre-flight report: MM 1D

- model_id: `mm_1d`
- dim: `1`
- xi0: `[0.8]`
- rate_time: `0.6`
- xs_density_at_xi0: `2.986354e-01`
- xs_density_min_on_eval_xi_grid: `2.532522e-02`
- D_star_min_on_eval_grid: `2.997735e-01`
- D_star_argmin: `[2.0]`
- truth_convergence_abs_err: `0.000000e+00`
- truth_coarse: `[-0.47091290307229977]`
- truth_fine: `[-0.47091290307229977]`
- sample_mean_xs: `[-0.14552015795368722]`
- sample_mean_xu: `[0.5888421321857961]`
- support_accepts_sampling: `True`

## Interpretation

Healthy runs should have: positive and non-tiny marginal density at the chosen conditioning point,
a denominator floor on the actual evaluation box that is comfortably away from zero,
and a small coarse-vs-fine truth error at the probe point.