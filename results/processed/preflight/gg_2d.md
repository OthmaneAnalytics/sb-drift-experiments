# Pre-flight report: GG 2D

- model_id: `gg_2d`
- dim: `2`
- xi0: `[0.0, 0.0]`
- rate_time: `0.6`
- xs_density_at_xi0: `1.785645e-01`
- xs_density_min_on_eval_xi_grid: `1.420651e-02`
- D_star_min_on_eval_grid: `2.285052e-03`
- D_star_argmin: `[-1.5, 1.5]`
- truth_convergence_abs_err: `1.110223e-16`
- truth_coarse: `[0.5463595323532792, -0.4525986804028248]`
- truth_fine: `[0.5463595323532792, -0.4525986804028247]`
- sample_mean_xs: `[0.00617392426084024, -0.0018881615496286032]`
- sample_mean_xu: `[0.2731095035662283, -0.17655299272100658]`
- support_accepts_sampling: `True`

## Interpretation

Healthy runs should have: positive and non-tiny marginal density at the chosen conditioning point,
a denominator floor on the actual evaluation box that is comfortably away from zero,
and a small coarse-vs-fine truth error at the probe point.