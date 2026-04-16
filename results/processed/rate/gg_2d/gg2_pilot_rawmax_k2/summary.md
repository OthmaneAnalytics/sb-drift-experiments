# Rate experiment summary: GG 2D

- model_id: `gg_2d`
- run_name: `gg2_pilot_rawmax_k2`
- dim: `2`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000]`
- reps: `10`
- t0: `0.6`
- xi0: `[0.0, 0.0]`
- x_grid_per_dim: `21`
- truth_grid_2d: `121`
- h0: `1.2`
- q: `0.7071067811865476`
- min_h_factor: `5.0`
- kappa_pair: `2.0`
- kappa_final: `2.0`
- selector_metric: `raw_max`
- trim_frac: `0.05`
- penalty_form: `one_sided`

## Theoretical benchmarks

- finite-range secant slope: `-0.289357`
- asymptotic slope: `-0.333333`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.149126`
- oracle: slope = `-0.163044`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `5.496392e-01`
- lepski at M=2000: mean sup-grid error = `4.745808e-01`
- lepski at M=4000: mean sup-grid error = `4.469870e-01`
- oracle at M=1000: mean sup-grid error = `4.227624e-01`
- oracle at M=2000: mean sup-grid error = `4.263903e-01`
- oracle at M=4000: mean sup-grid error = `3.372361e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `6.589439e-01`
- lepski at M=2000: mean ISE = `5.536371e-01`
- lepski at M=4000: mean ISE = `4.933762e-01`
- oracle at M=1000: mean ISE = `5.960964e-01`
- oracle at M=2000: mean ISE = `5.496502e-01`
- oracle at M=4000: mean ISE = `3.904330e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.392370`
- M=2000: mean gap = `1.125589`
- M=4000: mean gap = `1.270070`