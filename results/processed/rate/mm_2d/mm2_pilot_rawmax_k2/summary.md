# Rate experiment summary: MM 2D

- model_id: `mm_2d`
- run_name: `mm2_pilot_rawmax_k2`
- dim: `2`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000]`
- reps: `10`
- t0: `0.6`
- xi0: `[0.8, -0.8]`
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

- lepski: slope = `-0.526030`
- oracle: slope = `-0.499657`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `1.188429e+00`
- lepski at M=2000: mean sup-grid error = `7.980247e-01`
- lepski at M=4000: mean sup-grid error = `5.731546e-01`
- oracle at M=1000: mean sup-grid error = `7.270558e-01`
- oracle at M=2000: mean sup-grid error = `5.053605e-01`
- oracle at M=4000: mean sup-grid error = `3.637010e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `1.073578e+00`
- lepski at M=2000: mean ISE = `8.037839e-01`
- lepski at M=4000: mean ISE = `5.245439e-01`
- oracle at M=1000: mean ISE = `7.469576e-01`
- oracle at M=2000: mean ISE = `5.114215e-01`
- oracle at M=4000: mean ISE = `3.945471e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.586242`
- M=2000: mean gap = `1.582670`
- M=4000: mean gap = `1.615648`