# Rate experiment summary: GG 2D

- model_id: `gg_2d`
- run_name: `gg2_final_rawmax_k2`
- dim: `2`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `20`
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

- finite-range secant slope: `-0.291150`
- asymptotic slope: `-0.333333`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.273725`
- oracle: slope = `-0.243435`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `6.293090e-01`
- lepski at M=2000: mean sup-grid error = `4.560713e-01`
- lepski at M=4000: mean sup-grid error = `4.641652e-01`
- lepski at M=8000: mean sup-grid error = `3.323934e-01`
- oracle at M=1000: mean sup-grid error = `4.591569e-01`
- oracle at M=2000: mean sup-grid error = `3.845699e-01`
- oracle at M=4000: mean sup-grid error = `3.453263e-01`
- oracle at M=8000: mean sup-grid error = `2.711889e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `7.817425e-01`
- lepski at M=2000: mean ISE = `5.667704e-01`
- lepski at M=4000: mean ISE = `5.117813e-01`
- lepski at M=8000: mean ISE = `3.701507e-01`
- oracle at M=1000: mean ISE = `6.458255e-01`
- oracle at M=2000: mean ISE = `5.056751e-01`
- oracle at M=4000: mean ISE = `4.203511e-01`
- oracle at M=8000: mean ISE = `3.544674e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.510755`
- M=2000: mean gap = `1.177866`
- M=4000: mean gap = `1.467309`
- M=8000: mean gap = `1.229083`