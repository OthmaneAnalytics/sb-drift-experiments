# Rate experiment summary: GG 1D

- model_id: `gg_1d`
- run_name: `gg1_final_rawmax_k2`
- dim: `1`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `50`
- t0: `0.6`
- xi0: `[0.0]`
- x_grid_per_dim: `200`
- truth_grid_2d: `121`
- h0: `1.2`
- q: `0.7071067811865476`
- min_h_factor: `20.0`
- kappa_pair: `2.0`
- kappa_final: `2.0`
- selector_metric: `raw_max`
- trim_frac: `0.05`
- penalty_form: `one_sided`

## Theoretical benchmarks

- finite-range secant slope: `-0.349379`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.261710`
- oracle: slope = `-0.348899`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `2.440648e-01`
- lepski at M=2000: mean sup-grid error = `2.223155e-01`
- lepski at M=4000: mean sup-grid error = `1.963797e-01`
- lepski at M=8000: mean sup-grid error = `1.389487e-01`
- oracle at M=1000: mean sup-grid error = `1.620541e-01`
- oracle at M=2000: mean sup-grid error = `1.433734e-01`
- oracle at M=4000: mean sup-grid error = `1.139622e-01`
- oracle at M=8000: mean sup-grid error = `7.812654e-02`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `2.379710e-01`
- lepski at M=2000: mean ISE = `2.078020e-01`
- lepski at M=4000: mean ISE = `1.932629e-01`
- lepski at M=8000: mean ISE = `1.320401e-01`
- oracle at M=1000: mean ISE = `1.579346e-01`
- oracle at M=2000: mean ISE = `1.451162e-01`
- oracle at M=4000: mean ISE = `1.141242e-01`
- oracle at M=8000: mean ISE = `7.852177e-02`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.643389`
- M=2000: mean gap = `1.607961`
- M=4000: mean gap = `1.698153`
- M=8000: mean gap = `2.070325`