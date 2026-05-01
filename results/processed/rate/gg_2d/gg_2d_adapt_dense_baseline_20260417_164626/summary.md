# Rate experiment summary: GG 2D

- model_id: `gg_2d`
- run_name: `gg_2d_adapt_dense_baseline_20260417_164626`
- dim: `2`
- beta_effective: `2.0`
- sample_sizes: `[1000, 1500, 2000, 3000, 4000, 6000, 8000]`
- reps: `20`
- t0: `0.6`
- xi0: `[0.0, 0.0]`
- x_grid_per_dim: `21`
- truth_grid_2d: `121`
- h0: `1.2`
- q: `0.7071067811865476`
- min_h_factor: `1.0`
- kappa_pair: `2.0`
- kappa_final: `2.0`
- selector_metric: `raw_max`
- trim_frac: `0.05`
- penalty_form: `one_sided`

## Theoretical benchmarks

- finite-range secant slope: `-0.291150`
- asymptotic slope: `-0.333333`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.149103`
- oracle: slope = `-0.225770`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `6.187143e-01`
- lepski at M=1500: mean sup-grid error = `4.738997e-01`
- lepski at M=2000: mean sup-grid error = `5.412980e-01`
- lepski at M=3000: mean sup-grid error = `4.232790e-01`
- lepski at M=4000: mean sup-grid error = `4.863737e-01`
- lepski at M=6000: mean sup-grid error = `4.305921e-01`
- lepski at M=8000: mean sup-grid error = `4.277971e-01`
- oracle at M=1000: mean sup-grid error = `4.591314e-01`
- oracle at M=1500: mean sup-grid error = `4.192082e-01`
- oracle at M=2000: mean sup-grid error = `4.365679e-01`
- oracle at M=3000: mean sup-grid error = `3.550407e-01`
- oracle at M=4000: mean sup-grid error = `3.719877e-01`
- oracle at M=6000: mean sup-grid error = `3.041850e-01`
- oracle at M=8000: mean sup-grid error = `2.906416e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `7.619083e-01`
- lepski at M=1500: mean ISE = `5.942738e-01`
- lepski at M=2000: mean ISE = `6.398771e-01`
- lepski at M=3000: mean ISE = `5.011999e-01`
- lepski at M=4000: mean ISE = `6.038539e-01`
- lepski at M=6000: mean ISE = `4.825178e-01`
- lepski at M=8000: mean ISE = `5.269277e-01`
- oracle at M=1000: mean ISE = `6.358747e-01`
- oracle at M=1500: mean ISE = `5.942965e-01`
- oracle at M=2000: mean ISE = `5.287420e-01`
- oracle at M=3000: mean ISE = `4.724986e-01`
- oracle at M=4000: mean ISE = `4.577515e-01`
- oracle at M=6000: mean ISE = `3.904171e-01`
- oracle at M=8000: mean ISE = `3.684036e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.488388`
- M=1500: mean gap = `1.105093`
- M=2000: mean gap = `1.297492`
- M=3000: mean gap = `1.212026`
- M=4000: mean gap = `1.336898`
- M=6000: mean gap = `1.436428`
- M=8000: mean gap = `1.466056`