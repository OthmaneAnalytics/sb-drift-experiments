# Rate experiment summary: MM 2D

- model_id: `mm_2d`
- run_name: `mm_2d_adapt_dense_baseline_20260417_164626`
- dim: `2`
- beta_effective: `2.0`
- sample_sizes: `[1000, 1500, 2000, 3000, 4000, 6000, 8000]`
- reps: `20`
- t0: `0.6`
- xi0: `[0.8, -0.8]`
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

- lepski: slope = `-0.060662`
- oracle: slope = `-0.295673`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `2.041258e+00`
- lepski at M=1500: mean sup-grid error = `3.297987e+00`
- lepski at M=2000: mean sup-grid error = `2.544084e+00`
- lepski at M=3000: mean sup-grid error = `2.105695e+00`
- lepski at M=4000: mean sup-grid error = `1.537270e+00`
- lepski at M=6000: mean sup-grid error = `2.952984e+00`
- lepski at M=8000: mean sup-grid error = `2.072866e+00`
- oracle at M=1000: mean sup-grid error = `6.902773e-01`
- oracle at M=1500: mean sup-grid error = `5.407493e-01`
- oracle at M=2000: mean sup-grid error = `5.521955e-01`
- oracle at M=3000: mean sup-grid error = `5.262286e-01`
- oracle at M=4000: mean sup-grid error = `4.537312e-01`
- oracle at M=6000: mean sup-grid error = `3.703996e-01`
- oracle at M=8000: mean sup-grid error = `3.659636e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `3.285818e+00`
- lepski at M=1500: mean ISE = `6.329331e+00`
- lepski at M=2000: mean ISE = `4.518689e+00`
- lepski at M=3000: mean ISE = `3.686399e+00`
- lepski at M=4000: mean ISE = `2.517309e+00`
- lepski at M=6000: mean ISE = `5.582343e+00`
- lepski at M=8000: mean ISE = `3.759098e+00`
- oracle at M=1000: mean ISE = `7.450455e-01`
- oracle at M=1500: mean ISE = `6.341748e-01`
- oracle at M=2000: mean ISE = `5.075218e-01`
- oracle at M=3000: mean ISE = `5.483082e-01`
- oracle at M=4000: mean ISE = `4.840145e-01`
- oracle at M=6000: mean ISE = `3.644649e-01`
- oracle at M=8000: mean ISE = `3.809670e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `2.912289`
- M=1500: mean gap = `6.768727`
- M=2000: mean gap = `5.348724`
- M=3000: mean gap = `4.846168`
- M=4000: mean gap = `2.776397`
- M=6000: mean gap = `10.179754`
- M=8000: mean gap = `7.016285`