# Rate experiment summary: MM 2D

- model_id: `mm_2d`
- run_name: `mm2_final_rawmax_k2`
- dim: `2`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `20`
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

- finite-range secant slope: `-0.291150`
- asymptotic slope: `-0.333333`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.209368`
- oracle: slope = `-0.385828`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `1.178173e+00`
- lepski at M=2000: mean sup-grid error = `1.235535e+00`
- lepski at M=4000: mean sup-grid error = `1.043159e+00`
- lepski at M=8000: mean sup-grid error = `7.684646e-01`
- oracle at M=1000: mean sup-grid error = `6.990735e-01`
- oracle at M=2000: mean sup-grid error = `5.585509e-01`
- oracle at M=4000: mean sup-grid error = `4.373553e-01`
- oracle at M=8000: mean sup-grid error = `3.110137e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `1.076232e+00`
- lepski at M=2000: mean ISE = `1.079579e+00`
- lepski at M=4000: mean ISE = `9.383698e-01`
- lepski at M=8000: mean ISE = `7.303674e-01`
- oracle at M=1000: mean ISE = `7.403765e-01`
- oracle at M=2000: mean ISE = `5.951026e-01`
- oracle at M=4000: mean ISE = `4.500393e-01`
- oracle at M=8000: mean ISE = `3.387693e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.660514`
- M=2000: mean gap = `2.315097`
- M=4000: mean gap = `2.677067`
- M=8000: mean gap = `2.711902`