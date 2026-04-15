# Rate experiment summary: GG 1D

- model_id: `gg_1d`
- run_name: `gg1_trim05_k2`
- dim: `1`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000]`
- reps: `20`
- t0: `0.6`
- xi0: `[0.0]`
- x_grid_per_dim: `200`
- truth_grid_2d: `121`
- h0: `1.2`
- q: `0.7071067811865476`
- min_h_factor: `20.0`
- kappa: `2.0`
- selector_metric: `trimmed_max`
- trim_frac: `0.05`

## Theoretical benchmarks

- finite-range secant slope: `-0.347228`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.048979`
- oracle: slope = `-0.145385`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `2.204296e-01`
- lepski at M=2000: mean sup-grid error = `2.871961e-01`
- lepski at M=4000: mean sup-grid error = `2.059594e-01`
- oracle at M=1000: mean sup-grid error = `1.438550e-01`
- oracle at M=2000: mean sup-grid error = `1.327820e-01`
- oracle at M=4000: mean sup-grid error = `1.175965e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `2.142469e-01`
- lepski at M=2000: mean ISE = `2.465367e-01`
- lepski at M=4000: mean ISE = `1.902705e-01`
- oracle at M=1000: mean ISE = `1.309038e-01`
- oracle at M=2000: mean ISE = `1.309163e-01`
- oracle at M=4000: mean ISE = `1.124274e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.787431`
- M=2000: mean gap = `2.246871`
- M=4000: mean gap = `2.022163`