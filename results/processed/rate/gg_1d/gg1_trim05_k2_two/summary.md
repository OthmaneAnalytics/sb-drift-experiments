# Rate experiment summary: GG 1D

- model_id: `gg_1d`
- run_name: `gg1_trim05_k2_two`
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
- penalty_form: `two_sided`

## Theoretical benchmarks

- finite-range secant slope: `-0.347228`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.009076`
- oracle: slope = `-0.145385`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `2.954204e-01`
- lepski at M=2000: mean sup-grid error = `3.494212e-01`
- lepski at M=4000: mean sup-grid error = `2.917267e-01`
- oracle at M=1000: mean sup-grid error = `1.438550e-01`
- oracle at M=2000: mean sup-grid error = `1.327820e-01`
- oracle at M=4000: mean sup-grid error = `1.175965e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `2.879628e-01`
- lepski at M=2000: mean ISE = `3.126289e-01`
- lepski at M=4000: mean ISE = `2.859405e-01`
- oracle at M=1000: mean ISE = `1.309038e-01`
- oracle at M=2000: mean ISE = `1.309163e-01`
- oracle at M=4000: mean ISE = `1.124274e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `2.731617`
- M=2000: mean gap = `3.398623`
- M=4000: mean gap = `3.266171`