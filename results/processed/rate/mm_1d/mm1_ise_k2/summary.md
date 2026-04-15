# Rate experiment summary: MM 1D

- model_id: `mm_1d`
- run_name: `mm1_ise_k2`
- dim: `1`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000]`
- reps: `20`
- t0: `0.6`
- xi0: `[0.8]`
- x_grid_per_dim: `200`
- truth_grid_2d: `121`
- h0: `1.2`
- q: `0.7071067811865476`
- min_h_factor: `20.0`
- kappa: `2.0`
- selector_metric: `ise`
- trim_frac: `0.05`

## Theoretical benchmarks

- finite-range secant slope: `-0.347228`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `0.545303`
- oracle: slope = `-0.230709`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `4.079661e-01`
- lepski at M=2000: mean sup-grid error = `7.225983e-01`
- lepski at M=4000: mean sup-grid error = `8.688191e-01`
- oracle at M=1000: mean sup-grid error = `2.019592e-01`
- oracle at M=2000: mean sup-grid error = `2.278462e-01`
- oracle at M=4000: mean sup-grid error = `1.466774e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `4.808977e-01`
- lepski at M=2000: mean ISE = `8.779320e-01`
- lepski at M=4000: mean ISE = `1.048033e+00`
- oracle at M=1000: mean ISE = `2.187311e-01`
- oracle at M=2000: mean ISE = `2.446734e-01`
- oracle at M=4000: mean ISE = `1.706489e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `2.069613`
- M=2000: mean gap = `3.625427`
- M=4000: mean gap = `7.139106`