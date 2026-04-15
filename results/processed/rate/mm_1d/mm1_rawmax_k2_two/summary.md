# Rate experiment summary: MM 1D

- model_id: `mm_1d`
- run_name: `mm1_rawmax_k2_two`
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
- selector_metric: `raw_max`
- trim_frac: `0.05`
- penalty_form: `two_sided`

## Theoretical benchmarks

- finite-range secant slope: `-0.347228`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `0.397489`
- oracle: slope = `-0.230709`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `4.465361e-01`
- lepski at M=2000: mean sup-grid error = `6.381224e-01`
- lepski at M=4000: mean sup-grid error = `7.747626e-01`
- oracle at M=1000: mean sup-grid error = `2.019592e-01`
- oracle at M=2000: mean sup-grid error = `2.278462e-01`
- oracle at M=4000: mean sup-grid error = `1.466774e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `4.669061e-01`
- lepski at M=2000: mean ISE = `7.183495e-01`
- lepski at M=4000: mean ISE = `9.404964e-01`
- oracle at M=1000: mean ISE = `2.187311e-01`
- oracle at M=2000: mean ISE = `2.446734e-01`
- oracle at M=4000: mean ISE = `1.706489e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `2.262488`
- M=2000: mean gap = `3.672942`
- M=4000: mean gap = `6.270709`