# Rate experiment summary: MM 1D

- model_id: `mm_1d`
- run_name: `mm1_final_rawmax_k2`
- dim: `1`
- beta_effective: `2.0`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `50`
- t0: `0.6`
- xi0: `[0.8]`
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

- lepski: slope = `-0.110769`
- oracle: slope = `-0.390229`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `5.236430e-01`
- lepski at M=2000: mean sup-grid error = `4.757768e-01`
- lepski at M=4000: mean sup-grid error = `6.406517e-01`
- lepski at M=8000: mean sup-grid error = `3.671245e-01`
- oracle at M=1000: mean sup-grid error = `2.343623e-01`
- oracle at M=2000: mean sup-grid error = `2.027866e-01`
- oracle at M=4000: mean sup-grid error = `1.363926e-01`
- oracle at M=8000: mean sup-grid error = `1.085763e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `6.137988e-01`
- lepski at M=2000: mean ISE = `5.272210e-01`
- lepski at M=4000: mean ISE = `7.566674e-01`
- lepski at M=8000: mean ISE = `4.236832e-01`
- oracle at M=1000: mean ISE = `2.623407e-01`
- oracle at M=2000: mean ISE = `2.245352e-01`
- oracle at M=4000: mean ISE = `1.521577e-01`
- oracle at M=8000: mean ISE = `1.232540e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `2.380483`
- M=2000: mean gap = `2.403346`
- M=4000: mean gap = `5.628425`
- M=8000: mean gap = `3.215650`