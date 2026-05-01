# Rate experiment summary: MM 1D

- model_id: `mm_1d`
- run_name: `mm_1d_adapt_dense_baseline_20260417_164626`
- dim: `1`
- beta_effective: `2.0`
- sample_sizes: `[1000, 1500, 2000, 3000, 4000, 6000, 8000]`
- reps: `50`
- t0: `0.6`
- xi0: `[0.8]`
- x_grid_per_dim: `200`
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

- finite-range secant slope: `-0.349379`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `-0.184445`
- oracle: slope = `-0.360552`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `2.652680e+00`
- lepski at M=1500: mean sup-grid error = `3.292377e+00`
- lepski at M=2000: mean sup-grid error = `2.544307e+00`
- lepski at M=3000: mean sup-grid error = `2.169311e+00`
- lepski at M=4000: mean sup-grid error = `2.561961e+00`
- lepski at M=6000: mean sup-grid error = `2.384035e+00`
- lepski at M=8000: mean sup-grid error = `1.790280e+00`
- oracle at M=1000: mean sup-grid error = `2.351797e-01`
- oracle at M=1500: mean sup-grid error = `2.058071e-01`
- oracle at M=2000: mean sup-grid error = `2.096382e-01`
- oracle at M=3000: mean sup-grid error = `1.445384e-01`
- oracle at M=4000: mean sup-grid error = `1.464118e-01`
- oracle at M=6000: mean sup-grid error = `1.288016e-01`
- oracle at M=8000: mean sup-grid error = `1.134872e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `3.379935e+00`
- lepski at M=1500: mean ISE = `4.153210e+00`
- lepski at M=2000: mean ISE = `3.194982e+00`
- lepski at M=3000: mean ISE = `2.742993e+00`
- lepski at M=4000: mean ISE = `3.205988e+00`
- lepski at M=6000: mean ISE = `3.024335e+00`
- lepski at M=8000: mean ISE = `2.227928e+00`
- oracle at M=1000: mean ISE = `2.627439e-01`
- oracle at M=1500: mean ISE = `2.313618e-01`
- oracle at M=2000: mean ISE = `2.364807e-01`
- oracle at M=3000: mean ISE = `1.671622e-01`
- oracle at M=4000: mean ISE = `1.715720e-01`
- oracle at M=6000: mean ISE = `1.432337e-01`
- oracle at M=8000: mean ISE = `1.285752e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `12.256032`
- M=1500: mean gap = `20.138569`
- M=2000: mean gap = `13.796766`
- M=3000: mean gap = `16.894370`
- M=4000: mean gap = `20.291089`
- M=6000: mean gap = `21.029999`
- M=8000: mean gap = `21.881373`