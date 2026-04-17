# Rate experiment summary: GG 1D

- model_id: `gg_1d`
- run_name: `gg_1d_adapt_dense_baseline_20260417_164626`
- dim: `1`
- beta_effective: `2.0`
- sample_sizes: `[1000, 1500, 2000, 3000, 4000, 6000, 8000]`
- reps: `50`
- t0: `0.6`
- xi0: `[0.0]`
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

- lepski: slope = `-0.330285`
- oracle: slope = `-0.333658`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `2.442480e-01`
- lepski at M=1500: mean sup-grid error = `2.287218e-01`
- lepski at M=2000: mean sup-grid error = `2.115467e-01`
- lepski at M=3000: mean sup-grid error = `2.179639e-01`
- lepski at M=4000: mean sup-grid error = `1.751263e-01`
- lepski at M=6000: mean sup-grid error = `1.410353e-01`
- lepski at M=8000: mean sup-grid error = `1.225496e-01`
- oracle at M=1000: mean sup-grid error = `1.618823e-01`
- oracle at M=1500: mean sup-grid error = `1.361751e-01`
- oracle at M=2000: mean sup-grid error = `1.144949e-01`
- oracle at M=3000: mean sup-grid error = `1.275259e-01`
- oracle at M=4000: mean sup-grid error = `1.041024e-01`
- oracle at M=6000: mean sup-grid error = `9.376059e-02`
- oracle at M=8000: mean sup-grid error = `7.212805e-02`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `2.331669e-01`
- lepski at M=1500: mean ISE = `2.203946e-01`
- lepski at M=2000: mean ISE = `2.064334e-01`
- lepski at M=3000: mean ISE = `2.080068e-01`
- lepski at M=4000: mean ISE = `1.665750e-01`
- lepski at M=6000: mean ISE = `1.324814e-01`
- lepski at M=8000: mean ISE = `1.153186e-01`
- oracle at M=1000: mean ISE = `1.592781e-01`
- oracle at M=1500: mean ISE = `1.393137e-01`
- oracle at M=2000: mean ISE = `1.196620e-01`
- oracle at M=3000: mean ISE = `1.214047e-01`
- oracle at M=4000: mean ISE = `1.025378e-01`
- oracle at M=6000: mean ISE = `9.345876e-02`
- oracle at M=8000: mean ISE = `7.325276e-02`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `1.683205`
- M=1500: mean gap = `1.864515`
- M=2000: mean gap = `2.154159`
- M=3000: mean gap = `2.108757`
- M=4000: mean gap = `1.914822`
- M=6000: mean gap = `1.628153`
- M=8000: mean gap = `1.899475`