# Rate experiment summary: MM 1D

- model_id: `mm_1d`
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
- kappa: `1.75`

## Theoretical benchmarks

- finite-range secant slope: `-0.347228`
- asymptotic slope: `-0.400000`

## Fitted slopes: sup-grid error

- lepski: slope = `0.415828`
- oracle: slope = `-0.230709`

## Mean sup-grid error by sample size

- lepski at M=1000: mean sup-grid error = `4.375159e-01`
- lepski at M=2000: mean sup-grid error = `7.306031e-01`
- lepski at M=4000: mean sup-grid error = `7.786588e-01`
- oracle at M=1000: mean sup-grid error = `2.019592e-01`
- oracle at M=2000: mean sup-grid error = `2.278462e-01`
- oracle at M=4000: mean sup-grid error = `1.466774e-01`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `4.892749e-01`
- lepski at M=2000: mean ISE = `8.793129e-01`
- lepski at M=4000: mean ISE = `9.353825e-01`
- oracle at M=1000: mean ISE = `2.187311e-01`
- oracle at M=2000: mean ISE = `2.446734e-01`
- oracle at M=4000: mean ISE = `1.706489e-01`

## Adaptive/oracle gap (sup-grid error)

- M=1000: mean gap = `2.165259`
- M=2000: mean gap = `3.706144`
- M=4000: mean gap = `6.447463`