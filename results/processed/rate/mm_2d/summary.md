# Rate experiment summary: MM 2D

- model_id: `mm_2d`
- dim: `2`
- sample_sizes: `[1000, 2000, 4000]`
- reps: `10`
- t0: `0.6`
- xi0: `[0.8, -0.8]`
- x_grid_per_dim: `21`
- truth_grid_2d: `81`
- h0: `1.2`
- q: `0.7071067811865476`
- min_h_factor: `5.0`
- kappa: `2.5`

## Fitted slopes

- lepski: slope = `-0.212270`
- oracle: slope = `-0.465828`

## Mean ISE by sample size

- lepski at M=1000: mean ISE = `9.049201e-01`
- lepski at M=2000: mean ISE = `7.258615e-01`
- lepski at M=4000: mean ISE = `6.742344e-01`
- oracle at M=1000: mean ISE = `7.265266e-01`
- oracle at M=2000: mean ISE = `4.976644e-01`
- oracle at M=4000: mean ISE = `3.808862e-01`