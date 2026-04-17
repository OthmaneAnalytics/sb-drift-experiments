# CLT experiment summary: GG 1D

- model_id: `gg_1d`
- dim: `1`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `300`
- t0: `0.6`
- x0: `[0.2]`
- xi0: `[0.0]`
- alpha: `0.22`
- c: `1.0`
- out_tag: `clt_screen_20260416_234045_0p22`
- qq_ms: `[4000, 8000]`

## Diagnostics

- M=1000: mean Z = `0.048896`, var Z = `0.982646`, coverage = `94.00`%, Shapiro p = `0.006900`, Anderson stat = `1.080506`, 5% crit = `0.750000`, reject@5% = `1`
- M=2000: mean Z = `-0.010723`, var Z = `1.105202`, coverage = `94.33`%, Shapiro p = `0.784646`, Anderson stat = `0.279905`, 5% crit = `0.750000`, reject@5% = `0`
- M=4000: mean Z = `0.024292`, var Z = `1.141197`, coverage = `93.33`%, Shapiro p = `0.645006`, Anderson stat = `0.385282`, 5% crit = `0.750000`, reject@5% = `0`
- M=8000: mean Z = `-0.049192`, var Z = `0.966526`, coverage = `96.00`%, Shapiro p = `0.080571`, Anderson stat = `1.001824`, 5% crit = `0.750000`, reject@5% = `1`