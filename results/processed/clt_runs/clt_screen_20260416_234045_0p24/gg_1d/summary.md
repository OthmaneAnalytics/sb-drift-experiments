# CLT experiment summary: GG 1D

- model_id: `gg_1d`
- dim: `1`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `300`
- t0: `0.6`
- x0: `[0.2]`
- xi0: `[0.0]`
- alpha: `0.24`
- c: `1.0`
- out_tag: `clt_screen_20260416_234045_0p24`
- qq_ms: `[4000, 8000]`

## Diagnostics

- M=1000: mean Z = `0.049614`, var Z = `0.985824`, coverage = `94.33`%, Shapiro p = `0.003042`, Anderson stat = `1.316861`, 5% crit = `0.750000`, reject@5% = `1`
- M=2000: mean Z = `-0.010579`, var Z = `1.103481`, coverage = `94.33`%, Shapiro p = `0.978996`, Anderson stat = `0.167570`, 5% crit = `0.750000`, reject@5% = `0`
- M=4000: mean Z = `0.012835`, var Z = `1.136358`, coverage = `93.33`%, Shapiro p = `0.375248`, Anderson stat = `0.457155`, 5% crit = `0.750000`, reject@5% = `0`
- M=8000: mean Z = `-0.046546`, var Z = `1.003004`, coverage = `95.67`%, Shapiro p = `0.138975`, Anderson stat = `0.837183`, 5% crit = `0.750000`, reject@5% = `1`