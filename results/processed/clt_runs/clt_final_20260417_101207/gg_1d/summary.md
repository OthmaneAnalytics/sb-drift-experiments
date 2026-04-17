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
- out_tag: `clt_final_20260417_101207`
- qq_ms: `[4000, 8000]`

## Diagnostics

- M=1000: mean Z = `-0.029120`, var Z = `1.160508`, coverage = `92.67`%, Shapiro p = `0.139036`, Anderson stat = `0.689650`, 5% crit = `0.750000`, reject@5% = `0`
- M=2000: mean Z = `-0.083007`, var Z = `0.964855`, coverage = `96.33`%, Shapiro p = `0.929805`, Anderson stat = `0.253814`, 5% crit = `0.750000`, reject@5% = `0`
- M=4000: mean Z = `0.056018`, var Z = `0.828702`, coverage = `96.67`%, Shapiro p = `0.180917`, Anderson stat = `0.536864`, 5% crit = `0.750000`, reject@5% = `0`
- M=8000: mean Z = `-0.015311`, var Z = `0.974488`, coverage = `96.67`%, Shapiro p = `0.584810`, Anderson stat = `0.483753`, 5% crit = `0.750000`, reject@5% = `0`