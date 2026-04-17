# CLT experiment summary: MM 1D

- model_id: `mm_1d`
- dim: `1`
- sample_sizes: `[1000, 2000, 4000, 8000]`
- reps: `300`
- t0: `0.6`
- x0: `[0.3]`
- xi0: `[0.8]`
- alpha: `0.28`
- c: `1.0`
- out_tag: `clt_final_20260417_101207`
- qq_ms: `[4000, 8000]`

## Diagnostics

- M=1000: mean Z = `0.011673`, var Z = `1.180449`, coverage = `92.33`%, Shapiro p = `0.094434`, Anderson stat = `0.487901`, 5% crit = `0.750000`, reject@5% = `0`
- M=2000: mean Z = `-0.054668`, var Z = `0.964488`, coverage = `96.67`%, Shapiro p = `0.139434`, Anderson stat = `0.500433`, 5% crit = `0.750000`, reject@5% = `0`
- M=4000: mean Z = `0.133135`, var Z = `0.879781`, coverage = `96.33`%, Shapiro p = `0.905583`, Anderson stat = `0.173670`, 5% crit = `0.750000`, reject@5% = `0`
- M=8000: mean Z = `0.041255`, var Z = `0.894031`, coverage = `96.00`%, Shapiro p = `0.466130`, Anderson stat = `0.414661`, 5% crit = `0.750000`, reject@5% = `0`