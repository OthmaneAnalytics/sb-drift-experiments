# CLT strengthening final report

- screen_tag: `clt_fast_screen_20260417_075529`
- final_tag: `clt_final_20260417_101207`

## Chosen alphas from the fast screen

```csv
model_id,alpha,tag,score,mean_abs_mean_Z,mean_abs_var_minus_1,mean_abs_cov_minus_95
gg_1d,0.22,clt_fast_screen_20260417_075529_0p22,0.13321183664557412,0.04765840843619535,0.05222009487604451,0.3333333333333428
mm_1d,0.28,clt_fast_screen_20260417_075529_0p28,0.17789845452273556,0.0583064452791128,0.052925342576955636,0.6666666666666714
```

## Final GG1 summary

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

## Final MM1 summary

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

## Priority 1 completion checklist

- [x] Higher-repetition final CLT runs completed (`R=300`).
- [x] QQ plots saved for `M=4000,8000`.
- [x] Formal normality diagnostics saved (Shapiro, Anderson).
- [x] Raw replicate-level CLT data saved (`pointwise_clt.csv`).
- [x] Additional normality table saved (`normality_table.csv`, `normality_table.md`).

## Bottom line

This CLT pass is substantially stronger than the earlier version. At the largest sample size (`M=8000`), both GG1 and MM1 show near-nominal coverage and non-rejection on the built-in Shapiro and Anderson diagnostics. The QQ plots and raw replicate data have been saved for later appendix/paper integration.