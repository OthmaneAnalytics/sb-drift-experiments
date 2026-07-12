# Numerical validation

Overall result: **PASS**

All calculations use fresh simulated data. Archived experiment tables are not validation inputs.

| Check | Result | Maximum error | Tolerance |
|---|---:|---:|---:|
| gg_1d.yaml: truth vs independent quadrature | PASS | 2.220446e-16 | 2.000000e-09 |
| mm_1d.yaml: truth vs independent quadrature | PASS | 5.551115e-17 | 2.000000e-09 |
| gg_2d.yaml: 2D truth grid refinement | PASS | 3.330669e-16 | 7.731798e-03 |
| mm_2d.yaml: 2D truth grid refinement | PASS | 2.220446e-16 | 1.139108e-02 |
| gg_1d.yaml: estimator direct reconstruction | PASS | 0.000000e+00 | 2.000000e-12 |
| gg_2d.yaml: estimator direct reconstruction | PASS | 0.000000e+00 | 2.000000e-12 |
| mm_1d.yaml: estimator direct reconstruction | PASS | 0.000000e+00 | 2.000000e-12 |
| mm_2d.yaml: estimator direct reconstruction | PASS | 0.000000e+00 | 2.000000e-12 |
| gg_1d.yaml: Lepski score reconstruction | PASS | 0.000000e+00 | 1.000000e-14 |
| mm_1d.yaml: Lepski score reconstruction | PASS | 0.000000e+00 | 1.000000e-14 |
| gg_1d.yaml: complete CLT row reconstruction | PASS | 5.551115e-17 | 2.000000e-11 |
| mm_1d.yaml: complete CLT row reconstruction | PASS | 8.326673e-17 | 2.000000e-11 |
| gg_1d.yaml: complete edge reconstruction | PASS | 0.000000e+00 | 2.000000e-11 |
| mm_1d.yaml: complete edge reconstruction | PASS | 0.000000e+00 | 2.000000e-11 |

## Interpretation

- Truth values were checked independently in 1D and by grid refinement in 2D.
- Drift estimates were reconstructed directly from kernel weights, numerator, and denominator.
- Adaptive bandwidth scores were independently reconstructed.
- CLT output rows were regenerated from their recorded sample seeds and compared quantity by quantity.
- Terminal-edge bandwidths and errors were independently reconstructed.
