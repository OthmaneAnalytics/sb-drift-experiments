# Adaptivity decision note

The canonical baseline selector remains:

- selector_metric = raw_max
- penalty_form = one_sided
- trim_frac = 0.05
- kappa_pair = 2.0
- kappa_final = 2.0

Reason:
- the saved finalized runs (`*_final_rawmax_k2`) reproduce the paper numbers;
- the later `raw1s_k*` runs are not protocol-compatible with that baseline;
- among the saved nearby 1D variants, no alternative improves the hard MM1 case relative to the canonical baseline.

Practical conclusion:
- keep the canonical baseline selector for the paper;
- improve the presentation by adding selector-comparison evidence and uncertainty ribbons;
- explicitly acknowledge MM1 as the hard case.
