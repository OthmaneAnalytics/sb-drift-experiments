# Dense baseline adaptivity rerun

This rerun keeps the canonical baseline selector fixed:

- selector_metric = raw_max
- penalty_form = one_sided
- kappa_pair = 2.0
- kappa_final = 2.0

and evaluates it on a denser geometric sample-size grid

- M in {1000, 1500, 2000, 3000, 4000, 6000, 8000}

for the purpose of producing a smoother and more informative adaptivity-gap figure.

The visualization uses monotone interpolation (PCHIP) on the log-M scale together with interquartile ribbons.
