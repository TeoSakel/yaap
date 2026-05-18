# yaap

Yet Another Archetypes Package.

This R package implements different flavors of [Archetypal Analysis](https://arxiv.org/abs/2504.12392).

## Currently implemented

- [Principal archetypal analysis (PCHA)](https://doi.org/10.1016/j.neucom.2011.06.033),
  with gradient descent and NNLS-based fitting methods
- [Probabilistic archetypes](https://arxiv.org/abs/1312.7604)
  (Bernoulli, Poisson, Multinomial)
- [Kernel archetypes](https://doi.org/10.1016/j.neucom.2011.06.033)
- [Directional archetypes](https://doi.org/10.3389/fnins.2022.911034)
- Multiple initialization helpers, including `FurthestSum` and
  [`AA++`](https://arxiv.org/abs/2301.13748)
- Plotting methods, broom methods, and a tidymodels recipe step

## Work in progress

- Vignettes and examples
- Tools for tuning several archetype fits
- Tools for comparing and combining ensembles of archetype fits
- Publish on CRAN

The project is still a work in progress. Interfaces, documentation, and
examples may change as the package develops.
