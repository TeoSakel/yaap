# yaap

yaap fits archetypal analysis models for numeric matrix data. It includes Euclidean, probabilistic, kernel, and directional variants, plus initialization helpers, model selection paths, plotting methods, broom methods, and a tidymodels recipe step.

## Installation

Install the released version from CRAN:

```r
install.packages("yaap")
```

You can also install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("teosakel/yaap")
```

## Features

- [Principal archetypal analysis (PCHA)](https://doi.org/10.1016/j.neucom.2011.06.033),
  with gradient descent and NNLS-based fitting methods
- [Probabilistic archetypes](https://arxiv.org/abs/1312.7604)
  (Bernoulli, Poisson, Multinomial)
- [Kernel archetypes](https://doi.org/10.1016/j.neucom.2011.06.033)
- [Directional archetypes](https://doi.org/10.3389/fnins.2022.911034)
- Multiple initialization helpers, including `FurthestSum` and
  [`AA++`](https://arxiv.org/abs/2301.13748)
- Plotting methods, broom methods, and a tidymodels recipe step
