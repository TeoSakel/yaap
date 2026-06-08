# yaap

`yaap` is **Yet Another Archetypes Package**: a practical, matrix-first R
toolkit that brings together many variants and flavors of archetypal analysis
in one place. Some of these workflows are currently scattered across different
packages, while others have had little or no R implementation available.

The package centers on a shared `run_aa()` interface for fitting, inspecting,
comparing, plotting, and reusing archetypal analysis models. It supports
multivariate, functional, non-Gaussian, kernel, and directional AA, plus direct
solver wrappers when you want more control.

## Installation 🛠️

Install the released version from CRAN:

```r
install.packages("yaap")
```

You can install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("teosakel/yaap")
```

## Quick Start 🚀

The core workflow is matrix-based: rows are samples, columns are features, and
archetypes are fitted as extreme profiles in the convex hull of the data.

```r
library(yaap)

X <- as.matrix(iris[, 1:4])  # numeric columns

fit <- run_aa(X, K = 3, nrep = 5, scale = TRUE)

coordinates(fit)          # K x features: the archetype profiles
compositions(fit)[1:6, ]  # samples x K: each sample's archetype mixture

plot(fit, what = "profiles")
generics::glance(fit)
```

You can also fit by formula when that is more convenient:

```r
fit_iris <- run_aa(Species ~ ., data = iris, K = 3, scale = TRUE)
```

## What Can `yaap` Do? 🧰

- 📐 **Euclidean / Gaussian AA** with projected-gradient and NNLS solvers via
  `run_aa(..., method = "pgd")`, `run_aa(..., method = "nnls")`,
  `archetypes_pgd()`, and `archetypes_nnls()`.
- 🧮 **Probabilistic AA** for Gaussian, binomial, Poisson, and multinomial data
  via `run_aa(..., method = "paa", family = ...)`.
- 🌀 **Kernel AA** for nonlinear geometry via
  `run_aa(..., method = "kernel", kernel = ...)`, including precomputed kernels.
- 🧭 **Directional AA** for unit-vector, angular, or polarity-invariant data via
  `run_aa(..., method = "directional")`.
- 📈 **Metric and functional AA** using feature metrics through `scale = G` and
  direct support for `fda::fd` objects.
- 🧪 **Robust and missing-data workflows** for Gaussian AA, including automatic
  missing-data handling when `NA` values are present.
- 🧩 **Initialization helpers** through `aa_init()`, including random,
  Dirichlet, furthest-first, k-means++, FurthestSum, AA++, batched coreset-style
  initialization, and hull-outmost strategies.
- 📊 **Model comparison and diagnostics** with `archetypes_path()`,
  `screeplot()`, `AIC()`, loss tracking, plotting methods, and standard S3
  helpers such as `predict()`, `fitted()`, and `residuals()`.
- 🧹 **Workflow integration** through `generics::tidy()`, `generics::glance()`,
  `generics::augment()`, and the `recipes` step `step_archetypes()`.

## Vignette Tour 🗺️

The README is only the map. The vignettes are the snacks.

- [Introduction](vignettes/introduction.Rmd): the main `run_aa()` workflow,
  object structure, plotting, prediction, choosing `K`, robust fitting, missing
  data, and PGD vs NNLS.
- [Initialization](vignettes/initialization.Rmd): how `aa_init()` works, when
  initialization matters, and how FurthestSum, AA++, batched, and hull-based methods behave.
- [Non-Gaussian and Alternative Geometries](vignettes/non_gaussian_aa.Rmd):
  metric Gaussian AA, functional AA, kernel AA, probabilistic AA, and directional AA.
- [Tidymodels](vignettes/tidymodels.Rmd): using `step_archetypes()` inside
  `recipes` workflows and tuning archetype-related parameters.

## Methods and References 📚

`yaap` keeps the references close to the features they support:

- **Classical archetypal analysis**:
  [Cutler and Breiman (1994)](https://doi.org/10.1080/00401706.1994.10485840).
- **PCHA, FurthestSum, and kernel AA**:
  [Mørup and Hansen (2012)](https://doi.org/10.1016/j.neucom.2011.06.033).
- **Probabilistic archetypal analysis**:
  [Seth and Eugster (2016)](https://doi.org/10.1007/s10994-015-5498-8).
- **Directional archetypal analysis**:
  [Olsen et al. (2022)](https://doi.org/10.3389/fnins.2022.911034).
- **AA++ initialization**:
  [Mair and Sjölund (2023)](https://arxiv.org/abs/2301.13748).
- **Coreset-style initialization**:
  [Mair and Brefeld (2019)](https://proceedings.neurips.cc/paper_files/paper/2019/file/7f278ad602c7f47aa76d1bfc90f20263-Paper.pdf).
- **Validation and adapted AIC**:
  [Suleman (2017)](https://doi.org/10.1109/FUZZ-IEEE.2017.8015385).
- **A broader survey of AA methods**:
  [Alcacer et al. (2025)](https://arxiv.org/abs/2504.12392).
