# AGENTS.md

Guidance for Codex and other coding agents working on `yaap`.

`yaap` is an R package for archetypal analysis. Keep changes package-local,
matrix-centric, and compatible with ordinary R package development workflows.
Do not assume tidyverse-first data frame semantics; the core API is built
around numeric matrices, `Matrix` sparse matrices, S3 classes, and optional
formula or specialized object front ends.

## Project Shape

- Core fitting entry point: `R/run_aa.R`.
- Direct solver wrappers: `R/archetypes_*.R`.
- Base S3 object and methods: `R/archetypes_class.R`.
- Initialization helpers: `R/init.R`.
- Matrix/simplex utilities: `R/projections.R`, `R/utils.R`.
- Plotting, broom, and recipe interfaces: `R/plotting.R`, `R/broom.R`,
  `R/step_archetypes.R`.
- Tests: `tests/testthat/`, with shared helpers in
  `tests/testthat/helper-toy-data.R`.
- Vignettes: `vignettes/`.

Use existing package patterns before introducing new abstractions. Avoid broad
refactors unless they are required for the current change.

## R Package Conventions

- Write base-R-compatible package code. Use imported packages only when they
  are already listed in `DESCRIPTION` or when adding the dependency is clearly
  justified.
- Keep core computations matrix-based. Accept data frames only at API
  boundaries where the package already does so, then convert explicitly.
- Do not make tidyverse compatibility a design goal. `tibble` output is
  appropriate for existing broom-style methods, but not for internal solver
  state or core model objects.
- Prefer S3 methods and existing generics over ad hoc dispatch.
- Use `stop(..., call. = FALSE)` or existing validation helpers for user-facing
  errors. Error messages should name the bad argument and the expected shape or
  value.
- Preserve row names, column names, and archetype names where the surrounding
  code already does so.
- Use BLAS/LAPACK-backed matrix operations when practical: `crossprod()`,
  `tcrossprod()`, `chol()`, `qr()`, `backsolve()`, `forwardsolve()`, and matrix
  products are preferred over hand-written numeric loops when they are clear
  and stable.
- Be careful with dimensions. The package convention is rows = samples and
  columns = features.

## Archetype Object Invariants

Objects returned by fitters should be built through the existing constructors
and should satisfy these invariants:

- `coordinates`: `K x M`, unless a specialized class intentionally cannot
  store input-space coordinates.
- `coefficients`: `K x N`; rows are stochastic unless `slack > 0` permits a
  documented relaxation.
- `compositions`: `N x K`; rows are stochastic.
- `loss`: a data frame with per-iteration metrics. `.aa_final_loss()` indexes
  with `fit[["i"]] + 1L`.
- `call`: preserves the public entry point. Direct wrappers keep
  `archetypes_<method>`; `run_aa()` calls keep `run_aa`.
- `family`: a single non-empty string, usually `"gaussian"` unless the method
  uses a different observation family.
- `fit_info`: a list containing at least `method = "<method>"` plus meaningful
  metadata for kernel, directional, robust, missing-data, scaling, weighting,
  or family-specific behavior.

When adding a specialized class, add or update methods for the behaviors that
differ from the base `archetypes` class, such as `predict()`, `fitted()`,
`residuals()`, `AIC()`, `tidy()`, `glance()`, or plotting.

## Adding Or Changing Fitters

Before changing fitter behavior, inspect the closest existing implementation:

- Euclidean squared-error solvers: `R/archetypes_pgd.R` and
  `R/archetypes_nnls.R`.
- Alternative observation families: `R/archetypes_paa.R`.
- Directional geometry: `R/archetypes_directional.R`.
- Kernel-space methods: `R/archetypes_kernel_pgd.R`.
- Shared dispatch and preprocessing: `R/run_aa.R`.

Follow the engine-block pattern. A fitter block should return named functions:
`check`, `preprocess`, `edge_case`, `init`, `fit`, `final_loss`, and
`prepare_output`.

For a new solver:

1. Add or extend tests first.
2. Add the exported `archetypes_<method>()` wrapper with roxygen docs.
3. Implement `.aa_<method>_block(ctx, ...)`.
4. Validate unsupported inherited options in `check`.
5. Reuse shared Euclidean helpers when applicable:
   `.aa_euclidean_check()`, `.aa_euclidean_preprocess()`,
   `.aa_euclidean_edge_case()`, `.aa_euclidean_init()`, and
   `.aa_euclidean_output()`.
6. Wire `run_aa.default()` and `.aa_fit_engine()` method choices and switch dispatch.
7. Decide whether `tune_archetypes()`, `step_archetypes()`, broom methods,
   plotting, or vignettes need updates.
8. Regenerate documentation and `NAMESPACE` using `devtools::document()` if roxygen changed.

If translating an algorithm from a paper, add concise comments near the solver
that map code variables to paper notation, especially when matrices are
transposed relative to the paper.

## Tests

Use `testthat` edition 3 conventions already present in the package.

For fitter changes, cover:

- Direct wrapper returns the expected class and dimensions.
- `run_aa(..., method = "<method>")` dispatches correctly.
- User-facing call preservation.
- `coefficients` and `compositions` stochasticity.
- Finite loss or a documented exception.
- Invalid inputs and unsupported options produce clear errors.
- Edge cases such as `K = 1`, `K = nrow(x)`, sparse input, missing data, or
  method-specific incompatibilities when relevant.
- Specialized S3 methods when class behavior differs from base `archetypes`.

Use existing helpers such as `toy_matrix()`, `expect_archetypes_fit()`,
`expect_matrix_dim()`, and `expect_row_stochastic()`.

For optional dependencies in `Suggests`, use `testthat::skip_if_not_installed()`
or `requireNamespace(..., quietly = TRUE)` as appropriate. Do not make optional
packages hard runtime requirements unless `DESCRIPTION` is updated deliberately.

## Documentation

- Use roxygen2 comments for exported functions and S3 methods.
- Keep examples CRAN-safe: small, fast, deterministic, and without network access.
- Put slow examples, optional package examples, or illustrative workflows in
  `\dontrun{}` or vignettes as appropriate.
- After changing roxygen docs or exports, run:

```r
devtools::document()
```

Then review generated changes in `NAMESPACE` and `man/*.Rd`.

## CRAN Compatibility

Keep the package suitable for eventual CRAN submission:

- No network access in tests, examples, or vignettes.
- No writing outside temp directories during tests.
- Keep examples and tests reasonably fast.
- Guard optional dependencies from `Suggests`.
- Avoid relying on local files except through `system.file()` or checked-in
  test fixtures.
- Avoid changing user options, working directories, RNG state, or global state
  without restoring them. Prefer `withr` where helpful.
- Use deterministic seeds in stochastic tests.
- Keep package startup quiet.
- Do not use non-exported functions from other packages.
- Do not add compiled code, system dependencies, or external tools without
  documenting and checking the CRAN impact.

Before considering a change complete, run a verification level appropriate to
the scope:

```r
devtools::test()
```

For documentation/export changes:

```r
devtools::document()
devtools::test()
```

For release-like or CRAN-sensitive changes:

```r
devtools::check(args = "--as-cran")
```

The command-line equivalents are also valid:

```sh
R CMD build .
R CMD check --as-cran yaap_*.tar.gz
```

## Development Workflow

- Prefer focused edits and focused tests first, then broader checks.
- Keep generated documentation changes bundled with the code change that caused them.
- Do not manually edit generated `man/*.Rd` files when roxygen source should be changed instead.
- Do not reorder `NAMESPACE` or churn generated files without a roxygen reason.
- Check `git diff` before finishing and make sure unrelated user changes are left intact.
- If a command fails, inspect the failure before changing code. Fix the cause, not the symptom.

## Style Notes

- Follow the existing four-space indentation style in R files.
- Run `styler` with `strict=FALSE` to reformat the code after you write it.
- Use explicit `return()` only when exiting early; match the surrounding style.
- Prefer `[[` for list component access in internal code when consistency or
  partial-match avoidance matters.
- Keep comments factual and useful. Avoid comments that restate obvious code.
- Keep public argument names stable unless the user explicitly asks for an API break.
- Avoid writing tiny functions that are only used once. Prefer short comments
  headers for section to explain what they do or use local/anonymous functions when
  a function is required as an argument
