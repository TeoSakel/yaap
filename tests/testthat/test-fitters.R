test_that("archetypes_pgd fits toy data with expected invariants", {
    set.seed(1)
    X <- toy_matrix()
    fit <- archetypes_pgd(X, K = 3L, max_iter = 20L, tol_r2 = 0.95)

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["coordinates"]], 3L, 2L)
    expect_matrix_dim(fit[["coefficients"]], 3L, 250L)
    expect_matrix_dim(fit[["compositions"]], 250L, 3L)
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_equal(dim(fitted(fit)), dim(X))
    expect_equal(dim(residuals(fit)), dim(X))
})

test_that("archetypes_nnls fits toy data with expected invariants", {
    set.seed(1)
    X <- toy_matrix()
    fit <- archetypes_nnls(X, K = 3L, max_iter = 5L, tol_r2 = 0.95)

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["coordinates"]], 3L, 2L)
    expect_matrix_dim(fit[["coefficients"]], 3L, 250L)
    expect_matrix_dim(fit[["compositions"]], 250L, 3L)
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_equal(dim(fitted(fit)), dim(X))
    expect_equal(dim(residuals(fit)), dim(X))
})
