test_that("archetypes_directional fits spherical data with expected invariants", {
    set.seed(1)
    X <- directional_matrix()
    fit <- archetypes_directional(X, K = 3L, max_iter = 10L, tol_r2 = 0.9)

    expect_s3_class(fit, "directional_archetypes")
    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["coordinates"]], 3L, 3L)
    expect_matrix_dim(fit[["directions"]], 3L, 3L)
    expect_matrix_dim(fit[["coefficients"]], 3L, nrow(X))
    expect_matrix_dim(fit[["compositions"]], nrow(X), 3L)
    expect_matrix_dim(fit[["generator_data"]], nrow(X), ncol(X))
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_true(is_all_finite(fit[["loss"]][["loss"]]))
    expect_named(fit[["loss"]], c("loss", "r2"))
    expect_true(is_all_finite(fitted(fit)))
    expect_true(is_all_finite(residuals(fit)))
    S_pred <- predict(fit, X[1:5, ], max_iter = 3L)
    Y_pred <- predict(fit, X[1:5, ], type = "reconstruction", max_iter = 3L)
    expect_matrix_dim(S_pred, 5L, 3L)
    expect_matrix_dim(Y_pred, 5L, ncol(X))
    expect_equal(as.vector(rowSums(Y_pred^2)), rep(1, nrow(Y_pred)), tolerance = 1e-6)
    expect_error(AIC(fit), "not defined")
})

test_that("directional loss is monotone non-increasing", {
    set.seed(2)
    fit <- archetypes_directional(
        directional_matrix(),
        K = 3L,
        max_iter = 12L,
        tol_r2 = 1
    )
    expect_true(all(diff(fit[["loss"]][["loss"]]) <= 1e-8))
})

test_that("directional AA is stable to polarity flips", {
    X <- directional_matrix()
    signs <- rep(c(-1, 1), length.out = nrow(X))

    set.seed(3)
    fit <- archetypes_directional(X, K = 3L, max_iter = 12L, tol_r2 = 1)
    set.seed(3)
    fit_flip <- archetypes_directional(
        X * signs,
        K = 3L,
        max_iter = 12L,
        tol_r2 = 1
    )

    final_r2 <- tail(fit[["loss"]][["r2"]], 1L)
    final_r2_flip <- tail(fit_flip[["loss"]][["r2"]], 1L)
    expect_equal(final_r2, final_r2_flip, tolerance = 0.15)
})

test_that("run_aa dispatches to directional fitter", {
    set.seed(4)
    X <- directional_matrix()
    fit <- run_aa(
        X,
        K = 3L,
        method = "directional",
        max_iter = 5L,
        tol_r2 = 0.9
    )

    expect_s3_class(fit, "directional_archetypes")
    expect_identical(as.character(fit[["call"]][[1L]]), "run_aa")
    expect_matrix_dim(fit[["coordinates"]], 3L, 3L)
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
})

test_that("directional AA validates unsupported and invalid inputs", {
    X <- directional_matrix()
    X_zero <- X
    X_zero[1L, ] <- 0
    X_na <- X
    X_na[1L, 1L] <- NA_real_

    expect_error(archetypes_directional(X_zero, K = 3L), "zero-norm")
    expect_error(archetypes_directional(X_na, K = 3L), "missing")
    expect_error(
        archetypes_directional(Matrix::Matrix(X, sparse = TRUE), K = 3L),
        "dense numeric matrix"
    )
    expect_error(
        run_aa(X, K = 3L, method = "directional", missing = TRUE),
        "missing = TRUE"
    )
    expect_error(
        run_aa(X, K = 3L, method = "directional", robust = TRUE),
        "robust"
    )
    expect_error(
        run_aa(X, K = 3L, method = "directional", robust = "psi.huber"),
        "robust"
    )
    expect_warning(
        fit <- run_aa(X, K = 3L, method = "directional", scale = TRUE),
        "ignored"
    )
    expect_s3_class(fit, "directional_archetypes")
})
