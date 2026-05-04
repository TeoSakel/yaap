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

test_that("Tukey row weights downweight large row residuals", {
    weights <- .aa_bisquare_weights(c(0, 1, 4, 1e6))

    expect_equal(weights[1], 1)
    expect_true(weights[4] < weights[2])
    expect_equal(.aa_bisquare_weights(rep(0, 4)), rep(1, 4))
})

test_that("robust archetypes fitters keep expected invariants", {
    set.seed(1)
    X <- toy_matrix()
    X[1, ] <- X[1, ] + 25

    pgd <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        robust = TRUE,
        max_iter = 5L,
        tol_r2 = 0.95
    ))
    nnls <- suppressWarnings(archetypes_nnls(
        X,
        K = 3L,
        robust = TRUE,
        max_iter = 3L,
        tol_r2 = 0.95
    ))

    for (fit in list(pgd, nnls)) {
        expect_s3_class(fit, "archetypes")
        expect_matrix_dim(fit[["coordinates"]], 3L, 2L)
        expect_matrix_dim(fit[["coefficients"]], 3L, 250L)
        expect_matrix_dim(fit[["compositions"]], 250L, 3L)
        expect_row_stochastic(fit[["coefficients"]])
        expect_row_stochastic(fit[["compositions"]])
        expect_true(all(is.finite(fit[["loss"]][["rss"]])))
    }
})

test_that("archetypes fitters accept named coordinate matrix initialization", {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1,
            0.25, 0.25
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(paste0("x", 1:4), c("u", "v"))
    )
    init <- X[1:3, , drop = FALSE]
    rownames(init) <- c("left", "right", "top")

    pgd <- suppressWarnings(archetypes_pgd(X, K = 3L, init = init, max_iter = 1L))
    nnls <- suppressWarnings(archetypes_nnls(X, K = 3L, init = init, max_iter = 1L))

    for (fit in list(pgd, nnls)) {
        expect_named(as.data.frame(fit[["coordinates"]]), colnames(X))
        expect_equal(rownames(fit[["init"]]), rownames(init))
        expect_equal(rownames(fit[["coordinates"]]), rownames(init))
        expect_equal(rownames(fit[["coefficients"]]), rownames(init))
        expect_equal(colnames(fit[["compositions"]]), rownames(init))
        expect_equal(colnames(fit[["coefficients"]]), rownames(X))
        expect_equal(rownames(fit[["compositions"]]), rownames(X))
    }
})

test_that("coordinate matrix initialization validates dimensions and projects to convex hull", {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1,
            0.25, 0.25
        ),
        ncol = 2,
        byrow = TRUE
    )

    bad_nrow <- X[1:2, , drop = FALSE]
    bad_ncol <- cbind(X[1:3, , drop = FALSE], z = 0)
    outside <- X[1:3, , drop = FALSE]
    outside[1, ] <- c(2, 2)

    expect_error(archetypes_pgd(X, K = 3L, init = bad_nrow), "nrow\\(init\\)")
    expect_error(archetypes_pgd(X, K = 3L, init = bad_ncol), "ncol\\(init\\)")
    expect_warning(
        expect_warning(
            fit <- archetypes_pgd(X, K = 3L, init = outside, max_iter = 1L),
            "projected"
        ),
        "did not converge"
    )
    expect_equal(dim(fit[["init"]]), c(3L, 2L))
})

test_that("PGD coordinate matrix initialization honors delta-relaxed hull", {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1
        ),
        ncol = 2,
        byrow = TRUE
    )
    init <- X
    init[1, ] <- c(1.2, 0)

    expect_warning(
        strict <- .aa_matrix_init(X, K = 3L, init = init, eps = 0, delta = 0),
        "projected"
    )
    relaxed <- .aa_matrix_init(X, K = 3L, init = init, eps = 0, delta = 0.25)

    expect_equal(unname(strict[["A"]][1, ]), c(1, 0), tolerance = 1e-8)
    expect_equal(unname(relaxed[["A"]][1, ]), c(1.2, 0), tolerance = 1e-8)
    expect_true(rowSums(relaxed[["B"]])[1] <= 1.25 + 1e-8)
})
