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

test_that("run_aa dispatches to supported fitters", {
    set.seed(1)
    X <- toy_matrix()
    pgd <- run_aa(X, K = 3L, method = "pgd", max_iter = 20L, tol_r2 = 0.95)

    set.seed(1)
    nnls <- run_aa(X, K = 3L, method = "nnls", max_iter = 5L, tol_r2 = 0.95)

    for (fit in list(pgd, nnls)) {
        expect_s3_class(fit, "archetypes")
        expect_matrix_dim(fit[["coordinates"]], 3L, 2L)
        expect_matrix_dim(fit[["coefficients"]], 3L, 250L)
        expect_matrix_dim(fit[["compositions"]], 250L, 3L)
        expect_row_stochastic(fit[["coefficients"]])
        expect_row_stochastic(fit[["compositions"]])
    }
})

test_that("common fitter defaults are synchronized", {
    pgd_formals <- formals(archetypes_pgd)
    nnls_formals <- formals(archetypes_nnls)

    expect_identical(pgd_formals[["sd_threshold"]], nnls_formals[["sd_threshold"]])
    expect_identical(pgd_formals[["max_iter"]], nnls_formals[["max_iter"]])
    expect_identical(formals(run_aa)[["sd_threshold"]], nnls_formals[["sd_threshold"]])
    expect_identical(formals(run_aa)[["max_iter"]], nnls_formals[["max_iter"]])
})

test_that("fitters preserve the user-facing call", {
    X <- toy_matrix()

    pgd <- suppressWarnings(archetypes_pgd(X, K = 3L, max_iter = 1L))
    nnls <- suppressWarnings(archetypes_nnls(X, K = 3L, max_iter = 1L))
    entry <- suppressWarnings(run_aa(X, K = 3L, method = "nnls", max_iter = 1L))

    expect_identical(as.character(pgd[["call"]][[1L]]), "archetypes_pgd")
    expect_identical(as.character(nnls[["call"]][[1L]]), "archetypes_nnls")
    expect_identical(as.character(entry[["call"]][[1L]]), "run_aa")
})

test_that("run_aa validates method and method-specific arguments", {
    X <- toy_matrix()

    expect_error(run_aa(X, K = 3L, method = "bad"), "should be one of")
    expect_error(run_aa(X, K = 3L, method = "pgd", step_size = 0), "step_size")
    expect_error(run_aa(X, K = 3L, method = "nnls", bigM = 0), "bigM")
})

test_that("sparse preprocessing preserves sparse structure without centering", {
    X_dense <- matrix(
        c(
            1, 0, 2, 0,
            0, 3, 0, 1,
            4, 0, 0, 0,
            0, 5, 1, 0,
            2, 0, 0, 3,
            0, 1, 0, 0
        ),
        nrow = 6L,
        byrow = TRUE,
        dimnames = list(paste0("x", 1:6), paste0("v", 1:4))
    )
    X_sparse <- Matrix::Matrix(X_dense, sparse = TRUE)

    pre <- .aa_preprocess(X_sparse, sd_threshold = 0, weights = NULL, verbose = FALSE)

    expect_s4_class(pre[["X"]], "sparseMatrix")
    expect_null(attr(pre[["X"]], "scaled:center"))
    expect_equal(attr(pre[["X"]], "scaled:scale"), apply(X_dense, 2L, stats::sd))
    expect_equal(.dist2(X_sparse, center = TRUE), .dist2(X_dense, center = TRUE))
})

test_that("sparse preprocessing keeps NNLS bigM column sparse", {
    X <- Matrix::Matrix(
        matrix(c(1, 0, 0, 2, 3, 0, 0, 4), nrow = 4L, byrow = TRUE),
        sparse = TRUE
    )

    pre <- .aa_preprocess(X, sd_threshold = 0, weights = NULL, verbose = FALSE, bigM = 200)

    expect_s4_class(pre[["X"]], "sparseMatrix")
    expect_equal(attr(pre[["X"]], "bigM"), 1L)
    expect_equal(as.numeric(pre[["X"]][, 1L]), rep(200, nrow(X)))
})

test_that("archetypes fitters accept sparse input with expected invariants", {
    X <- Matrix::Matrix(
        matrix(
            c(
                1, 0, 2, 0,
                0, 3, 0, 1,
                4, 0, 0, 0,
                0, 5, 1, 0,
                2, 0, 0, 3,
                0, 1, 0, 0
            ),
            nrow = 6L,
            byrow = TRUE,
            dimnames = list(paste0("x", 1:6), paste0("v", 1:4))
        ),
        sparse = TRUE
    )

    set.seed(2)
    pgd <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "uniform_archetypes",
        sd_threshold = 0,
        max_iter = 2L
    ))
    set.seed(2)
    nnls <- suppressWarnings(archetypes_nnls(
        X,
        K = 2L,
        init = "uniform_archetypes",
        sd_threshold = 0,
        max_iter = 1L
    ))

    for (fit in list(pgd, nnls)) {
        expect_s3_class(fit, "archetypes")
        expect_matrix_dim(fit[["coordinates"]], 2L, 4L)
        expect_equal(dim(fit[["coefficients"]]), c(2L, 6L))
        expect_equal(dim(fit[["compositions"]]), c(6L, 2L))
        expect_row_stochastic(fit[["coefficients"]])
        expect_row_stochastic(fit[["compositions"]])
        expect_true(all(is.finite(fit[["loss"]][["rss"]])))
    }
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
