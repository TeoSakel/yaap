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

test_that("fitters integrate with hull_outmost initialization", {
    set.seed(1)
    X <- toy_matrix()

    fits <- list(
        suppressWarnings(run_aa(
            X,
            K = 3L,
            method = "pgd",
            init = "hull_outmost",
            init_args = list(hull_method = "projected", projected_dim = 2L),
            max_iter = 5L,
            tol_r2 = 0.9
        )),
        suppressWarnings(run_aa(
            X,
            K = 3L,
            method = "nnls",
            init = "hull_outmost",
            init_args = list(hull_method = "partitioned", n_partitions = 5L),
            max_iter = 3L,
            tol_r2 = 0.9
        )),
        suppressWarnings(archetypes_pgd(
            X,
            K = 3L,
            init = "hull_outmost",
            init_args = list(hull_method = "full"),
            max_iter = 5L,
            tol_r2 = 0.9
        )),
        suppressWarnings(archetypes_nnls(
            X,
            K = 3L,
            init = "hull_outmost",
            init_args = list(hull_method = "projected", projected_dim = 1L),
            max_iter = 3L,
            tol_r2 = 0.9
        ))
    )

    for (fit in fits) {
        expect_s3_class(fit, "archetypes")
        expect_matrix_dim(fit[["coordinates"]], 3L, 2L)
        expect_matrix_dim(fit[["coefficients"]], 3L, nrow(X))
        expect_matrix_dim(fit[["compositions"]], nrow(X), 3L)
        expect_row_stochastic(fit[["coefficients"]])
        expect_row_stochastic(fit[["compositions"]])
    }
})

test_that("common fitter defaults are synchronized", {
    pgd_formals <- formals(archetypes_pgd)
    nnls_formals <- formals(archetypes_nnls)

    expect_identical(pgd_formals[["sd_threshold"]], nnls_formals[["sd_threshold"]])
    expect_identical(pgd_formals[["max_iter"]], nnls_formals[["max_iter"]])
    expect_identical(formals(run_aa.default)[["sd_threshold"]], nnls_formals[["sd_threshold"]])
    expect_identical(formals(run_aa.default)[["max_iter"]], nnls_formals[["max_iter"]])
})

test_that("run_aa is an S3 generic", {
    expect_true(isS3stdGeneric(run_aa))
    expect_true("run_aa.default" %in% methods("run_aa"))
    expect_true("run_aa.fd" %in% methods("run_aa"))
})

test_that("run_aa handles fda fd objects through basis coefficients", {
    testthat::skip_if_not_installed("fda")

    set.seed(1)
    basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = 5L)
    coefs <- matrix(
        seq_len(25L) / 25,
        nrow = 5L,
        ncol = 5L,
        dimnames = list(paste0("b", 1:5), paste0("x", 1:5))
    )
    fd <- fda::fd(coefs, basis)

    fit <- suppressWarnings(run_aa(fd, K = 2L, max_iter = 1L, sd_threshold = 0))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["coordinates"]], 2L, 5L)
    expect_matrix_dim(fit[["init"]], 2L, 5L)
    expect_s3_class(fit[["data"]], "fd")
    expect_matrix_dim(fit[["coefficients"]], 2L, 5L)
    expect_matrix_dim(fit[["compositions"]], 5L, 2L)
    expect_equal(anames(fit), rownames(fit[["coefficients"]]))
    expect_s3_class(fitted(fit), "fd")
    expect_s3_class(residuals(fit), "fd")
    expect_matrix_dim(predict(fit, fd), 5L, 2L)

    A_fd <- coordinates_fd(fit)
    expect_s3_class(A_fd, "fd")
    expect_equal(dim(stats::coef(A_fd)), c(5L, 2L))

    anames(fit) <- c("left", "right")
    expect_equal(anames(fit), c("left", "right"))
    expect_equal(colnames(stats::coef(coordinates_fd(fit))), c("left", "right"))
    expect_s3_class(coordinates_fd(fit, basis = basis), "fd")
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
    expect_error(run_aa(X, K = 3L, method = "pgd", max_no_update = 0), "max_no_update")
    expect_error(run_aa(X, K = 3L, method = "nnls", bigM = 0), "bigM")
    expect_error(run_aa(X, K = 3L, method = "nnls", max_no_update = 0), "max_no_update")
})

test_that("PGD stalls instead of converging when no updates are accepted", {
    set.seed(1)
    X <- toy_matrix()

    expect_warning(
        withCallingHandlers(
            fit <- run_aa(
                X,
                K = 3L,
                method = "pgd",
                max_iter = 10L,
                step_size = 1e100,
                max_iter_optimizer = 1L,
                max_no_update = 2L
            ),
            warning = function(w) {
                if (grepl("Algorithm did not converge", conditionMessage(w)))
                    invokeRestart("muffleWarning")
            }
        ),
        "PGD stalled"
    )

    expect_false(fit[["converged"]])
    expect_equal(nrow(fit[["loss"]]), 3L)
    expect_equal(diff(fit[["loss"]][["loss"]]), c(0, 0))
})

test_that("scale preprocessing supports TRUE, FALSE, vector, and matrix transforms", {
    X <- toy_matrix()
    pre_default <- .aa_preprocess(X, sd_threshold = 0, weights = NULL, verbose = FALSE)
    pre_raw <- .aa_preprocess(X, sd_threshold = 0, weights = NULL, verbose = FALSE, scale = FALSE)
    sd <- apply(X, 2L, stats::sd)
    pre_vector <- .aa_preprocess(
        X,
        sd_threshold = 0,
        weights = NULL,
        verbose = FALSE,
        scale = sd
    )
    pre_matrix <- .aa_preprocess(
        X,
        sd_threshold = 0,
        weights = NULL,
        verbose = FALSE,
        scale = diag(1 / sd^2)
    )

    expect_equal(unname(colMeans(pre_default[["X"]])), rep(0, ncol(X)), tolerance = 1e-12)
    expect_equal(unname(apply(pre_default[["X"]], 2L, stats::sd)), rep(1, ncol(X)))
    expect_equal(pre_raw[["X"]], X, ignore_attr = TRUE)
    expect_equal(
        as.matrix(dist(pre_default[["X"]]))^2,
        as.matrix(dist(pre_matrix[["X"]]))^2,
        tolerance = 1e-10
    )
    expect_equal(pre_vector[["X"]], pre_matrix[["X"]], ignore_attr = TRUE)
    expect_equal(
        pre_vector[["undo_scale"]](pre_vector[["X"]][1:2, ], pre_vector[["X"]]),
        pre_matrix[["undo_scale"]](pre_matrix[["X"]][1:2, ], pre_matrix[["X"]])
    )
    expect_equal(pre_raw[["undo_scale"]](pre_raw[["X"]][1:2, ], pre_raw[["X"]]), X[1:2, ])
})

test_that("automatic bigM preserves old default on z-scored data and scales raw data", {
    X <- toy_matrix()

    pre_default <- .aa_preprocess(X, sd_threshold = 0, weights = NULL, verbose = FALSE, bigM = NULL)
    pre_raw <- .aa_preprocess(
        10 * X,
        sd_threshold = 0,
        weights = NULL,
        verbose = FALSE,
        bigM = NULL,
        scale = FALSE
    )

    expect_equal(attr(pre_default[["X"]], "bigM.value"), 200)
    expect_gt(attr(pre_raw[["X"]], "bigM.value"), 200)
})

test_that("scale validation rejects invalid inputs", {
    X <- toy_matrix()

    expect_error(archetypes_pgd(X, K = 3L, scale = diag(3)), "one row and column per feature")
    expect_error(archetypes_pgd(X, K = 3L, scale = matrix(c(1, 2, 0, 1), 2)), "symmetric")
    bad_scale <- diag(2)
    bad_scale[1, 1] <- NA_real_
    expect_error(archetypes_pgd(X, K = 3L, scale = bad_scale), "missing or non-finite")
    expect_error(archetypes_pgd(X, K = 3L, scale = c(1, 1, 1)), "one value per feature")
    expect_error(archetypes_pgd(X, K = 3L, scale = c(1, 0)), "positive")
    expect_error(archetypes_pgd(X, K = 3L, scale = diag(c(1, 0))), "positive definite")
})

test_that("PGD and NNLS accept matrix scale and return original-unit coordinates", {
    set.seed(5)
    X <- toy_matrix()
    scale <- matrix(c(2, 0.3, 0.3, 1), nrow = 2L)

    pgd <- suppressWarnings(archetypes_pgd(X, K = 3L, scale = scale, max_iter = 2L))
    nnls <- suppressWarnings(archetypes_nnls(X, K = 3L, scale = scale, max_iter = 2L))

    for (fit in list(pgd, nnls)) {
        expect_s3_class(fit, "archetypes")
        expect_matrix_dim(fit[["coordinates"]], 3L, ncol(X))
        expect_equal(colnames(fit[["coordinates"]]), colnames(X))
        expect_true(all(is.finite(fit[["coordinates"]])))
        expect_true(all(is.finite(fit[["loss"]][["loss"]])))
    }
})

test_that("NNLS matrix scale keeps bigM outside returned coordinates", {
    set.seed(6)
    X <- toy_matrix()
    scale <- Matrix::Diagonal(n = ncol(X), x = c(2, 1))

    fit <- suppressWarnings(archetypes_nnls(X, K = 3L, scale = scale, max_iter = 1L, bigM = 200))

    expect_matrix_dim(fit[["coordinates"]], 3L, ncol(X))
    expect_false("bigM" %in% colnames(fit[["coordinates"]]))
    expect_true(all(is.finite(fit[["loss"]][["loss"]])))
})

test_that("NNLS keeps previous iterate when candidate loss does not improve", {
    set.seed(1)
    X <- toy_matrix()

    fit <- suppressWarnings(archetypes_nnls(
        X,
        K = 3L,
        max_iter = 20L,
        bigM = 5,
        max_no_update = 2L,
        init_args = list(refinement_steps = 0L)
    ))
    loss <- fit[["loss"]]

    expect_false(fit[["converged"]])
    expect_true(all(diff(loss[["loss"]]) <= 0))
    expect_equal(tail(diff(loss[["loss"]]), 2L), c(0, 0))
})

test_that("NNLS warns when raw coefficients are far from simplex", {
    set.seed(1)
    X <- toy_matrix()

    expect_warning(
        withCallingHandlers(
            archetypes_nnls(X, K = 3L, max_iter = 1L, bigM = 1),
            warning = function(w) {
                if (grepl("Algorithm did not converge", conditionMessage(w)))
                    invokeRestart("muffleWarning")
            }
        ),
        "not close to simplex"
    )
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
        expect_true(all(is.finite(fit[["loss"]][["loss"]])))
    }
})

test_that("missing-data PGD defaults on for dense NA input", {
    X <- matrix(
        c(
            1, NA, 2, 0,
            0, 3, NA, 1,
            4, 0, 0, NA,
            NA, 5, 1, 0,
            2, NA, 0, 3,
            0, 1, NA, 0
        ),
        nrow = 6L,
        byrow = TRUE,
        dimnames = list(paste0("x", 1:6), paste0("v", 1:4))
    )

    set.seed(3)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "uniform_archetypes",
        sd_threshold = 0,
        max_iter = 3L
    ))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["coordinates"]], 2L, 4L)
    expect_equal(dim(fit[["coefficients"]]), c(2L, 6L))
    expect_equal(dim(fit[["compositions"]]), c(6L, 2L))
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_true(all(is.finite(fit[["loss"]][["loss"]])))
    expect_true(all(diff(fit[["loss"]][["loss"]]) <= 1e-8))
})

test_that("missing-data PGD treats sparse structural zeros as missing", {
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
    X <- Matrix::Matrix(X_dense, sparse = TRUE)

    set.seed(4)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "uniform_archetypes",
        missing = TRUE,
        sd_threshold = 0,
        max_iter = 3L
    ))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["coordinates"]], 2L, 4L)
    expect_equal(dim(fit[["coefficients"]]), c(2L, 6L))
    expect_equal(dim(fit[["compositions"]]), c(6L, 2L))
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_true(all(is.finite(fit[["loss"]][["loss"]])))
    expect_true(all(diff(fit[["loss"]][["loss"]]) <= 1e-8))
})

test_that("missing preprocessing scales observed entries", {
    X <- matrix(
        c(
            1, NA,
            2, 10,
            NA, 20,
            4, 30
        ),
        nrow = 4L,
        byrow = TRUE,
        dimnames = list(NULL, c("a", "b"))
    )
    pre <- .aa_preprocess(X, sd_threshold = 0, weights = NULL, verbose = FALSE, missing = TRUE)
    M <- pre[["M"]]

    expect_equal(unname(colMeans(pre[["X"]][M[, 1L], 1L, drop = FALSE])), 0, tolerance = 1e-12)
    expect_equal(stats::sd(pre[["X"]][M[, 1L], 1L]), 1)
    expect_equal(unname(colMeans(pre[["X"]][M[, 2L], 2L, drop = FALSE])), 0, tolerance = 1e-12)
    expect_equal(stats::sd(pre[["X"]][M[, 2L], 2L]), 1)

    X_sparse <- Matrix::Matrix(
        matrix(c(1, 0, 2, 10, 0, 20, 4, 30), nrow = 4L, byrow = TRUE),
        sparse = TRUE
    )
    pre_sparse <- .aa_preprocess(
        X_sparse,
        sd_threshold = 0,
        weights = NULL,
        verbose = FALSE,
        missing = TRUE
    )
    entries <- Matrix::summary(pre_sparse[["M"]])
    observed <- pre_sparse[["X"]][pre_sparse[["M"]]]

    expect_s4_class(pre_sparse[["M"]], "sparseMatrix")
    expect_equal(length(observed), length(entries[["i"]]))
    expect_true(all(is.finite(observed)))
})

test_that("missing preprocessing sparsifies very sparse dense masks", {
    X <- matrix(NA_real_, nrow = 10L, ncol = 10L)
    X[cbind(seq_len(9L), seq_len(9L))] <- seq_len(9L)

    pre <- .aa_preprocess(
        X,
        sd_threshold = 0,
        weights = NULL,
        verbose = FALSE,
        scale = FALSE,
        missing = TRUE
    )

    expect_s4_class(pre[["X"]], "sparseMatrix")
    expect_s4_class(pre[["M"]], "sparseMatrix")
    expect_equal(length(Matrix::summary(pre[["M"]])[["i"]]), 9L)
})

test_that("missing-data PGD validates unsupported combinations", {
    X <- toy_matrix()
    X[1, 1] <- NA_real_

    expect_error(archetypes_pgd(X, K = 3L, robust = TRUE), "robust")
    expect_error(archetypes_pgd(X, K = 3L, weights = rep(1, nrow(X))), "weights")
    expect_error(archetypes_pgd(X, K = 3L, scale = diag(ncol(X))), "matrix `scale`")
    expect_error(run_aa(X, K = 3L, method = "nnls"), "missing = TRUE")
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
        expect_true(all(is.finite(fit[["loss"]][["loss"]])))
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
