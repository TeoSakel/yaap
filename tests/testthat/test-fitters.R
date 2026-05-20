test_that("Euclidean fitters and run_aa share core toy-data invariants", {
    X <- toy_matrix()
    cases <- list(
        pgd = function() archetypes_pgd(X, K = 3L, max_iter = 20L, tol_r2 = 0.95),
        nnls = function() archetypes_nnls(X, K = 3L, max_iter = 5L, tol_r2 = 0.95),
        run_aa_pgd = function() run_aa(X, K = 3L, method = "pgd", max_iter = 20L, tol_r2 = 0.95),
        run_aa_nnls = function() run_aa(X, K = 3L, method = "nnls", max_iter = 5L, tol_r2 = 0.95)
    )

    for (case in cases) {
        set.seed(1)
        fit <- suppressWarnings(case())
        expect_archetypes_fit(
            fit,
            K = 3L,
            n = nrow(X),
            p = ncol(X),
            fitted_dim = dim(X),
            residual_dim = dim(X)
        )
    }
})

test_that("PGD and NNLS integrate with hull_outmost initialization", {
    X <- toy_matrix()
    cases <- list(
        pgd = function() {
            run_aa(
                X,
                K = 3L,
                method = "pgd",
                init = "hull_outmost",
                init_args = list(hull_method = "projected", projected_dim = 2L),
                max_iter = 5L,
                tol_r2 = 0.9
            )
        },
        nnls = function() {
            run_aa(
                X,
                K = 3L,
                method = "nnls",
                init = "hull_outmost",
                init_args = list(hull_method = "partitioned", n_partitions = 5L),
                max_iter = 3L,
                tol_r2 = 0.9
            )
        }
    )

    for (case in cases) {
        set.seed(1)
        expect_archetypes_fit(suppressWarnings(case()), K = 3L, n = nrow(X), p = ncol(X))
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
    expect_true("run_aa.formula" %in% methods("run_aa"))
    expect_true("run_aa.fd" %in% methods("run_aa"))
})

test_that("run_aa handles formula input", {
    set.seed(1)
    fit <- suppressWarnings(run_aa(
        Species ~ .,
        data = iris,
        K = 3L,
        max_iter = 1L,
        sd_threshold = 0
    ))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(coordinates(fit), 3L, 4L)
    expect_matrix_dim(fit[["coefficients"]], 3L, nrow(iris))
    expect_matrix_dim(fit[["compositions"]], nrow(iris), 3L)
    expect_equal(colnames(fit[["data"]]), names(iris)[1:4])
    expect_equal(fit[["formula"]], Species ~ .)
    expect_identical(as.character(fit[["call"]][[1L]]), "run_aa")
})

test_that("run_aa formula input can use the formula environment", {
    set.seed(1)
    df <- iris[1:25, 1:4]
    fit <- with(df, suppressWarnings(run_aa(
        ~ Sepal.Length + Sepal.Width,
        K = 2L,
        max_iter = 1L
    )))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["data"]], nrow(df), 2L)
    expect_matrix_dim(coordinates(fit), 2L, 2L)
    expect_equal(colnames(fit[["data"]]), c("Sepal.Length", "Sepal.Width"))
})

test_that("run_aa formula input requires K", {
    expect_error(
        run_aa(Species ~ ., data = iris, max_iter = 1L),
        "argument \"K\" is missing"
    )
})

test_that("run_aa formula input handles missing data like model.frame", {
    df <- iris[1:25, ]
    df[1L, "Sepal.Length"] <- NA_real_
    fit <- suppressWarnings(run_aa(
        Species ~ .,
        data = df,
        K = 3L,
        max_iter = 1L,
        sd_threshold = 0
    ))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["data"]], nrow(df) - 1L, 4L)
    expect_matrix_dim(fit[["compositions"]], nrow(df) - 1L, 3L)
    expect_error(
        run_aa(Species ~ ., data = df, K = 3L, na.action = stats::na.fail),
        "missing values"
    )
})

test_that("run_aa formula input ignores the response", {
    df <- iris[1:25, ]
    df[1L, "Species"] <- NA
    fit <- suppressWarnings(run_aa(
        Species ~ .,
        data = df,
        K = 3L,
        max_iter = 1L,
        sd_threshold = 0
    ))

    expect_s3_class(fit, "archetypes")
    expect_matrix_dim(fit[["data"]], nrow(df), 4L)
    expect_equal(colnames(fit[["data"]]), names(iris)[1:4])
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
    expect_s3_class(coordinates(fit), "fd")
    expect_matrix_dim(t(stats::coef(coordinates(fit))), 2L, 5L)
    expect_matrix_dim(fit[["init"]], 2L, 5L)
    expect_s3_class(fit[["data"]], "fd")
    expect_matrix_dim(fit[["coefficients"]], 2L, 5L)
    expect_matrix_dim(fit[["compositions"]], 5L, 2L)
    expect_equal(anames(fit), rownames(fit[["coefficients"]]))
    expect_s3_class(fitted(fit), "fd")
    expect_s3_class(residuals(fit), "fd")
    rec_fd <- predict(fit, fd)
    expect_s3_class(rec_fd, "fd")
    expect_equal(predict(fit, fd, type = "compositions"), compositions(fit), tolerance = 1e-6, ignore_attr = TRUE)
    rec_fd_explicit <- predict(fit, fd, type = "reconstruction")
    expect_equal(stats::coef(rec_fd_explicit), stats::coef(rec_fd), tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(stats::coef(rec_fd), stats::coef(fitted(fit)), tolerance = 1e-6, ignore_attr = TRUE)

    A_fd <- coordinates(fit)
    expect_s3_class(A_fd, "fd")
    expect_equal(dim(stats::coef(A_fd)), c(5L, 2L))

    anames(fit) <- c("left", "right")
    expect_equal(anames(fit), c("left", "right"))
    expect_equal(colnames(stats::coef(coordinates(fit))), c("left", "right"))
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
    expect_error(run_aa(X, K = 3L, method = "nnls", max_kappa = 0), "max_kappa")
    expect_no_error(suppressWarnings(run_aa(X, K = 3L, method = "nnls", max_kappa = Inf, max_iter = 1L)))
    expect_error(run_aa(X, K = 3L, method = "pgd", max_kappa = Inf), "unused")
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
                if (grepl("Algorithm did not converge", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        "PGD stalled"
    )

    expect_equal(nrow(fit[["loss"]]), 3L)
    expect_equal(diff(fit[["loss"]][["loss"]]), c(0, 0))
})

.aa_test_euclidean_preprocess <- function(x,
                                          sd_threshold = 0,
                                          weights = NULL,
                                          verbose = FALSE,
                                          scale = TRUE,
                                          missing = FALSE,
                                          bigM = 0) {
    .aa_euclidean_preprocess(
        list(
            x = x,
            sd_threshold = sd_threshold,
            weights = weights,
            verbose = verbose,
            scale = scale,
            missing = missing
        ),
        bigM = bigM
    )
}

test_that("scale preprocessing supports TRUE, FALSE, vector, and matrix transforms", {
    X <- toy_matrix()
    pre_default <- .aa_test_euclidean_preprocess(X)
    pre_raw <- .aa_test_euclidean_preprocess(X, scale = FALSE)
    sd <- apply(X, 2L, stats::sd)
    pre_vector <- .aa_test_euclidean_preprocess(
        X,
        scale = sd
    )
    pre_matrix <- .aa_test_euclidean_preprocess(
        X,
        scale = diag(1 / sd^2)
    )

    matrix_scale_factor <- attr(pre_matrix[["X"]], "scale:factor")
    expect_type(matrix_scale_factor, "double")
    expect_equal(matrix_scale_factor[[1L]], ncol(pre_matrix[["X"]]))
    expect_equal(length(matrix_scale_factor) - 1L, ncol(pre_matrix[["X"]]) * (ncol(pre_matrix[["X"]]) + 1L) / 2)
    expect_equal(.aa_unpack_lower_tri(matrix_scale_factor), t(chol(diag(1 / sd^2))))
    expect_equal(unname(colMeans(pre_default[["X"]])), rep(0, ncol(X)), tolerance = 1e-12)
    expect_equal(unname(apply(pre_default[["X"]], 2L, stats::sd)), rep(1, ncol(X)))
    expect_equal(pre_raw[["X"]], X, ignore_attr = TRUE)
    expect_equal(
        as.matrix(dist(pre_default[["X"]]))^2,
        as.matrix(dist(pre_matrix[["X"]]))^2,
        tolerance = 1e-10
    )
    expect_equal(pre_vector[["X"]], pre_matrix[["X"]], ignore_attr = TRUE)
    vector_map <- .aa_euclidean_feature_map(pre_vector[["X"]])
    matrix_map <- .aa_euclidean_feature_map(pre_matrix[["X"]])
    raw_map <- .aa_euclidean_feature_map(pre_raw[["X"]])
    expect_equal(
        .aa_feature_map_inverse(vector_map, pre_vector[["X"]][1:2, ]),
        .aa_feature_map_inverse(matrix_map, pre_matrix[["X"]][1:2, ])
    )
    expect_equal(.aa_feature_map_inverse(raw_map, pre_raw[["X"]][1:2, ]), X[1:2, ])
})

test_that("euclidean preprocessing applies sqrt-normalized sample weights", {
    X <- toy_matrix()
    weights <- rep(1, nrow(X))
    weights[2] <- 4

    pre <- .aa_test_euclidean_preprocess(X, scale = FALSE, weights = weights)
    w_norm <- weights / mean(weights)
    X_expected <- X * sqrt(w_norm)

    expect_equal(pre[["X"]], X_expected, ignore_attr = TRUE)
    expect_equal(attr(pre[["X"]], "weights"), w_norm)
})

test_that("euclidean preprocessing warns when sample weights include zeros", {
    X <- toy_matrix()
    weights <- rep(1, nrow(X))
    weights[1] <- 0

    pre <- NULL
    expect_warning(
        pre <- .aa_test_euclidean_preprocess(X, scale = FALSE, weights = weights),
        "Some sample weights are zero"
    )

    w_norm <- weights / mean(weights)
    X_expected <- X * sqrt(w_norm)
    expect_equal(pre[["X"]], X_expected, ignore_attr = TRUE)
})

test_that("euclidean preprocessing rejects all-zero sample weights", {
    X <- toy_matrix()
    expect_error(
        .aa_test_euclidean_preprocess(X, weights = rep(0, nrow(X))),
        "at least one weight must be positive"
    )
})

test_that("automatic bigM preserves old default on z-scored data and scales raw data", {
    X <- toy_matrix()

    pre_default <- .aa_test_euclidean_preprocess(X, bigM = NULL)
    pre_raw <- .aa_test_euclidean_preprocess(
        10 * X,
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
        expect_archetypes_fit(fit, K = 3L, n = nrow(X), p = ncol(X))
        expect_equal(colnames(coordinates(fit)), colnames(X))
    }
})

test_that("NNLS matrix scale keeps bigM outside returned coordinates", {
    set.seed(6)
    X <- toy_matrix()
    scale <- Matrix::Diagonal(n = ncol(X), x = c(2, 1))

    fit <- suppressWarnings(archetypes_nnls(X, K = 3L, scale = scale, max_iter = 1L, bigM = 200))

    expect_matrix_dim(coordinates(fit), 3L, ncol(X))
    expect_false("bigM" %in% colnames(coordinates(fit)))
    expect_true(is_all_finite(fit[["loss"]][["loss"]]))
})

test_that("NNLS reports best loss when final candidate does not improve", {
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

    expect_named(loss, c("loss", "r2", "rloss", "k_S", "k_A"))
    expect_true(all(diff(loss[["loss"]]) <= 1e-8))
    expect_true(any(loss[["rloss"]] > loss[["loss"]]))
    expect_equal(loss[["loss"]][nrow(loss)], min(loss[["loss"]]))
})

test_that("NNLS max_no_update records rejected candidate before stalling", {
    set.seed(1)
    X <- toy_matrix()

    expect_warning(
        withCallingHandlers(
            fit <- archetypes_nnls(
                X,
                K = 3L,
                max_iter = 20L,
                bigM = 5,
                max_no_update = 1L,
                tol = 1e-12,
                tol_r2 = 1,
                init_args = list(refinement_steps = 0L)
            ),
            warning = function(w) {
                if (grepl("Algorithm did not converge", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        "NNLS stalled"
    )
    loss <- fit[["loss"]]

    expect_false(fit[["converged"]])
    expect_true(tail(loss[["rloss"]], 1L) > tail(loss[["loss"]], 1L))
})


test_that("PGD loss history contains only shared metrics", {
    set.seed(1)
    fit <- suppressWarnings(archetypes_pgd(
        toy_matrix(),
        K = 3L,
        max_iter = 2L,
        tol_r2 = 1
    ))

    expect_named(fit[["loss"]], c("loss", "r2"))
})

test_that(".aa_check_convergence only uses loss and r2", {
    loss <- data.frame(loss = c(10, 9), r2 = c(0.1, 0.2))

    expect_false(.aa_check_convergence(loss, 1L, tol = 1e-6, tol_r2 = 0.9, verbose = FALSE))
})

test_that("non-NNLS direct fitters do not accept max_kappa", {
    X <- toy_matrix()[1:12, , drop = FALSE]

    expect_error(archetypes_pgd(X, K = 3L, max_kappa = Inf), "unused")
    expect_error(archetypes_kernel_pgd(X, K = 3L, max_kappa = Inf), "unused")
    expect_error(archetypes_paa(X, K = 3L, max_kappa = Inf), "unused")
    expect_error(archetypes_directional(directional_matrix(12L), K = 3L, max_kappa = Inf), "unused")
})


test_that("NNLS warns when raw coefficients are far from simplex", {
    set.seed(1)
    X <- toy_matrix()

    expect_warning(
        withCallingHandlers(
            archetypes_nnls(X, K = 3L, max_iter = 1L, bigM = 1),
            warning = function(w) {
                if (grepl("Algorithm did not converge", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        "not close to simplex"
    )
})

test_that("non-convergence warning reports realized iteration count", {
    X <- matrix(seq_len(40), nrow = 10L, ncol = 4L)
    A <- matrix(seq_len(12), nrow = 3L, ncol = 4L)
    B <- matrix(1 / nrow(X), nrow = 3L, ncol = nrow(X))
    S <- matrix(1 / 3, nrow = nrow(X), ncol = 3L)
    loss <- data.frame(
        loss = seq(400, 369),
        r2 = seq(0.2, 0.7, length.out = 32L),
        k_S = NA_real_,
        k_A = NA_real_
    )

    ctx <- list(
        call = quote(run_aa(X, K = 3L, method = "nnls")),
        x = X,
        max_iter = 100L,
        verbose = FALSE
    )
    prep <- list(
        X = X,
        family = "gaussian"
    )
    fit <- list(
        A0 = A,
        A = A,
        B = B,
        S = S,
        delta = 0,
        i = 31L,
        loss = loss,
        converged = FALSE
    )

    expect_warning(
        out <- .aa_euclidean_output(ctx, prep, fit),
        "Algorithm did not converge after 31 iterations",
        fixed = TRUE
    )
    expect_equal(nrow(out[["loss"]]), 32L)
})

test_that("sparse preprocessing preserves sparse structure without centering", {
    X_sparse <- sparse_test_matrix()
    X_dense <- as.matrix(X_sparse)

    pre <- .aa_test_euclidean_preprocess(X_sparse)

    expect_s4_class(pre[["X"]], "sparseMatrix")
    expect_null(attr(pre[["X"]], "scaled:center"))
    expect_equal(attr(pre[["X"]], "scaled:scale"), apply(X_dense, 2L, stats::sd))
    expect_equal(.aa_dist2(X_sparse, center = TRUE), .aa_dist2(X_dense, center = TRUE))
})

test_that("sparse preprocessing keeps NNLS bigM column sparse", {
    X <- Matrix::Matrix(
        matrix(c(1, 0, 0, 2, 3, 0, 0, 4), nrow = 4L, byrow = TRUE),
        sparse = TRUE
    )

    pre <- .aa_test_euclidean_preprocess(X, bigM = 200)

    expect_s4_class(pre[["X"]], "sparseMatrix")
    expect_equal(attr(pre[["X"]], "bigM"), 1L)
    expect_equal(as.numeric(pre[["X"]][, 1L]), rep(200, nrow(X)))
})

test_that("archetypes fitters accept sparse input with expected invariants", {
    X <- sparse_test_matrix()

    set.seed(2)
    pgd <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "random",
        sd_threshold = 0,
        max_iter = 2L
    ))
    set.seed(2)
    nnls <- suppressWarnings(archetypes_nnls(
        X,
        K = 2L,
        init = "random",
        sd_threshold = 0,
        max_iter = 1L
    ))

    for (fit in list(pgd, nnls)) {
        expect_archetypes_fit(fit, K = 2L, n = nrow(X), p = ncol(X))
    }
})

test_that("missing-data PGD defaults on for dense NA input", {
    X <- missing_test_matrix()

    set.seed(3)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "random",
        sd_threshold = 0,
        max_iter = 3L
    ))

    expect_archetypes_fit(fit, K = 2L, n = nrow(X), p = ncol(X))
    expect_true(all(diff(fit[["loss"]][["loss"]]) <= 1e-8))
})

test_that("missing-data PGD supports zero optimizer iterations", {
    X <- missing_test_matrix()

    set.seed(3)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "random",
        sd_threshold = 0,
        max_iter = 0L
    ))

    expect_s3_class(fit, "archetypes")
    expect_true(fit[["converged"]])
    expect_equal(nrow(fit[["loss"]]), 1L)
    expect_true(is_number(fit[["loss"]][["loss"]]))
})

test_that("missing-data PGD handles K edge cases directly", {
    X <- missing_test_matrix()

    fit_mean <- archetypes_pgd(
        X,
        K = 1L,
        scale = FALSE,
        sd_threshold = 0,
        max_iter = 0L
    )
    expect_equal(unname(coordinates(fit_mean)[1L, ]), unname(colMeans(X, na.rm = TRUE)))
    expect_equal(nrow(fit_mean[["loss"]]), 1L)
    expect_true(is_number(AIC(fit_mean)))

    fit_identity <- archetypes_pgd(
        X,
        K = nrow(X),
        scale = FALSE,
        sd_threshold = 0,
        max_iter = 0L
    )
    expect_equal(unname(crossprod(as.matrix(fit_identity[["coefficients"]]))), diag(nrow(X)), ignore_attr = TRUE)
    expect_equal(fit_identity[["loss"]][["loss"]], 0)
    expect_equal(nrow(fit_identity[["loss"]]), 1L)
})

test_that("missing-data PGD treats sparse structural zeros as missing", {
    X <- sparse_test_matrix()

    set.seed(4)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "random",
        missing = TRUE,
        sd_threshold = 0,
        max_iter = 3L
    ))

    expect_archetypes_fit(fit, K = 2L, n = nrow(X), p = ncol(X))
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
    pre <- .aa_test_euclidean_preprocess(X, missing = TRUE)
    M <- pre[["M"]]

    expect_equal(unname(colMeans(pre[["X"]][M[, 1L], 1L, drop = FALSE])), 0, tolerance = 1e-12)
    expect_equal(stats::sd(pre[["X"]][M[, 1L], 1L]), 1)
    expect_equal(unname(colMeans(pre[["X"]][M[, 2L], 2L, drop = FALSE])), 0, tolerance = 1e-12)
    expect_equal(stats::sd(pre[["X"]][M[, 2L], 2L]), 1)

    X_sparse <- Matrix::Matrix(
        matrix(c(1, 0, 2, 10, 0, 20, 4, 30), nrow = 4L, byrow = TRUE),
        sparse = TRUE
    )
    pre_sparse <- .aa_test_euclidean_preprocess(
        X_sparse,
        missing = TRUE
    )
    entries <- Matrix::summary(pre_sparse[["M"]])
    observed <- pre_sparse[["X"]][pre_sparse[["M"]]]

    expect_s4_class(pre_sparse[["M"]], "sparseMatrix")
    expect_equal(length(observed), length(entries[["i"]]))
    expect_true(is_all_finite(observed))
})

test_that("missing preprocessing sparsifies very sparse dense masks", {
    X <- matrix(NA_real_, nrow = 10L, ncol = 10L)
    X[cbind(seq_len(9L), seq_len(9L))] <- seq_len(9L)

    pre <- .aa_test_euclidean_preprocess(
        X,
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
    expect_error(archetypes_pgd(X, K = 3L, robust = "psi.huber"), "robust")
    expect_error(archetypes_pgd(X, K = 3L, weights = rep(1, nrow(X))), "weights")
    expect_error(archetypes_pgd(X, K = 3L, scale = diag(ncol(X))), "matrix `scale`")
    expect_error(run_aa(X, K = 3L, method = "nnls"), "missing = TRUE")
})

test_that("robust fitters accept MASS psi selectors", {
    testthat::skip_if_not_installed("MASS")
    set.seed(1)
    X <- toy_matrix()
    X[1, ] <- X[1, ] + 25

    pgd <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        robust = "psi.huber",
        robust_args = list(k = 1.345),
        max_iter = 5L,
        tol_r2 = 0.95
    ))
    nnls <- suppressWarnings(archetypes_nnls(
        X,
        K = 3L,
        robust = "psi.hampel",
        robust_args = list(a = 2, b = 4, c = 8),
        max_iter = 3L,
        tol_r2 = 0.95
    ))
    kernel <- suppressWarnings(archetypes_kernel_pgd(
        X,
        K = 3L,
        kernel = "linear",
        robust = "psi.huber",
        robust_args = list(k = 1.345),
        max_iter = 3L,
        tol_r2 = 0.95
    ))

    for (fit in list(pgd, nnls, kernel)) {
        expect_archetypes_fit(fit, K = 3L, n = nrow(X), p = ncol(X))
        expect_true(fit[["fit_info"]][["robust"]])
        expect_match(fit[["fit_info"]][["robust_psi"]], "psi\\.(huber|hampel)")
    }
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
        expect_archetypes_fit(fit, K = 3L, n = nrow(X), p = ncol(X))
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

    expected_names <- rownames(init)

    for (fit in list(pgd, nnls)) {
        fit_names <- rownames(coordinates(fit))
        expect_named(as.data.frame(coordinates(fit)), colnames(X))
        expect_setequal(fit_names, expected_names)
        expect_equal(rownames(fit[["init"]]), fit_names)
        expect_equal(rownames(fit[["coefficients"]]), fit_names)
        expect_equal(colnames(fit[["compositions"]]), fit_names)
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

# nrep -----------------------------------------------------------------

test_that("run_aa with nrep = 1 returns same structure as default", {
    set.seed(42)
    X <- toy_matrix()
    fit_default <- suppressWarnings(run_aa(X, K = 3L, max_iter = 10L))
    set.seed(42)
    fit_nrep1 <- suppressWarnings(run_aa(X, K = 3L, max_iter = 10L, nrep = 1L))

    expect_s3_class(fit_nrep1, "archetypes")
    expect_equal(coordinates(fit_default), coordinates(fit_nrep1))
})

test_that("run_aa with nrep > 1 returns the best sequential restart", {
    set.seed(7)
    X <- toy_matrix()
    final_loss <- function(fit) tail(fit[["loss"]][["loss"]], 1L)

    fit_nrep <- suppressWarnings(run_aa(X, K = 3L, max_iter = 5L, nrep = 3L))
    nrep_loss <- final_loss(fit_nrep)

    set.seed(7)
    single_losses <- suppressWarnings(replicate(3L, {
        final_loss(run_aa(X, K = 3L, max_iter = 5L))
    }))
    expect_archetypes_fit(fit_nrep, K = 3L, n = nrow(X), p = ncol(X))
    expect_equal(nrep_loss, min(single_losses), tolerance = 1e-8)
})

test_that("run_aa nrep validation rejects invalid values", {
    X <- toy_matrix()
    for (nrep in list(0L, -1L, 1.5)) {
        expect_error(run_aa(X, K = 3L, nrep = nrep), "`nrep`")
    }
})

test_that("run_aa nrep works with non-default methods", {
    X <- toy_matrix()
    cases <- list(
        nnls = function() run_aa(X, K = 3L, method = "nnls", max_iter = 5L, nrep = 2L),
        kernel = function() {
            run_aa(
                tcrossprod(X),
                K = 3L,
                method = "kernel",
                kernel = "precomputed",
                data = X,
                max_iter = 5L,
                nrep = 2L
            )
        },
        paa = function() run_aa(X, K = 3L, method = "paa", max_iter = 5L, nrep = 2L)
    )

    for (case in cases) {
        set.seed(3)
        fit <- suppressWarnings(case())
        expect_true(inherits(fit, "archetypes") || inherits(fit, "kernel_archetypes"))
    }
})
