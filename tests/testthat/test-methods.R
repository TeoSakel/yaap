test_that("archetypes constructor validates dimensions and stochastic rows", {
    fit <- manual_fit()

    expect_s3_class(fit, "archetypes")
    expect_s3_class(fit[["feature_map"]], "euclidean_feature_map")
    default_map_fit <- archetypes(
        A = coordinates(fit),
        coefficients = fit[["coefficients"]],
        compositions = fit[["compositions"]]
    )
    expect_s3_class(default_map_fit[["feature_map"]], "euclidean_feature_map")
    expect_false("archetypes" %in% names(default_map_fit[["feature_map"]]))
    expect_equal(default_map_fit[["A"]], coordinates(fit))
    expect_error(
        archetypes(
            A = coordinates(fit),
            coefficients = fit[["coefficients"]],
            compositions = fit[["compositions"]],
            feature_map = list()
        ),
        "aa_feature_map"
    )
    expect_error(
        archetypes(
            coordinates(fit),
            fit[["coefficients"]][-1, , drop = FALSE],
            fit[["compositions"]],
            feature_map = fit[["feature_map"]]
        ),
        "nrow\\(coefficients\\)"
    )
    bad_coefficients <- fit[["coefficients"]]
    bad_coefficients[1, 1] <- 2
    expect_error(
        archetypes(coordinates(fit), bad_coefficients, fit[["compositions"]], feature_map = fit[["feature_map"]])
    )
    bad_compositions <- fit[["compositions"]]
    bad_compositions[1, 1] <- 2
    expect_error(
        archetypes(coordinates(fit), fit[["coefficients"]], bad_compositions, feature_map = fit[["feature_map"]])
    )
})

test_that("S3 methods return expected values", {
    fit <- manual_fit()
    X <- fit[["data"]]

    expect_equal(anames(fit), rownames(coordinates(fit)))
    expect_true(all(c("A", "coefficients", "compositions") %in% names(fit)))
    expect_equal(coordinates(fit), coordinates(fit))
    expect_equal(compositions(fit), fit[["compositions"]])
    expect_equal(coefficients(fit), fit[["coefficients"]])
    expect_equal(fitted(fit), compositions(fit) %*% coordinates(fit))
    expect_equal(residuals(fit), X - fitted(fit))
    expect_equal(residuals(fit, type = "response"), X - fitted(fit))
    expect_equal(residuals(fit, type = "pearson"), X - fitted(fit))
    fit_no_data <- fit
    fit_no_data[["data"]] <- NULL
    expect_error(residuals(fit_no_data), "Original data")
    expect_output(print(fit), "Archetypal Analysis")
    expect_true(is.numeric(AIC(fit)))
})

test_that("print.archetypes reports compact fitting method labels", {
    X <- toy_matrix()[1:20, , drop = FALSE]

    set.seed(1)
    robust <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        robust = TRUE,
        delta = 0.2,
        max_iter = 1L
    ))
    robust_output <- paste(capture.output(print(robust)), collapse = "\n")
    expect_match(robust_output, "Robust PGD Archetypal Analysis")
    expect_false(grepl("Relaxed", robust_output, fixed = TRUE))
    expect_match(robust_output, "Number of Archetypes: 3")

    X_missing <- missing_test_matrix()
    set.seed(2)
    missing_fit <- suppressWarnings(archetypes_pgd(
        X_missing,
        K = 2L,
        init = "random",
        sd_threshold = 0,
        max_iter = 1L
    ))
    missing_output <- paste(capture.output(print(missing_fit)), collapse = "\n")
    expect_match(missing_output, "PGD Archetypal Analysis")
    expect_false(grepl("Missing-data PGD", missing_output, fixed = TRUE))

    paa <- suppressWarnings(archetypes_paa(
        matrix(c(1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE),
        K = 2L,
        family = "bernoulli",
        max_iter = 1L
    ))
    paa_output <- paste(capture.output(print(paa)), collapse = "\n")
    expect_match(paa_output, "PAA Archetypal Analysis \\(bernoulli family\\)")
})

test_that("summary.archetypes reports fit details, loss, and coordinates", {
    X <- toy_matrix()[1:20, , drop = FALSE]

    set.seed(1)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        robust = TRUE,
        delta = 0.2,
        max_iter = 1L
    ))

    smry <- summary(fit)
    expect_s3_class(smry, "summary.archetypes")
    expect_equal(smry[["fit_info"]][["method"]], "pgd")
    expect_true(smry[["fit_info"]][["robust"]])
    expect_equal(smry[["fit_info"]][["delta"]], 0.2)
    expect_equal(smry[["fit_info"]][["scaling"]], "none")
    expect_equal(smry[["final_loss"]], fit[["loss"]][nrow(fit[["loss"]]), , drop = FALSE])
    expect_equal(smry[["coordinates"]], coordinates(fit))

    smry_output <- paste(capture.output(print(smry)), collapse = "\n")
    expect_match(smry_output, "Fit Details:")
    expect_match(smry_output, "none")
    expect_match(smry_output, "Final Loss Metrics:")
    expect_match(smry_output, "Coordinates:")
})

test_that("summary.archetypes reports scaling mode", {
    X <- toy_matrix()[1:20, , drop = FALSE]

    set.seed(1)
    unscaled <- suppressWarnings(archetypes_pgd(X, K = 3L, max_iter = 1L))
    expect_equal(summary(unscaled)[["fit_info"]][["scaling"]], "none")

    set.seed(1)
    zscored <- suppressWarnings(archetypes_pgd(X, K = 3L, scale = TRUE, max_iter = 1L))
    expect_equal(summary(zscored)[["fit_info"]][["scaling"]], "z-score")

    set.seed(1)
    zscored_nnls <- suppressWarnings(archetypes_nnls(X, K = 3L, scale = TRUE, max_iter = 1L))
    expect_equal(summary(zscored_nnls)[["fit_info"]][["scaling"]], "z-score")

    set.seed(1)
    custom <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        scale = c(x = 1, y = 2),
        max_iter = 1L
    ))
    expect_equal(summary(custom)[["fit_info"]][["scaling"]], "custom")

    set.seed(1)
    metric <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        scale = diag(c(1, 2)),
        max_iter = 1L
    ))
    expect_equal(summary(metric)[["fit_info"]][["scaling"]], "metric")
})

test_that("summary.kernel_archetypes omits coordinates", {
    X <- toy_matrix()[1:12, , drop = FALSE]
    fit <- suppressWarnings(archetypes_kernel_pgd(
        X,
        K = 3L,
        kernel = "linear",
        init = c(1L, 2L, 3L),
        max_iter = 1L
    ))

    smry <- summary(fit)
    expect_s3_class(smry, "summary.archetypes")
    expect_equal(smry[["fit_info"]][["method"]], "kernel")
    expect_equal(smry[["coordinates"]], coordinates(fit))
    expect_true(any(grepl("Coordinates:", capture.output(print(smry)), fixed = TRUE)))
})

test_that("pearson residuals apply stored sample weights", {
    fit <- manual_fit()
    w <- c(1, 4, 1, 0)
    fit[["weights"]] <- w
    res <- residuals(fit, type = "response")
    expect_equal(residuals(fit, type = "pearson"), sqrt(w) * res)
})

test_that("adapted AIC returns NA outside covariance assumptions", {
    high_dim <- archetypes(
        A = matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2L, byrow = TRUE),
        coefficients = matrix(c(1, 0, 0, 0, 1, 0), nrow = 2L, byrow = TRUE),
        compositions = matrix(c(1, 0, 0, 1, 0.5, 0.5), nrow = 3L, byrow = TRUE),
        data = matrix(seq_len(12), nrow = 3L),
        feature_map = .aa_identity_feature_map(matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2L, byrow = TRUE))
    )
    expect_warning(
        expect_true(is.na(AIC(high_dim))),
        "not larger than the number of features"
    )

    singular <- archetypes(
        A = matrix(c(1, 2, 2, 4), nrow = 2L, byrow = TRUE),
        coefficients = matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2L, byrow = TRUE),
        compositions = matrix(c(1, 0, 0, 1, 0.5, 0.5, 0.25, 0.75), nrow = 4L, byrow = TRUE),
        data = cbind(seq_len(4), 2 * seq_len(4)),
        feature_map = .aa_identity_feature_map(matrix(c(1, 2, 2, 4), nrow = 2L, byrow = TRUE))
    )
    expect_warning(
        aic <- AIC(singular),
        "pseudo-inverse"
    )
    expect_true(is_number(aic))
})

test_that("predict returns row-stochastic compositions for new data", {
    set.seed(1)
    X <- toy_data()
    fit <- archetypes_nnls(as.matrix(X), K = 3L, max_iter = 5L, tol_r2 = 0.95)

    rec <- predict(fit, X[1:5, ])
    rec_explicit <- predict(fit, X[1:5, ], type = "reconstruction")
    pred <- predict(fit, X[1:5, ], type = "compositions")

    expect_matrix_dim(pred, 5L, 3L)
    expect_row_stochastic(pred)
    expect_equal(rec_explicit, rec)
    expect_matrix_dim(rec, 5L, ncol(X))
    expect_equal(rec, pred %*% coordinates(fit), tolerance = 1e-8)
})

test_that("Gaussian predict uses training feature map for compositions and reconstructions", {
    set.seed(2)
    X <- toy_matrix()
    X[, 2L] <- 100 * X[, 2L]
    rownames(X) <- paste0("x", seq_len(nrow(X)))

    cases <- list(
        pgd_zscore = function() archetypes_pgd(X, K = 3L, scale = TRUE, max_iter = 5L, tol_r2 = 0.99),
        nnls_zscore = function() archetypes_nnls(X, K = 3L, scale = TRUE, max_iter = 5L, tol_r2 = 0.99),
        pgd_vector = function() archetypes_pgd(X, K = 3L, scale = c(x = 1, y = 50), max_iter = 5L, tol_r2 = 0.99),
        nnls_vector = function() archetypes_nnls(X, K = 3L, scale = c(x = 1, y = 50), max_iter = 5L, tol_r2 = 0.99),
        pgd_metric = function() archetypes_pgd(X, K = 3L, scale = diag(c(1, 1 / 2500)), max_iter = 5L, tol_r2 = 0.99),
        nnls_metric = function() archetypes_nnls(X, K = 3L, scale = diag(c(1, 1 / 2500)), max_iter = 5L, tol_r2 = 0.99)
    )

    for (case in cases) {
        set.seed(2)
        fit <- suppressWarnings(case())
        rec <- predict(fit, X)
        pred <- predict(fit, X, type = "compositions")

        expect_equal(pred, fit[["compositions"]], tolerance = 1e-6, ignore_attr = TRUE)
        expect_equal(rec, fitted(fit), tolerance = 1e-6, ignore_attr = TRUE)
    }
})

test_that("Gaussian predict applies training low-variance mask", {
    X <- cbind(toy_matrix(), constant = 7)

    set.seed(3)
    fit <- suppressWarnings(archetypes_pgd(X, K = 3L, scale = TRUE, max_iter = 3L))
    rec <- predict(fit, X)
    pred <- predict(fit, X, type = "compositions")

    expect_equal(pred, fit[["compositions"]], tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(rec, fitted(fit), tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(unname(rec[, "constant"]), rep(7, nrow(X)))
})

test_that("Gaussian predict accepts sample-weighted fits", {
    X <- toy_matrix()
    weights <- seq_len(nrow(X))

    set.seed(4)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        weights = weights,
        max_iter = 3L
    ))

    pred <- predict(fit, X, type = "compositions")
    rec <- predict(fit, X)

    expect_true(fit[["fit_info"]][["sample_weights"]])
    expect_matrix_dim(pred, nrow(X), 3L)
    expect_row_stochastic(pred)
    expect_matrix_dim(rec, nrow(X), ncol(X))
    expect_true(all(is.finite(rec)))
})

test_that("Gaussian predict accepts robust fits", {
    X <- toy_matrix()

    set.seed(5)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        robust = TRUE,
        max_iter = 3L
    ))

    pred <- predict(fit, X, type = "compositions")
    rec <- predict(fit, X)

    expect_true(fit[["fit_info"]][["robust"]])
    expect_matrix_dim(pred, nrow(X), 3L)
    expect_row_stochastic(pred)
    expect_matrix_dim(rec, nrow(X), ncol(X))
    expect_true(all(is.finite(rec)))
})

test_that("Gaussian predict rejects missing-data PGD fits", {
    X <- missing_test_matrix()

    set.seed(5)
    fit <- suppressWarnings(archetypes_pgd(
        X,
        K = 2L,
        init = "random",
        sd_threshold = 0,
        max_iter = 1L
    ))

    expect_error(predict(fit, X), "missing-data")
})
