test_that("archetypes constructor validates dimensions and stochastic rows", {
    fit <- manual_fit()

    expect_s3_class(fit, "archetypes")
    expect_error(
        archetypes(
            fit[["coordinates"]],
            fit[["coefficients"]][-1, , drop = FALSE],
            fit[["compositions"]]
        ),
        "nrow\\(coefficients\\)"
    )
    bad_coefficients <- fit[["coefficients"]]
    bad_coefficients[1, 1] <- 2
    expect_error(
        archetypes(fit[["coordinates"]], bad_coefficients, fit[["compositions"]])
    )
    bad_compositions <- fit[["compositions"]]
    bad_compositions[1, 1] <- 2
    expect_error(
        archetypes(fit[["coordinates"]], fit[["coefficients"]], bad_compositions)
    )
})

test_that("S3 methods return expected values", {
    fit <- manual_fit()
    X <- fit[["data"]]

    expect_equal(anames(fit), rownames(fit[["coordinates"]]))
    expect_true(all(c("coordinates", "coefficients", "compositions") %in% names(fit)))
    expect_equal(coefficients(fit), fit[["coefficients"]])
    expect_equal(fitted(fit), fit[["compositions"]] %*% fit[["coordinates"]])
    expect_equal(residuals(fit), X - fitted(fit))
    expect_equal(residuals(fit, type = "response"), X - fitted(fit))
    expect_equal(residuals(fit, type = "pearson"), X - fitted(fit))
    fit_no_data <- fit
    fit_no_data[["data"]] <- NULL
    expect_error(residuals(fit_no_data), "Original data")
    expect_output(print(fit), "Archetypes Summary")
    expect_true(is.numeric(AIC(fit)))
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
        coordinates = matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2L, byrow = TRUE),
        coefficients = matrix(c(1, 0, 0, 0, 1, 0), nrow = 2L, byrow = TRUE),
        compositions = matrix(c(1, 0, 0, 1, 0.5, 0.5), nrow = 3L, byrow = TRUE),
        data = matrix(seq_len(12), nrow = 3L)
    )
    expect_warning(
        expect_true(is.na(AIC(high_dim))),
        "not larger than the number of features"
    )

    singular <- archetypes(
        coordinates = matrix(c(1, 2, 2, 4), nrow = 2L, byrow = TRUE),
        coefficients = matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2L, byrow = TRUE),
        compositions = matrix(c(1, 0, 0, 1, 0.5, 0.5, 0.25, 0.75), nrow = 4L, byrow = TRUE),
        data = cbind(seq_len(4), 2 * seq_len(4))
    )
    expect_warning(
        aic <- AIC(singular),
        "pseudo-inverse"
    )
    expect_true(is.finite(aic))
})

test_that("predict returns row-stochastic compositions for new data", {
    set.seed(1)
    X <- toy_data()
    fit <- archetypes_nnls(as.matrix(X), K = 3L, max_iter = 5L, tol_r2 = 0.95)

    pred <- predict(fit, X[1:5, ])
    pred_default <- predict(fit, X[1:5, ], type = "compositions")
    rec <- predict(fit, X[1:5, ], type = "reconstruction")

    expect_matrix_dim(pred, 5L, 3L)
    expect_row_stochastic(pred)
    expect_equal(pred, pred_default)
    expect_matrix_dim(rec, 5L, ncol(X))
    expect_equal(rec, pred %*% fit[["coordinates"]], tolerance = 1e-8)
})
