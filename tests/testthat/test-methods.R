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
    fit_no_data <- fit
    fit_no_data[["data"]] <- NULL
    expect_error(residuals(fit_no_data), "Original data")
    expect_output(print(fit), "Archetypes Summary")
    expect_true(is.numeric(AIC(fit)))
})

test_that("predict returns row-stochastic compositions for new data", {
    set.seed(1)
    X <- toy_data()
    fit <- archetypes_nnls(as.matrix(X), K = 3L, max_iter = 5L, tol_r2 = 0.95)

    pred <- predict(fit, X[1:5, ])

    expect_matrix_dim(pred, 5L, 3L)
    expect_row_stochastic(pred)
})
