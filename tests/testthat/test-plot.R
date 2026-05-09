test_that("plot.archetypes smoke tests supported plot modes", {
    fit <- manual_fit()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_identical(plot(fit, "loss"), fit)
    expect_identical(plot(fit, "compositions"), fit)
    expect_identical(plot(fit, "composition"), fit)
    expect_identical(plot(fit, "composision"), fit)
    expect_identical(plot(fit, "profiles"), fit)
    expect_identical(plot(fit, "coordinates"), fit)
})

test_that("plot.archetypes handles higher-dimensional coordinate projections", {
    fit <- manual_fit()
    X <- cbind(fit[["data"]], z = c(0, 1, 1, 0.2), w = c(1, 0, 1, 0.4))
    fit[["data"]] <- X
    fit[["coordinates"]] <- X[1:3, , drop = FALSE]
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_identical(plot(fit, "coordinates"), fit)
    expect_identical(plot(fit, "coordinates", projection = "pca"), fit)
})
