test_that("composition consistency is symmetric with unit diagonal", {
    fit <- manual_fit()

    expect_equal(
        .aa_nmi(fit[["compositions"]], fit[["compositions"]]),
        1,
        tolerance = 1e-8
    )
    expect_equal(
        consistency(fit, fit, what = "compositions"),
        1,
        tolerance = 1e-8
    )

    fit2 <- fit
    fit2[["compositions"]] <- fit2[["compositions"]][, c(2, 3, 1), drop = FALSE]
    fit2[["coefficients"]] <- fit2[["coefficients"]][c(2, 3, 1), , drop = FALSE]
    fit2[["A"]] <- coordinates(fit2)[c(2, 3, 1), , drop = FALSE]

    score12 <- consistency(fit, fit2, what = "compositions")
    score21 <- consistency(fit2, fit, what = "compositions")

    expect_true(is.finite(score12))
    expect_equal(score12, score21, tolerance = 1e-8)
    expect_error(.aa_nmi(fit[["compositions"]] * 2, fit[["compositions"]]), "row-stochastic")
})

test_that("coordinate consistency returns NA for decreasing K and uses mean pairwise distance denominator", {
    fit3 <- manual_fit()
    X <- fit3[["data"]]
    A4_base <- coordinates(fit3)
    A4_base[, 1L] <- A4_base[, 1L] + 0.1
    A4 <- rbind(A4_base, c(0.5, 0.5))
    B4 <- diag(4)
    S4 <- diag(4)
    loss <- data.frame(loss = 0, r2 = 1, k_S = 1, k_A = 1)
    fit4 <- archetypes(
        A = A4,
        coefficients = B4,
        compositions = S4,
        loss = loss,
        data = X,
        feature_map = .aa_identity_feature_map(A4)
    )

    expected_d2 <- mean(.aa_greedy_coordinate_d2(coordinates(fit3), coordinates(fit4)))
    expected <- 1 - expected_d2 / (2 * sum(matrixStats::colVars(X)))

    expect_equal(consistency(fit3, fit4, what = "coordinates"), expected)
    expect_true(is.na(consistency(fit4, fit3, what = "coordinates")))
})

test_that("coordinate matching greedily selects the closest remaining pair", {
    ax <- matrix(c(
        0, 0,
        1, 0
    ), ncol = 2L, byrow = TRUE)
    ay <- matrix(c(
        0.9, 0,
        0.2, 0
    ), ncol = 2L, byrow = TRUE)

    expect_equal(
        .aa_greedy_coordinate_d2(ax, ay),
        c(0.01, 0.04),
        tolerance = 1e-12
    )
})
