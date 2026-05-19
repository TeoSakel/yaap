paa_binary_data <- function() {
    matrix(
        c(
            1, 0, 0, 1,
            1, 1, 0, 1,
            0, 1, 1, 0,
            0, 0, 1, 0,
            1, 0, 1, 1,
            0, 1, 0, 0
        ),
        ncol = 4,
        byrow = TRUE
    )
}

paa_count_data <- function() {
    matrix(
        c(
            4, 0, 1, 0,
            5, 1, 0, 0,
            0, 3, 5, 1,
            1, 4, 4, 0,
            0, 1, 0, 6,
            1, 0, 1, 5
        ),
        ncol = 4,
        byrow = TRUE
    )
}

test_that("archetypes_paa fits supported families with shared invariants", {
    cases <- list(
        gaussian = toy_matrix(),
        bernoulli = paa_binary_data(),
        poisson = paa_count_data(),
        multinomial = paa_count_data()
    )

    for (family in names(cases)) {
        fit <- suppressWarnings(archetypes_paa(
            cases[[family]],
            K = 3L,
            family = family,
            max_iter = 6L,
            tol_r2 = 1
        ))

        expect_s3_class(fit, "archetypes")
        expect_identical(fit[["family"]], family)
        expect_matrix_dim(coordinates(fit), 3L, ncol(cases[[family]]))
        expect_matrix_dim(fit[["coefficients"]], 3L, nrow(cases[[family]]))
        expect_matrix_dim(fit[["compositions"]], nrow(cases[[family]]), 3L)
        expect_row_stochastic(fit[["coefficients"]])
        expect_row_stochastic(fit[["compositions"]])
        expect_true(is_all_finite(fit[["loss"]][["loss"]]))
        expect_named(fit[["loss"]], c("loss", "r2"))
        expect_true(all(diff(fit[["loss"]][["loss"]]) <= 1e-8))
    }
})

test_that("run_aa dispatches to PAA", {
    fit <- suppressWarnings(run_aa(
        paa_binary_data(),
        K = 3L,
        method = "paa",
        family = "bernoulli",
        max_iter = 4L,
        tol_r2 = 1
    ))

    expect_s3_class(fit, "archetypes")
    expect_identical(fit[["family"]], "bernoulli")
    expect_true(is_all_finite(fit[["loss"]][["loss"]]))
})

test_that("PAA warns and ignores explicit run_aa scaling", {
    X <- paa_count_data()

    expect_warning(
        fit <- run_aa(X, K = 2L, method = "paa", family = "poisson", scale = TRUE),
        "ignored"
    )
    expect_s3_class(fit, "archetypes")

    expect_warning(
        fit <- run_aa(X, K = 2L, method = "paa", family = "poisson", scale = c(1, 1, 1, 1)),
        "ignored"
    )
    expect_s3_class(fit, "archetypes")

    df <- as.data.frame(X)
    expect_warning(
        fit <- run_aa(V1 ~ ., data = df, K = 2L, method = "paa", family = "poisson", scale = TRUE),
        "ignored"
    )
    expect_s3_class(fit, "archetypes")
})

test_that("PAA validates family-specific inputs", {
    expect_error(archetypes_paa(matrix(c(0, 2, 1, 0), ncol = 2), 2L, family = "bernoulli"),
                 "Bernoulli")
    expect_error(archetypes_paa(matrix(c(0, -1, 1, 0), ncol = 2), 2L, family = "poisson"),
                 "Poisson")
    expect_error(archetypes_paa(matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE),
                                2L, family = "multinomial"),
                 "positive totals")
})

test_that("PAA fitted and predict are family-aware", {
    X <- paa_count_data()
    fit <- suppressWarnings(archetypes_paa(
        X,
        K = 3L,
        family = "multinomial",
        max_iter = 6L,
        tol_r2 = 1
    ))

    X_hat <- fitted(fit)
    expect_matrix_dim(X_hat, nrow(X), ncol(X))
    expect_equal(unname(rowSums(X_hat)), unname(rowSums(X)), tolerance = 1e-6)

    rec <- predict(fit, X[1:3, , drop = FALSE], max_iter = 3L)
    rec_explicit <- predict(fit, X[1:3, , drop = FALSE], type = "reconstruction", max_iter = 3L)
    pred <- predict(fit, X[1:3, , drop = FALSE], type = "compositions", max_iter = 3L)
    expect_matrix_dim(pred, 3L, 3L)
    expect_row_stochastic(pred)
    expect_equal(rec, rec_explicit)
    expect_matrix_dim(rec, 3L, ncol(X))
    expect_equal(unname(rowSums(rec)), unname(rowSums(X[1:3, , drop = FALSE])), tolerance = 1e-6)
})

test_that("archetypes default family is gaussian", {
    fit <- manual_fit()
    expect_identical(fit[["family"]], "gaussian")
})

test_that("PAA coordinate plots reject observation-space data for non-Gaussian families", {
    fit <- suppressWarnings(archetypes_paa(
        paa_binary_data(),
        K = 3L,
        family = "bernoulli",
        max_iter = 3L,
        tol_r2 = 1
    ))

    expect_error(
        plot(fit, "coordinates"),
        "parameter space"
    )
})

test_that("PAA profile plots use parameter-space coordinates", {
    fit <- suppressWarnings(archetypes_paa(
        paa_binary_data(),
        K = 3L,
        family = "bernoulli",
        max_iter = 3L,
        tol_r2 = 1
    ))

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    prof <- plot(fit, "profiles")

    expect_named(prof, c("archetype", "feature", "value"))
    expect_equal(
        prof[["value"]],
        as.vector(coordinates(fit)),
        tolerance = 1e-8
    )
})
