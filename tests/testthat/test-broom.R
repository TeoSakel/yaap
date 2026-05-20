test_that("tidy.archetypes returns correct long-form for 'coordinates'", {
    fit <- manual_fit()
    X <- fit[["data"]]

    td <- tidy(fit) # default: matrix = "coordinates"

    expect_s3_class(td, "tbl_df")
    expect_named(td, c("archetype", "term", "value"))
    expect_equal(nrow(td), 3L * ncol(X))
    expect_true(all(td[["archetype"]] %in% anames(fit)))
    expect_true(all(td[["term"]] %in% colnames(X)))
    expect_type(td[["value"]], "double")
})

test_that("tidy.archetypes returns correct long-form for 'coefficients'", {
    fit <- manual_fit()
    X <- fit[["data"]]

    td <- tidy(fit, matrix = "coefficients")

    expect_s3_class(td, "tbl_df")
    expect_named(td, c("archetype", "sample", "value"))
    expect_equal(nrow(td), 3L * nrow(X))
    expect_true(all(td[["archetype"]] %in% anames(fit)))
})

test_that("tidy.archetypes returns correct long-form for 'compositions'", {
    fit <- manual_fit()
    X <- fit[["data"]]

    td <- tidy(fit, matrix = "compositions")

    expect_s3_class(td, "tbl_df")
    expect_named(td, c("sample", "archetype", "value"))
    expect_equal(nrow(td), nrow(X) * 3L)
    # reconstructed compositions should still sum to 1 per sample
    sums <- tapply(td[["value"]], td[["sample"]], sum)
    expect_equal(as.numeric(sums), rep(1, nrow(X)), tolerance = 1e-6)
})

test_that("tidy.archetypes uses anames as archetype labels", {
    fit <- manual_fit()
    anames(fit) <- c("Alpha", "Beta", "Gamma")

    td <- tidy(fit)
    expect_true(all(td[["archetype"]] %in% c("Alpha", "Beta", "Gamma")))
})

test_that("glance.archetypes returns a one-row tibble with correct fields", {
    fit <- manual_fit()

    gl <- glance(fit)

    expect_s3_class(gl, "tbl_df")
    expect_equal(nrow(gl), 1L)
    expect_named(gl, c("K", "converged", "loss", "r2", "n_iter", "aic", "family"))
    expect_equal(gl[["K"]], 3L)
    expect_type(gl[["converged"]], "logical")
    expect_type(gl[["loss"]], "double")
    expect_true(gl[["r2"]] >= 0 && gl[["r2"]] <= 1)
    expect_true(gl[["n_iter"]] >= 0L)
    expect_equal(gl[["family"]], "gaussian")
    # AIC requires stored data; fit above has data stored so should be finite
    expect_true(is.numeric(gl[["aic"]]))
})

test_that("augment.archetypes uses stored data when data = NULL", {
    fit <- manual_fit()
    X <- fit[["data"]]

    aug <- augment(fit)

    expect_s3_class(aug, "tbl_df")
    expect_equal(nrow(aug), nrow(X))
    comp_cols <- paste0(".", anames(fit))
    expect_true(all(comp_cols %in% colnames(aug)))
    expect_equal(ncol(aug), ncol(X) + 3L)
})

test_that("augment.archetypes errors when no data available", {
    X <- toy_matrix()
    fit <- archetypes(
        A            = matrix(1:6, 3, 2),
        coefficients = matrix(1 / 3, 3, 3),
        compositions = matrix(1 / 3, 3, 3),
        loss         = data.frame(loss = 0, r2 = 1, k_S = 1, k_A = 1),
        feature_map  = .aa_identity_feature_map(matrix(1:6, 3, 2))
    )
    expect_error(augment(fit), "Data must be provided")
})

test_that("augment.archetypes accepts explicit new data", {
    fit <- manual_fit()
    X <- fit[["data"]]

    X_new <- X[2:4, ]
    aug <- augment(fit, data = X_new)

    expect_equal(nrow(aug), 3L)
    comp_cols <- paste0(".", anames(fit))
    expect_true(all(comp_cols %in% colnames(aug)))
    # compositions for new data must be non-negative and sum to 1
    S_new <- as.matrix(aug[, comp_cols])
    expect_true(all(S_new >= -1e-8))
    expect_equal(rowSums(S_new), rep(1, 3L), tolerance = 1e-6)
})

test_that("tidy/glance/augment dispatch without attaching broom", {
    # generics package re-exports the S3 generics; dispatch must work with
    # only yaap loaded (broom is in Suggests only)
    set.seed(1)
    X <- toy_matrix()
    fit <- suppressWarnings(run_aa(X, K = 2L, max_iter = 5L))
    expect_s3_class(tidy(fit), "tbl_df")
    expect_s3_class(glance(fit), "tbl_df")
    expect_s3_class(augment(fit), "tbl_df")
})

# kernel_archetypes -----------------------------------------------------------

test_that("tidy.kernel_archetypes works for coordinates", {
    set.seed(1)
    X <- toy_matrix()
    fit <- suppressWarnings(run_aa(X, K = 3L, method = "kernel", kernel = "linear", max_iter = 10L))

    td <- tidy(fit)
    expect_s3_class(td, "tbl_df")
    expect_named(td, c("archetype", "term", "value"))
    expect_equal(nrow(td), 3L * ncol(X))
})

test_that("tidy.kernel_archetypes warns and returns empty tibble when coordinates is NULL", {
    set.seed(1)
    X <- toy_matrix()
    G <- exp(-as.matrix(dist(X))^2 / median(as.matrix(dist(X)))^2)
    fit <- suppressWarnings(run_aa(G, K = 3L, method = "kernel", kernel = "precomputed", max_iter = 5L))

    expect_false("coordinates" %in% names(fit))
    expect_null(coordinates(fit))
    expect_warning(td <- tidy(fit), "coordinates")
    expect_equal(nrow(td), 0L)
})

test_that("tidy.kernel_archetypes returns correct compositions", {
    set.seed(1)
    X <- toy_matrix()
    fit <- suppressWarnings(run_aa(X, K = 3L, method = "kernel", kernel = "linear", max_iter = 10L))

    td <- tidy(fit, matrix = "compositions")
    expect_named(td, c("sample", "archetype", "value"))
    expect_equal(nrow(td), nrow(X) * 3L)
})

test_that("glance.kernel_archetypes returns a one-row tibble without family", {
    set.seed(1)
    X <- toy_matrix()
    fit <- suppressWarnings(run_aa(X, K = 3L, method = "kernel", kernel = "linear", max_iter = 10L))

    gl <- glance(fit)
    expect_s3_class(gl, "tbl_df")
    expect_equal(nrow(gl), 1L)
    expect_named(gl, c("K", "converged", "loss", "r2", "n_iter"))
    expect_false("family" %in% colnames(gl))
    expect_equal(gl[["K"]], 3L)
})

test_that("augment.kernel_archetypes uses stored compositions", {
    set.seed(1)
    X <- toy_matrix()
    fit <- suppressWarnings(run_aa(X, K = 3L, method = "kernel", kernel = "linear", max_iter = 10L))

    aug <- augment(fit)
    expect_s3_class(aug, "tbl_df")
    expect_equal(nrow(aug), nrow(X))
    comp_cols <- paste0(".", anames(fit))
    expect_true(all(comp_cols %in% colnames(aug)))
})
