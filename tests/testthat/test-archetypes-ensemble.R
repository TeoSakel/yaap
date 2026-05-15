test_that("tune_archetypes expands K, replicates, and atomic tuning arguments", {
    set.seed(1)
    X <- toy_matrix()[1:10, , drop = FALSE]

    ens <- tune_archetypes(
        X,
        K = 2:3,
        nrep = 2L,
        delta = c(0, 0.05),
        max_iter = 2L,
        tol_r2 = 0,
        max_kappa = Inf,
        eval = function(fit) utils::tail(fit[["loss"]][["loss"]], 1L),
        eval_name = "score"
    )

    expect_s3_class(ens, "archetypes_ensemble")
    expect_equal(nrow(ens[["grid"]]), 8L)
    expect_equal(length(ens[["fits"]]), 8L)
    expect_equal(sort(unique(ens[["grid"]][["K"]])), 2:3)
    expect_equal(sort(unique(ens[["grid"]][["replicate"]])), 1:2)
    expect_equal(sort(unique(ens[["grid"]][["delta"]])), c(0, 0.05))
    expect_equal(ens[["prefer_metric"]], "score")
    expect_equal(ens[["prefer_direction"]], "minimize")
    expect_true(all(is.finite(ens[["metrics"]][["score"]])))
})

test_that("tune_archetypes expands multiple initialization methods", {
    set.seed(2)
    X <- toy_matrix()[1:10, , drop = FALSE]

    ens <- tune_archetypes(
        X,
        K = 2L,
        init = c("random", "furthest_sum"),
        max_iter = 1L,
        tol_r2 = 0,
        max_kappa = Inf
    )

    expect_equal(nrow(ens[["grid"]]), 2L)
    expect_equal(sort(ens[["grid"]][["init"]]), c("furthest_sum", "random"))
    expect_equal(ens[["prefer_metric"]], "AIC")
    expect_true(all(is.finite(ens[["metrics"]][["AIC"]])))
})

test_that("summary and best use existing ensemble metrics", {
    set.seed(3)
    X <- toy_matrix()[1:10, , drop = FALSE]
    calls <- 0L
    score <- function(fit) {
        calls <<- calls + 1L
        nrow(fit[["coordinates"]])
    }

    ens <- tune_archetypes(
        X,
        K = 2:3,
        eval = score,
        eval_name = "score",
        direction = "maximize",
        max_iter = 1L,
        tol_r2 = 0,
        max_kappa = Inf
    )
    expect_equal(calls, 2L)
    expect_equal(nrow(summary(ens)), 2L)

    selected <- best(ens)
    expect_equal(selected[["metrics"]][["score"]], 3)
    expect_equal(calls, 2L)

    selected_loss <- best(ens, metric = "loss", direction = "minimize")
    expect_s3_class(selected_loss[["fit"]], "archetypes")
    expect_error(best(ens, metric = "not_a_metric"), "column")
    expect_equal(calls, 2L)
})

test_that("composition consistency is symmetric with unit diagonal", {
    fit <- manual_fit()
    expect_equal(
        .nmi(fit[["compositions"]], fit[["compositions"]]),
        1,
        tolerance = 1e-8
    )
    expect_equal(
        consistency(fit, fit, what = "compositions"),
        1,
        tolerance = 1e-8
    )

    ens <- new_archetypes_ensemble(
        data = fit[["data"]],
        fits = list(a = fit, b = fit),
        grid = data.frame(model_id = c("a", "b"), K = 3L, replicate = 1L),
        metrics = data.frame(model_id = c("a", "b"), K = 3L, replicate = 1L, AIC = c(1, 2)),
        prefer_metric = "AIC",
        prefer_direction = "minimize",
        eval_fun = AIC,
        eval_name = "AIC"
    )

    score <- consistency(ens, "compositions")

    expect_true(all(is.finite(score)))
    expect_equal(score, t(score))
    expect_equal(unname(diag(score)), c(1, 1), tolerance = 1e-8)
    expect_error(.nmi(fit[["compositions"]] * 2, fit[["compositions"]]), "row-stochastic")
})

test_that("coordinate consistency returns NA for decreasing K and uses column variance denominator", {
    fit3 <- manual_fit()
    X <- fit3[["data"]]
    A4 <- rbind(fit3[["coordinates"]], c(0.5, 0.5))
    B4 <- diag(4)
    S4 <- diag(4)
    loss <- data.frame(loss = 0, r2 = 1, k_S = 1, k_A = 1)
    fit4 <- archetypes(
        coordinates = A4,
        coefficients = B4,
        compositions = S4,
        loss = loss,
        data = X
    )
    ens <- new_archetypes_ensemble(
        data = X,
        fits = list(k3 = fit3, k4 = fit4),
        grid = data.frame(model_id = c("k3", "k4"), K = c(3L, 4L), replicate = 1L),
        metrics = data.frame(model_id = c("k3", "k4"), K = c(3L, 4L), replicate = 1L, AIC = c(1, 2)),
        prefer_metric = "AIC",
        prefer_direction = "minimize",
        eval_fun = AIC,
        eval_name = "AIC"
    )

    score <- consistency(ens, "coordinates")
    expected_d2 <- mean(.aa_greedy_coordinate_d2(fit3[["coordinates"]], fit4[["coordinates"]]))
    expected <- 1 - expected_d2 / mean(matrixStats::colVars(X))

    expect_equal(consistency(fit3, fit4, what = "coordinates"), expected)
    expect_true(is.na(consistency(fit4, fit3, what = "coordinates")))
    expect_equal(score[["k3", "k4"]], expected)
    expect_true(is.na(score[["k4", "k3"]]))
})
