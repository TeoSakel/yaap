test_that("aa_init covers every initialization method", {
    X <- scale(toy_matrix())
    methods <- c(
        "uniform_archetypes",
        "furthest_first",
        "kmeans_pp",
        "furthest_sum",
        "coreset_initfn",
        "aa_pp",
        "aa_pp_mc"
    )

    for (method in methods) {
        set.seed(1)
        args <- list(X = X, K = 3L, method = method)
        if (method == "coreset_initfn")
            args[["m"]] <- 25L
        if (method == "aa_pp_mc")
            args[["batch_size"]] <- 25L

        init <- do.call(aa_init, args)

        expect_named(init, c("A", "B"))
        expect_matrix_dim(init[["A"]], 3L, 2L)
        expect_matrix_dim(init[["B"]], 3L, 250L)
        expect_row_stochastic(init[["B"]])
        expect_equal(init[["A"]], init[["B"]] %*% X, tolerance = 1e-8)
    }
})

test_that("aa_init validates required method arguments", {
    X <- scale(toy_matrix())

    expect_error(aa_init(X, K = 3L, method = "coreset_initfn"), "m")
    suppressWarnings(expect_error(aa_init(X, K = 3L, method = "aa_pp_mc")))
})
