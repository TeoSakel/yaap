test_that("aa_init covers every initialization method", {
    X <- scale(toy_matrix())
    methods <- c(
        "uniform_archetypes",
        "furthest_first",
        "kmeans_pp",
        "furthest_sum",
        "coreset_initfn",
        "aa_pp",
        "aa_pp_mc",
        "hull_outmost"
    )

    for (method in methods) {
        set.seed(1)
        args <- list(X = X, K = 3L, method = method)
        if (method == "coreset_initfn")
            args[["m"]] <- 25L
        if (method == "aa_pp_mc")
            args[["batch_size"]] <- 25L
        if (method == "hull_outmost")
            args[["hull_method"]] <- "projected"

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
    expect_error(aa_init(X, K = 3L, method = "hull_outmost", hull_method = "bad"),
                 "should be one of")
    expect_error(aa_init(X, K = 3L, method = "hull_outmost", projected_dim = 0),
                 "projected_dim")
    expect_error(aa_init(X, K = 3L, method = "hull_outmost", n_partitions = 0),
                 "n_partitions")
    expect_error(aa_init(X, K = 3L, method = "hull_outmost", n_projection_max = 0),
                 "n_projection_max")
    expect_error(aa_init(X, K = 3L, method = "hull_outmost", use_unique_candidates = NA),
                 "use_unique_candidates")
})

test_that("aa_init uses stable archetype names instead of selected data row names", {
    X <- scale(toy_matrix())
    rownames(X) <- paste0("sample_", seq_len(nrow(X)))

    init <- aa_init(X, K = 3L, method = "uniform_archetypes")

    expect_equal(rownames(init[["A"]]), c("A1", "A2", "A3"))
    expect_equal(rownames(init[["B"]]), c("A1", "A2", "A3"))

    named <- .ind_to_init(X, stats::setNames(c(2L, 4L), c("left", "right")), sparse = FALSE)

    expect_equal(rownames(named[["A"]]), c("left", "right"))
    expect_equal(rownames(named[["B"]]), c("left", "right"))
})
