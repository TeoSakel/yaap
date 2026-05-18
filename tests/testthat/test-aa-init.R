test_that("aa_init covers every initialization method", {
    X <- scale(toy_matrix())
    methods <- c(
        "random",
        "furthest_first",
        "kmeans_pp",
        "furthest_sum",
        "aa_pp",
        "hull_outmost"
    )

    for (method in methods) {
        set.seed(1)
        args <- list(X = X, K = 3L, method = method)
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

    expect_error(aa_init(X, K = 3L, method = "coreset_initfn"), "should be one of")
    expect_error(aa_init(X, K = 3L, method = "aa_pp_mc", batch_size = 25L), "should be one of")
    expect_error(aa_init(X, K = 3L, method = "random", batch_size = 2L), "batch_size")
    expect_error(aa_init(X, K = 3L, method = "random", batch_size = nrow(X) + 1L), "batch_size")
    expect_error(aa_init(X, K = 3L, method = "random", batch_type = "bad"), "should be one of")
    expect_error(aa_init(X, K = 3L, method = "random", batch_replace = NA), "batch_replace")
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

test_that("furthest_sum refinement returns valid unique indices", {
    X <- scale(toy_matrix())[1:12, , drop = FALSE]

    set.seed(1)
    ind <- furthest_sum(X, K = 4L)

    expect_length(ind, 4L)
    expect_equal(length(unique(ind)), 4L)
    expect_true(all(ind >= 1L & ind <= nrow(X)))
})

test_that("furthest_sum accepts no refinement and caps excess refinement", {
    X <- scale(toy_matrix())[1:5, , drop = FALSE]

    set.seed(2)
    no_refinement <- furthest_sum(X, K = 4L, refinement_steps = 0L)
    set.seed(2)
    capped_refinement <- furthest_sum(X, K = 4L, refinement_steps = 100L)

    for (ind in list(no_refinement, capped_refinement)) {
        expect_length(ind, 4L)
        expect_equal(length(unique(ind)), 4L)
        expect_true(all(ind >= 1L & ind <= nrow(X)))
    }
})

test_that("furthest_sum handles K equal to the number of samples", {
    X <- scale(toy_matrix())[1:5, , drop = FALSE]

    set.seed(2)
    ind <- furthest_sum(X, K = nrow(X), refinement_steps = 100L)

    expect_length(ind, nrow(X))
    expect_equal(sort(ind), seq_len(nrow(X)))
})

test_that("furthest_sum validates refinement_steps", {
    X <- scale(toy_matrix())[1:8, , drop = FALSE]

    expect_error(furthest_sum(X, K = 3L, refinement_steps = -1L), "refinement_steps")
    expect_error(furthest_sum(X, K = 3L, refinement_steps = 1.5), "refinement_steps")
    expect_error(furthest_sum(X, K = 3L, refinement_steps = NA_integer_), "refinement_steps")
    expect_error(furthest_sum(X, K = 3L, refinement_steps = c(1L, 2L)), "refinement_steps")
})

test_that("aa_init passes refinement_steps to furthest_sum", {
    X <- scale(toy_matrix())[1:20, , drop = FALSE]

    set.seed(3)
    init <- aa_init(X, K = 4L, method = "furthest_sum", refinement_steps = 100L)

    expect_named(init, c("A", "B"))
    expect_matrix_dim(init[["A"]], 4L, 2L)
    expect_matrix_dim(init[["B"]], 4L, 20L)
    expect_row_stochastic(init[["B"]])
    expect_equal(init[["A"]], init[["B"]] %*% X, tolerance = 1e-8)
})

test_that("aa_init applies distal batching by default", {
    X <- rbind(
        matrix(0, nrow = 8L, ncol = 2L),
        c(10, 0),
        c(-5, 8.660254),
        c(-5, -8.660254)
    )

    set.seed(1)
    init <- aa_init(X, K = 3L, method = "random", batch_size = 3L)

    selected <- which(colSums(init[["B"]]) > 0)
    expect_equal(sort(selected), 9:11)
    expect_equal(init[["A"]], init[["B"]] %*% X, tolerance = 1e-8)
})

test_that("batched dirichlet uses only sampled candidate rows", {
    X <- scale(toy_matrix())[1:20, , drop = FALSE]

    set.seed(11)
    init <- aa_init(X, K = 4L, method = "dirichlet", batch_size = 6L)

    active <- which(colSums(init[["B"]]) > 0)
    expect_length(active, 6L)
    expect_row_stochastic(init[["B"]])
    expect_equal(init[["A"]], init[["B"]] %*% X, tolerance = 1e-8)
})

test_that("aa_pp accepts batch_size as its mini-batch approximation", {
    X <- scale(toy_matrix())[1:30, , drop = FALSE]

    set.seed(12)
    init <- aa_init(X, K = 4L, method = "aa_pp", batch_size = 10L, batch_type = "uniform")

    expect_named(init, c("A", "B"))
    expect_matrix_dim(init[["A"]], 4L, 2L)
    expect_matrix_dim(init[["B"]], 4L, 30L)
    expect_row_stochastic(init[["B"]])
    expect_equal(init[["A"]], init[["B"]] %*% X, tolerance = 1e-8)
})

test_that("aa_init uses stable archetype names instead of selected data row names", {
    X <- scale(toy_matrix())
    rownames(X) <- paste0("sample_", seq_len(nrow(X)))

    init <- aa_init(X, K = 3L, method = "random")

    expect_equal(rownames(init[["A"]]), c("A1", "A2", "A3"))
    expect_equal(rownames(init[["B"]]), c("A1", "A2", "A3"))

    named <- .aa_ind_to_init(X, stats::setNames(c(2L, 4L), c("left", "right")), sparse = FALSE)

    expect_equal(rownames(named[["A"]]), c("left", "right"))
    expect_equal(rownames(named[["B"]]), c("left", "right"))
})
