test_that("kernel PGD matches Euclidean PGD for the linear kernel", {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1,
            0.4, 0.3,
            0.7, 0.2
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(paste0("x", 1:5), c("u", "v"))
    )
    ind <- c(1L, 2L, 3L)
    B0 <- onehot(ind, sparse = FALSE, nc = nrow(X))
    rownames(B0) <- paste0("A", seq_along(ind))
    colnames(B0) <- rownames(X)

    A0 <- X[ind, , drop = FALSE]
    rownames(A0) <- rownames(B0)
    pgd <- suppressWarnings(archetypes_pgd(
        X,
        K = 3L,
        init = A0,
        scale = FALSE,
        max_iter = 3L,
        tol_r2 = 1,
        eps = 1e-8
    ))
    kernel <- suppressWarnings(archetypes_kernel_pgd(
        data = X,
        K = 3L,
        gram = tcrossprod(X),
        init = B0,
        max_iter = 3L,
        tol_r2 = 1,
        eps = 1e-8
    ))

    expect_s3_class(kernel, "kernel_archetypes")
    expect_equal(kernel[["loss"]][["rss"]], pgd[["loss"]][["rss"]], tolerance = 1e-6)
    expect_equal(kernel[["loss"]][["r2"]], pgd[["loss"]][["r2"]], tolerance = 1e-7)
    expect_equal(unname(kernel[["coefficients"]]), unname(pgd[["coefficients"]]), tolerance = 1e-7)
    expect_equal(unname(kernel[["compositions"]]), unname(pgd[["compositions"]]), tolerance = 1e-7)
})

test_that("kernel PGD fits an RBF kernel and returns proxy coordinates", {
    set.seed(1)
    X <- toy_matrix()
    X <- X[seq_len(30L), , drop = FALSE]

    fit <- suppressWarnings(archetypes_kernel_pgd(
        data = X,
        K = 3L,
        kernel = "rbf",
        kernel_args = list(gamma = 5),
        max_iter = 5L,
        tol_r2 = 0.95
    ))

    expect_s3_class(fit, "kernel_archetypes")
    expect_matrix_dim(fit[["coefficients"]], 3L, nrow(X))
    expect_matrix_dim(fit[["compositions"]], nrow(X), 3L)
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_true(all(is.finite(fit[["loss"]][["rss"]])))
    expect_true(all(diff(fit[["loss"]][["rss"]]) <= 1e-8))
    expect_equal(fit[["coordinates_proxy"]], fit[["coefficients"]] %*% X)
})

test_that("kernel PGD passes refinement_steps to furthest_sum initialization", {
    set.seed(1)
    X <- toy_matrix()[1:20, , drop = FALSE]

    fit <- suppressWarnings(archetypes_kernel_pgd(
        data = X,
        K = 4L,
        kernel = "linear",
        init = "furthest_sum",
        init_args = list(refinement_steps = 100L),
        max_iter = 1L,
        tol_r2 = 1
    ))

    expect_s3_class(fit, "kernel_archetypes")
    expect_matrix_dim(fit[["init"]], 4L, nrow(X))
    expect_row_stochastic(fit[["init"]])
})

test_that("run_aa dispatches to kernel PGD", {
    X <- toy_matrix()[1:20, , drop = FALSE]

    fit <- suppressWarnings(run_aa(
        X,
        K = 3L,
        method = "kernel",
        kernel = "rbf",
        kernel_args = list(sigma = 0.5),
        max_iter = 2L
    ))

    expect_s3_class(fit, "kernel_archetypes")
    expect_identical(as.character(fit[["call"]][[1L]]), "run_aa")
    expect_matrix_dim(fit[["coefficients"]], 3L, nrow(X))
    expect_equal(fit[["coordinates_proxy"]], fit[["coefficients"]] %*% X)
})

test_that("kernel PGD validates Gram and kernel inputs", {
    X <- toy_matrix()[1:8, , drop = FALSE]
    G <- tcrossprod(X)

    expect_error(archetypes_kernel_pgd(data = X, K = 2L), "Supply either")
    expect_error(
        archetypes_kernel_pgd(data = X, K = 2L, gram = G, kernel = "linear"),
        "exactly one"
    )
    expect_error(archetypes_kernel_pgd(data = X[-1L, ], K = 2L, gram = G), "rows")
    expect_error(archetypes_kernel_pgd(K = 2L, gram = G[1:7, ]), "square")

    G_bad <- G
    G_bad[1, 2] <- G_bad[1, 2] + 1
    expect_error(archetypes_kernel_pgd(K = 2L, gram = G_bad), "symmetric")

    G_bad <- G
    G_bad[1, 1] <- NA_real_
    expect_error(archetypes_kernel_pgd(K = 2L, gram = G_bad), "missing or non-finite")

    custom <- function(X) tcrossprod(X)
    fit <- suppressWarnings(archetypes_kernel_pgd(
        data = X,
        K = 2L,
        kernel = custom,
        init = c(1L, 2L),
        max_iter = 1L
    ))
    expect_s3_class(fit, "kernel_archetypes")
})

test_that("kernel archetypes methods expose names, residuals, and proxy plots", {
    X <- toy_matrix()[1:12, , drop = FALSE]
    fit <- suppressWarnings(archetypes_kernel_pgd(
        data = X,
        K = 3L,
        kernel = "linear",
        init = c(1L, 2L, 3L),
        max_iter = 1L
    ))

    expect_equal(names(fit), c("A1", "A2", "A3"))
    names(fit) <- c("left", "right", "top")
    expect_equal(names(fit), c("left", "right", "top"))
    expect_equal(rownames(coefficients(fit)), names(fit))

    res <- residuals(fit)
    expect_length(res, nrow(X))
    expect_true(all(is.finite(res)))
    expect_error(fitted(fit), "not defined")

    grDevices::pdf(tempfile(fileext = ".pdf"))
    expect_silent(plot(fit, what = "loss"))
    expect_silent(plot(fit, what = "coordinates"))
    grDevices::dev.off()
})
