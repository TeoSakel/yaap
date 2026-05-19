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
        x = tcrossprod(X),
        K = 3L,
        kernel = "precomputed",
        data = X,
        init = B0,
        max_iter = 3L,
        tol_r2 = 1,
        eps = 1e-8
    ))

    expect_s3_class(kernel, "kernel_archetypes")
    expect_equal(kernel[["loss"]][["loss"]], pgd[["loss"]][["loss"]], tolerance = 1e-6)
    expect_equal(kernel[["loss"]][["r2"]], pgd[["loss"]][["r2"]], tolerance = 1e-7)
    row_order <- function(A) do.call(order, as.data.frame(round(A, 7L)))
    kernel_order <- row_order(coordinates(kernel))
    pgd_order <- row_order(coordinates(pgd))
    expect_equal(unname(coordinates(kernel)[kernel_order, , drop = FALSE]),
                 unname(coordinates(pgd)[pgd_order, , drop = FALSE]),
                 tolerance = 1e-7)
    expect_equal(unname(kernel[["coefficients"]][kernel_order, , drop = FALSE]),
                 unname(pgd[["coefficients"]][pgd_order, , drop = FALSE]),
                 tolerance = 1e-7)
})

test_that("kernel PGD fits an RBF kernel and returns coordinates", {
    set.seed(1)
    X <- toy_matrix()
    X <- X[seq_len(30L), , drop = FALSE]

    fit <- suppressWarnings(archetypes_kernel_pgd(
        x = X,
        K = 3L,
        kernel = "rbf",
        kernel_args = list(sigma = 0.5),
        max_iter = 5L,
        tol_r2 = 0.95
    ))

    expect_s3_class(fit, "kernel_archetypes")
    expect_matrix_dim(fit[["coefficients"]], 3L, nrow(X))
    expect_matrix_dim(fit[["compositions"]], nrow(X), 3L)
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_true(is_all_finite(fit[["loss"]][["loss"]]))
    expect_named(fit[["loss"]], c("loss", "r2"))
    expect_true(all(diff(fit[["loss"]][["loss"]]) <= 1e-8))
    expect_equal(coordinates(fit), fit[["coefficients"]] %*% X)
    expect_equal(coordinates(fit), coordinates(fit))
    expect_error(
        archetypes_kernel_pgd(
            x = X,
            K = 3L,
            kernel = "rbf",
            kernel_args = list(gamma = 5)
        ),
        "RBF kernels use `sigma`"
    )
    expect_error(
        archetypes_kernel_pgd(
            x = X,
            K = 3L,
            kernel = "laplace",
            kernel_args = list(gamma = 5)
        ),
        "Laplace kernels use `sigma`"
    )
})
test_that("coordinates.kernel_archetypes computes the input-space proxy on demand", {
    X <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 3L, byrow = TRUE)
    rownames(X) <- c("s1", "s2", "s3")
    colnames(X) <- c("x", "y")
    B <- matrix(c(0.5, 0.5, 0, 0, 0.25, 0.75), nrow = 2L, byrow = TRUE)
    rownames(B) <- c("A1", "A2")
    colnames(B) <- rownames(X)
    S <- matrix(c(1, 0, 0, 1, 0.5, 0.5), nrow = 3L, byrow = TRUE)
    rownames(S) <- rownames(X)
    colnames(S) <- rownames(B)
    fit <- kernel_archetypes(
        coefficients = B,
        compositions = S,
        gram = diag(3L),
        data = X
    )

    expect_false("coordinates" %in% names(fit))
    expect_equal(coordinates(fit), fit[["coefficients"]] %*% X)
})

test_that("coordinates.kernel_archetypes preserves Matrix dispatch for sparse data", {
    X <- Matrix::Matrix(
        matrix(c(0, 0, 1, 0, 0, 1), nrow = 3L, byrow = TRUE),
        sparse = TRUE
    )
    rownames(X) <- c("s1", "s2", "s3")
    colnames(X) <- c("x", "y")
    B <- matrix(c(0.5, 0.5, 0, 0, 0.25, 0.75), nrow = 2L, byrow = TRUE)
    rownames(B) <- c("A1", "A2")
    colnames(B) <- rownames(X)
    S <- matrix(c(1, 0, 0, 1, 0.5, 0.5), nrow = 3L, byrow = TRUE)
    rownames(S) <- rownames(X)
    colnames(S) <- rownames(B)

    fit <- kernel_archetypes(
        coefficients = B,
        compositions = S,
        gram = diag(3L),
        data = X
    )

    expect_true(inherits(coordinates(fit), "Matrix"))
})


test_that("kernel PGD passes refinement_steps to furthest_sum initialization", {
    set.seed(1)
    X <- toy_matrix()[1:20, , drop = FALSE]

    fit <- suppressWarnings(archetypes_kernel_pgd(
        x = X,
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

test_that("kernel PGD supports dirichlet initialization", {
    set.seed(4)
    X <- toy_matrix()[1:20, , drop = FALSE]

    fit <- suppressWarnings(archetypes_kernel_pgd(
        x = X,
        K = 3L,
        kernel = "rbf",
        kernel_args = list(sigma = 0.5),
        init = "dirichlet",
        init_args = list(alpha = 0.5, batch_size = 6L),
        max_iter = 1L,
        tol_r2 = 1
    ))

    expect_s3_class(fit, "kernel_archetypes")
    expect_matrix_dim(fit[["init"]], 3L, nrow(X))
    expect_row_stochastic(fit[["init"]])
    expect_length(which(colSums(fit[["init"]]) > 0), 6L)
})

test_that("kernel PGD accepts batching for row-index initializers", {
    X <- rbind(
        matrix(0, nrow = 8L, ncol = 2L),
        c(10, 0),
        c(-5, 8.660254),
        c(-5, -8.660254)
    )

    set.seed(5)
    fit <- suppressWarnings(archetypes_kernel_pgd(
        x = X,
        K = 3L,
        kernel = "linear",
        init = "random",
        init_args = list(batch_size = 3L),
        max_iter = 1L,
        tol_r2 = 1
    ))

    expect_equal(sort(which(colSums(fit[["init"]]) > 0)), 9:11)
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
    expect_equal(coordinates(fit), fit[["coefficients"]] %*% X)
})

test_that("kernel PGD warns and ignores explicit run_aa scaling", {
    X <- toy_matrix()[1:8, , drop = FALSE]

    expect_warning(
        fit <- run_aa(X, K = 2L, method = "kernel", kernel = "linear", scale = TRUE),
        "ignored"
    )
    expect_s3_class(fit, "kernel_archetypes")

    expect_warning(
        fit <- run_aa(X, K = 2L, method = "kernel", kernel = "linear", scale = c(1, 1)),
        "ignored"
    )
    expect_s3_class(fit, "kernel_archetypes")
})

test_that("kernel PGD validates Gram and kernel inputs", {
    X <- toy_matrix()[1:8, , drop = FALSE]
    G <- tcrossprod(X)

    expect_error(
        archetypes_kernel_pgd(x = X, K = 2L, kernel = "linear", data = X),
        "only used with `kernel = 'precomputed'`"
    )
    expect_error(
        archetypes_kernel_pgd(x = G, data = X[-1L, ], K = 2L, kernel = "precomputed"),
        "rows"
    )
    expect_error(archetypes_kernel_pgd(x = G[1:7, ], K = 2L, kernel = "precomputed"), "square")

    G_bad <- G
    G_bad[1, 2] <- G_bad[1, 2] + 1
    expect_error(archetypes_kernel_pgd(x = G_bad, K = 2L, kernel = "precomputed"), "symmetric")

    G_bad <- G
    G_bad[1, 1] <- NA_real_
    expect_error(archetypes_kernel_pgd(x = G_bad, K = 2L, kernel = "precomputed"),
                 "missing or non-finite")

    G_bad <- G
    G_bad[1, 1] <- -10
    expect_error(archetypes_kernel_pgd(x = G_bad, K = 2L, kernel = "precomputed"),
                 "positive semidefinite")

    custom <- function(X) tcrossprod(X)
    fit <- suppressWarnings(archetypes_kernel_pgd(
        x = X,
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
        x = X,
        K = 3L,
        kernel = "linear",
        init = c(1L, 2L, 3L),
        max_iter = 1L
    ))

    expect_equal(sort(anames(fit)), c("A1", "A2", "A3"))
    anames(fit) <- c("left", "right", "top")
    expect_equal(anames(fit), c("left", "right", "top"))
    expect_equal(rownames(coefficients(fit)), anames(fit))

    res <- residuals(fit)
    expect_length(res, nrow(X))
    expect_true(is_all_finite(res))
    expect_error(fitted(fit), "not defined")
    expect_error(predict(fit, X), "not currently defined")
    expect_error(predict(fit, X, type = "compositions"), "not currently defined")
    expect_error(predict(fit, X, type = "reconstruction"), "not currently defined")

    grDevices::pdf(tempfile(fileext = ".pdf"))
    expect_silent(plot(fit, what = "loss"))
    expect_silent(plot(fit, what = "compositions", legend = FALSE))
    expect_silent(plot(fit, what = "ternary"))
    expect_silent(plot(fit, what = "coordinates"))
    grDevices::dev.off()
})
