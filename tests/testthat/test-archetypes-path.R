test_that("archetypes_path fits a shared-data sequence over K", {
    set.seed(42)
    X <- toy_matrix()

    path <- archetypes_path(X, K = 3L, max_iter = 3L, tol_r2 = 0.95)

    expect_s3_class(path, "archetypes_path")
    expect_equal(length(path), 3L)
    expect_equal(names(path), paste0("K", 1:3))
    expect_equal(path$K, 1:3)
    expect_identical(path$data, X)
    expect_equal(as.character(path$call[[1L]]), "archetypes_path")
    expect_true(all(vapply(path$models, function(fit) is.null(fit[["data"]]), logical(1L))))

    fit1 <- path[[1L]]
    fit2 <- path[["K2"]]
    expect_archetypes_fit(fit1, K = 1L, n = nrow(X), p = ncol(X))
    expect_archetypes_fit(fit2, K = 2L, n = nrow(X), p = ncol(X))
    expect_identical(fit1[["data"]], X)
    expect_identical(fit2[["data"]], X)
    expect_equal(as.integer(fit2[["call"]][["K"]]), 2L)
    expect_error(path[["K"]], "Unknown archetypes path model")
})

test_that("single bracket subsets archetypes_path models", {
    set.seed(42)
    X <- toy_matrix()
    path <- archetypes_path(X, K = 2:4, max_iter = 3L, tol_r2 = 0.95)

    first <- path[1L]
    named <- path[c("K2", "K4")]

    expect_s3_class(first, "archetypes_path")
    expect_equal(first$K, 2L)
    expect_equal(names(first), "K2")
    expect_archetypes_fit(first[[1L]], K = 2L, n = nrow(X), p = ncol(X))
    expect_identical(first$data, X)

    expect_s3_class(named, "archetypes_path")
    expect_equal(named$K, c(2L, 4L))
    expect_equal(names(named), c("K2", "K4"))
    expect_archetypes_fit(named[["K4"]], K = 4L, n = nrow(X), p = ncol(X))
    expect_error(path["K"], "Unknown archetypes path model")
})

test_that("run_aa returns a path only for vector K", {
    set.seed(42)
    X <- toy_matrix()

    scalar <- run_aa(X, K = 2L, max_iter = 3L, tol_r2 = 0.95)
    path <- run_aa(X, K = c(2L, 1L, 2L), max_iter = 3L, tol_r2 = 0.95)

    expect_s3_class(scalar, "archetypes")
    expect_s3_class(path, "archetypes_path")
    expect_equal(path$K, 1:2)
    expect_equal(as.character(path$call[[1L]]), "run_aa")
    expect_equal(as.character(path[[1L]][["call"]][[1L]]), "run_aa")
})

test_that("archetypes_path supports formula input", {
    set.seed(42)
    path <- archetypes_path(Species ~ ., data = iris, K = 2L, max_iter = 3L, tol_r2 = 0.95)

    expect_s3_class(path, "archetypes_path")
    expect_equal(path$K, 1:2)
    expect_s3_class(path[["K2"]], "archetypes")
    expect_equal(path$formula, Species ~ .)
})

test_that("archetypes_path supports fda fd input", {
    testthat::skip_if_not_installed("fda")

    basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = 5L)
    coefs <- matrix(
        seq_len(25L) / 25,
        nrow = 5L,
        ncol = 5L,
        dimnames = list(paste0("b", 1:5), paste0("x", 1:5))
    )
    fd <- fda::fd(coefs, basis)

    path <- suppressWarnings(archetypes_path(fd, K = 2L, max_iter = 1L, sd_threshold = 0))
    fit <- path[["K2"]]

    expect_s3_class(path, "archetypes_path")
    expect_s3_class(path$data, "fd")
    expect_s3_class(fit, "archetypes")
    expect_s3_class(fit[["data"]], "fd")
    expect_s3_class(coordinates(fit), "fd")
    expect_matrix_dim(fit[["coefficients"]], 2L, 5L)
    expect_matrix_dim(fit[["compositions"]], 5L, 2L)
})


test_that("archetypes_path supports Frank-Wolfe method", {
    set.seed(42)
    X <- toy_matrix()

    path <- suppressWarnings(archetypes_path(
        X,
        K = 2:3,
        method = "fw",
        max_iter = 2L,
        tol_r2 = 0.95
    ))

    expect_s3_class(path, "archetypes_path")
    expect_equal(path$method, "fw")
    expect_archetypes_fit(path[["K2"]], K = 2L, n = nrow(X), p = ncol(X))
    expect_archetypes_fit(path[["K3"]], K = 3L, n = nrow(X), p = ncol(X))
})

test_that("archetypes_path validates K and ambiguous matrix initialization", {
    X <- toy_matrix()

    expect_error(archetypes_path(X, K = c(1L, NA_integer_)), "`K`")
    expect_error(archetypes_path(X, K = 0L), "`K`")
    expect_error(archetypes_path(X, K = nrow(X) + 1L), "number of samples")
    expect_error(
        archetypes_path(X, K = 1:2, init = X[1:2, , drop = FALSE]),
        "`init`"
    )
    expect_no_error(suppressWarnings(archetypes_path(
        X,
        K = 2:3,
        method = "nnls",
        max_kappa = Inf,
        max_iter = 1L
    )))
    expect_error(
        archetypes_path(X, K = 2:3, method = "pgd", max_kappa = Inf, max_iter = 1L),
        "unused"
    )
})

test_that("screeplot.archetypes_path scores final metrics", {
    set.seed(42)
    X <- toy_matrix()
    path <- archetypes_path(X, K = 1:2, max_iter = 3L, tol_r2 = 0.95)

    default <- screeplot(path, plot = FALSE)
    expect_s3_class(default, "data.frame")
    expect_named(default, c("K", "metric", "value", "loss", "r2", "n_iter", "converged"))
    expect_equal(default[["metric"]], rep("AIC", 2L))
    expect_equal(default[["K"]], 1:2)
    expect_true(is.na(default[["value"]][[1L]]))

    loss <- screeplot(path, y = "loss", plot = FALSE)
    expect_equal(loss[["metric"]], rep("loss", 2L))
    expect_equal(loss[["value"]], vapply(1:2, function(i) tail(path[[i]][["loss"]][["loss"]], 1L), numeric(1L)))

    fn <- screeplot(path, y = function(fit) ncol(compositions(fit)), plot = FALSE)
    expect_equal(fn[["metric"]], rep("function", 2L))
    expect_equal(fn[["value"]], 1:2)

    expect_error(screeplot(path, y = "NULL", plot = FALSE), "`y`")
    expect_error(screeplot(path, y = function(fit) c(1, 2), plot = FALSE), "single numeric")
})
