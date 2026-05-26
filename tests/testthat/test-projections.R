test_that("projection helpers return row-stochastic matrices", {
    mat <- matrix(
        c(
            -1, 2, 0.5,
            0.2, 0.2, 0.2,
            4, 1, -2
        ),
        nrow = 3,
        byrow = TRUE
    )

    simplex <- proj_simplex(mat)
    l1 <- proj_l1(mat)

    expect_matrix_dim(simplex, 3L, 3L)
    expect_matrix_dim(l1, 3L, 3L)
    expect_row_stochastic(simplex)
    expect_row_stochastic(l1)
})

test_that("fit_simplex projects data onto archetype simplex", {
    A <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("A1", "A2", "A3"), c("x", "y"))
    )
    X <- matrix(
        c(
            0.2, 0.2,
            0.8, 0.1
        ),
        ncol = 2,
        byrow = TRUE
    )

    S <- fit_simplex(A, X)

    expect_matrix_dim(S, 2L, 3L)
    expect_named(as.data.frame(S), rownames(A))
    expect_row_stochastic(S)
})

test_that("fit_simplex returns affine simplex coordinates away from the origin", {
    A <- matrix(
        c(
            2, 1,
            5, 1,
            3, 4
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("A1", "A2", "A3"), c("x", "y"))
    )
    S_expected <- matrix(
        c(
            0.2, 0.3, 0.5,
            0.1, 0.8, 0.1
        ),
        ncol = 3L,
        byrow = TRUE,
        dimnames = list(NULL, rownames(A))
    )
    X <- S_expected %*% A

    S <- fit_simplex(A, X)

    expect_equal(S, S_expected, tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(S %*% A, X, tolerance = 1e-8)
})

test_that("fit_simplex nnls doubles bigM when raw row sums drift", {
    A <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2L, byrow = TRUE)
    X <- matrix(c(0.2, 0.2), ncol = 2L)
    seen_bigM <- numeric()
    calls <- 0L

    testthat::local_mocked_bindings(
        .aa_solve_nnls = function(Y, X, use_svd = FALSE) {
            calls <<- calls + 1L
            seen_bigM <<- c(seen_bigM, X[1L, 1L])
            if (calls < 4L) {
                return(matrix(c(0.2, 0.2, 0.1), nrow = 1L))
            }
            matrix(c(0.2, 0.3, 0.5), nrow = 1L)
        },
        .package = "yaap"
    )

    S <- fit_simplex(A, X, bigM = 10)

    expect_equal(unname(seen_bigM), c(10, 20, 40, 80))
    expect_equal(calls, 4L)
    expect_row_stochastic(S)
})

test_that("fit_simplex nnls warns when bigM retries do not reach simplex row sums", {
    A <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2L, byrow = TRUE)
    X <- matrix(c(0.2, 0.2), ncol = 2L)
    calls <- 0L

    testthat::local_mocked_bindings(
        .aa_solve_nnls = function(Y, X, use_svd = FALSE) {
            calls <<- calls + 1L
            matrix(c(0.2, 0.2, 0.1), nrow = 1L)
        },
        .package = "yaap"
    )

    expect_warning(
        S <- fit_simplex(A, X, bigM = 10),
        "did not satisfy the raw simplex constraint"
    )

    expect_equal(calls, 4L)
    expect_row_stochastic(S)
})

test_that("fit_simplex QP path errors when quadprog is unavailable", {
    A <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2L, byrow = TRUE)
    X <- matrix(c(0.2, 0.2), ncol = 2L)

    testthat::local_mocked_bindings(
        .aa_require_namespace_for = function(pkg, feature) {
            stop("mocked missing quadprog", call. = FALSE)
        },
        .package = "yaap"
    )

    expect_error(
        fit_simplex(A, X, method = "QP"),
        "mocked missing quadprog"
    )
})

test_that("matrix initialization reports missing quadprog clearly", {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1
        ),
        ncol = 2L,
        byrow = TRUE
    )
    init <- X

    testthat::local_mocked_bindings(
        .aa_require_namespace_for = function(pkg, feature) {
            err <- simpleError(sprintf("missing %s for %s", pkg, feature))
            err[["pkg"]] <- pkg
            err[["feature"]] <- feature
            class(err) <- c("aa_missing_namespace_error", class(err))
            stop(err)
        },
        .package = "yaap"
    )

    expect_error(
        .aa_matrix_init(X, K = 3L, init = init, eps = 0),
        "Matrix-valued `init` requires the `quadprog` package",
        fixed = TRUE
    )
})
