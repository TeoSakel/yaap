test_that("dependency metadata keeps optional solvers optional", {
    desc_path <- system.file("DESCRIPTION", package = "YAAAP")
    if (!nzchar(desc_path))
        desc_path <- testthat::test_path("..", "..", "DESCRIPTION")
    desc <- read.dcf(desc_path)[1L, ]
    imports <- desc[["Imports"]]
    suggests <- desc[["Suggests"]]

    expect_false(grepl("\\barchetypes\\b", suggests))
    expect_false(grepl("\\bquadprog\\b", imports))
    expect_true(grepl("\\bquadprog\\b", suggests))
    expect_false(grepl("\\bMASS\\b", imports))
    expect_false(grepl("\\bMASS\\b", suggests))
})

test_that("internal pseudoinverse handles singular matrices without MASS", {
    ginv <- get(".aa_ginv", asNamespace("YAAAP"))
    S <- matrix(c(1, 2, 2, 4), nrow = 2L)
    G <- ginv(S)

    expect_equal(S %*% G %*% S, S, tolerance = 1e-10)
    expect_equal(G %*% S %*% G, G, tolerance = 1e-10)
})

test_that("row weight validation returns NULL for trivial weights", {
    check_row_weights <- get(".aa_check_row_weights", asNamespace("YAAAP"))

    expect_null(check_row_weights(NULL, 3L))
    expect_null(check_row_weights(c(1, 1, 1), 3L))
    expect_equal(check_row_weights(c(1, 0.5, 1), 3L), c(1, 0.5, 1))
})

test_that("effic uses a pseudo-inverse when cov(X) is singular", {
    effic <- get("effic", asNamespace("YAAAP"))
    ginv <- get(".aa_ginv", asNamespace("YAAAP"))
    X <- matrix(seq_len(12), nrow = 3L)
    Y <- X / 2

    expect_warning(
        eta <- effic(X, Y),
        "pseudo-inverse"
    )
    expect_equal(eta, sum(diag(ginv(stats::cov(X)) %*% stats::cov(Y))))
    expect_true(is_number(eta))
})

test_that("effic keeps Cholesky path for positive definite covariance", {
    effic <- get("effic", asNamespace("YAAAP"))
    X <- as.matrix(iris[1:50, 1:4])
    Y <- X / 2

    expect_warning(eta <- effic(X, Y), NA)
    expect_equal(eta, sum(stats::cov(Y) * chol2inv(chol(stats::cov(X)))))
})
