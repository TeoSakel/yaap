test_that("effic uses a pseudo-inverse when cov(X) is singular", {
    effic <- get("effic", asNamespace("YAAAP"))
    X <- matrix(seq_len(12), nrow = 3L)
    Y <- X / 2

    expect_warning(
        eta <- effic(X, Y),
        "pseudo-inverse"
    )
    expect_equal(eta, sum(diag(MASS::ginv(stats::cov(X)) %*% stats::cov(Y))))
    expect_true(is_number(eta))
})

test_that("effic keeps Cholesky path for positive definite covariance", {
    effic <- get("effic", asNamespace("YAAAP"))
    X <- as.matrix(iris[1:50, 1:4])
    Y <- X / 2

    expect_warning(eta <- effic(X, Y), NA)
    expect_equal(eta, sum(stats::cov(Y) * chol2inv(chol(stats::cov(X)))))
})
