test_that("dependency metadata keeps optional solvers optional", {
    desc_path <- system.file("DESCRIPTION", package = "yaap")
    if (!nzchar(desc_path))
        desc_path <- testthat::test_path("..", "..", "DESCRIPTION")
    desc <- read.dcf(desc_path)[1L, ]
    imports <- desc[["Imports"]]
    suggests <- desc[["Suggests"]]

    expect_false(grepl("\\barchetypes\\b", suggests))
    expect_false(grepl("\\bquadprog\\b", imports))
    expect_true(grepl("\\bquadprog\\b", suggests))
    expect_false(grepl("\\bMASS\\b", imports))
    expect_true(grepl("\\bMASS\\b", suggests))
})

test_that("robust weight function resolves MASS psi functions", {
    testthat::skip_if_not_installed("MASS")
    weight_fun <- get(".aa_weight_fun", asNamespace("yaap"))
    row_rss <- c(0, 1, 4, 1e6)

    expect_null(weight_fun(FALSE, list()))
    expect_equal(
        weight_fun(TRUE, list())(row_rss),
        weight_fun("psi.bisquare", list())(row_rss)
    )
    expect_equal(weight_fun(TRUE, list())(rep(0, 4)), rep(1, 4))
    expect_true(weight_fun("psi.huber", list(k = 1.345))(row_rss)[4] <
                    weight_fun("psi.huber", list(k = 1.345))(row_rss)[2])
    expect_true(weight_fun("psi.hampel", list(a = 2, b = 4, c = 8))(row_rss)[4] <
                    weight_fun("psi.hampel", list(a = 2, b = 4, c = 8))(row_rss)[2])
})

test_that("robust weight function accepts custom MASS-contract psi functions", {
    weight_fun <- get(".aa_weight_fun", asNamespace("yaap"))
    psi <- function(u, deriv = 0) {
        if (deriv != 0)
            return(rep(1, length(u)))
        ifelse(abs(u) <= 1, 1, 0.25)
    }

    expect_equal(weight_fun(psi, list())(c(0, 1, 100)), c(1, 1, 0.25))
})

test_that("robust weight function validates inputs and outputs", {
    weight_fun <- get(".aa_weight_fun", asNamespace("yaap"))
    bad_length <- function(u, deriv = 0) 1
    bad_values <- function(u, deriv = 0) rep(NA_real_, length(u))

    expect_error(weight_fun(NA, list()), "robust")
    expect_error(weight_fun("psi.huber", 1), "robust_args")
    expect_error(weight_fun("not_a_psi", list()), "robust")
    expect_error(weight_fun(bad_length, list())(c(0, 1)), "row_weights")
    expect_error(weight_fun(bad_values, list())(c(0, 1)), "row_weights")
})

test_that("row weight validation returns validated weights", {
    check_row_weights <- get(".aa_check_row_weights", asNamespace("yaap"))

    expect_null(check_row_weights(NULL, 3L))
    expect_equal(check_row_weights(c(1, 1, 1), 3L), c(1, 1, 1))
    expect_equal(check_row_weights(c(1, 0.5, 1), 3L), c(1, 0.5, 1))
})

test_that(".aa_effic uses a pseudo-inverse when cov(X) is singular", {
    testthat::skip_if_not_installed("MASS")
    .aa_effic <- get(".aa_effic", asNamespace("yaap"))
    X <- matrix(seq_len(12), nrow = 3L)
    Y <- X / 2

    expect_warning(
        eta <- .aa_effic(X, Y),
        "pseudo-inverse"
    )
    expect_equal(eta, sum(diag(MASS::ginv(stats::cov(X)) %*% stats::cov(Y))))
    expect_true(is_number(eta))
})

test_that(".aa_effic keeps Cholesky path for positive definite covariance", {
    .aa_effic <- get(".aa_effic", asNamespace("yaap"))
    X <- as.matrix(iris[1:50, 1:4])
    Y <- X / 2

    expect_warning(eta <- .aa_effic(X, Y), NA)
    expect_equal(eta, sum(stats::cov(Y) * chol2inv(chol(stats::cov(X)))))
})
