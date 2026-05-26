truncated_simplex_data <- function() {
    N <- 1000L
    threshold <- 0.8
    noise <- 0
    K <- 3L
    M <- 2L

    set.seed(20250804)
    A <- c(
        cos(0), cos(2 * pi / 3), cos(2 * pi / 3 * 2),
        sin(0), sin(2 * pi / 3), sin(2 * pi / 3 * 2)
    )
    A <- matrix(A, nrow = K, ncol = M)
    S <- -log(matrix(runif(K * N), nrow = N, ncol = K))
    S <- S / rowSums(S)
    S <- S[rowSums(S > threshold) == 0, ]

    S %*% A + noise * rnorm(nrow(S) * M)
}

test_that("pgd delta relaxes coefficient row sums within bounds", {
    X <- truncated_simplex_data()
    deltas <- c(0, 0.25, 0.5)
    rowsums <- vector("list", length(deltas))

    for (i in seq_along(deltas)) {
        set.seed(1)
        delta <- deltas[[i]]
        fit <- archetypes_pgd(
            X,
            K = 3L,
            delta = delta,
            max_iter = 30L,
            tol_r2 = 0.98
        )
        rowsums[[i]] <- rowSums(fit[["coefficients"]])

        expect_true(all(rowsums[[i]] >= 1 - delta - 1e-8))
        expect_true(all(rowsums[[i]] <= 1 + delta + 1e-8))
    }

    expect_equal(unname(rowsums[[1]]), rep(1, 3L), tolerance = 1e-8)
    expect_true(any(abs(rowsums[[2]] - 1) > 1e-4) ||
        any(abs(rowsums[[3]] - 1) > 1e-4))
})

test_that("archetypes constructor accepts coefficient row sums at slack boundary", {
    coefficients <- diag(rep(1 + 0.35 + sqrt(.Machine$double.eps) / 2, 3L))
    compositions <- diag(3L)
    A <- matrix(0, nrow = 3L, ncol = 2L)

    expect_s3_class(
        archetypes(
            A = A,
            coefficients = coefficients,
            compositions = compositions,
            slack = 0.35
        ),
        "archetypes"
    )
})
