toy_path <- function() {
    path <- system.file("extdata", "toy.csv", package = "YAAAP")
    if (!nzchar(path))
        path <- file.path("..", "..", "inst", "extdata", "toy.csv")
    path
}

toy_data <- function() {
    read.csv(toy_path())
}

toy_matrix <- function() {
    as.matrix(toy_data())
}

expect_row_stochastic <- function(x, tolerance = 1e-8) {
    expect_true(all(x >= -tolerance))
    expect_equal(unname(rowSums(x)), rep(1, nrow(x)), tolerance = tolerance)
}

expect_matrix_dim <- function(x, nrow, ncol) {
    expect_true(is.matrix(x))
    expect_equal(dim(x), c(nrow, ncol))
}

manual_fit <- function() {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1,
            0.4, 0.3
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(NULL, c("x", "y"))
    )
    A <- X[1:3, , drop = FALSE]
    B <- rbind(
        c(1, 0, 0, 0),
        c(0, 1, 0, 0),
        c(0, 0, 1, 0)
    )
    S <- matrix(
        c(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            0.3, 0.4, 0.3
        ),
        ncol = 3,
        byrow = TRUE
    )
    loss <- data.frame(rss = c(2, 1, 0.5), r2 = c(0, 0.5, 0.75), k_S = 1, k_A = 1)

    archetypes(
        coordinates = A,
        coefficients = B,
        compositions = S,
        loss = loss,
        data = X
    )
}

directional_matrix <- function(n = 90L) {
    theta <- seq(0.05, pi / 2 - 0.05, length.out = n)
    X <- cbind(cos(theta), sin(theta), 0.35 * sin(2 * theta))
    X <- X / sqrt(rowSums(X * X))
    scale <- 0.5 + seq_len(n) / n
    X * scale
}
