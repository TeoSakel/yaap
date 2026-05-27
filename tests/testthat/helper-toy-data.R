toy_path <- function() {
    path <- system.file("extdata", "toy.csv", package = "yaap")
    if (!nzchar(path)) {
        path <- file.path("..", "..", "inst", "extdata", "toy.csv")
    }
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

expect_matrix_like_dim <- function(x, nrow, ncol) {
    expect_true(is.matrix(x) || inherits(x, "Matrix"))
    expect_equal(dim(x), c(nrow, ncol))
}

expect_archetypes_fit <- function(fit,
                                  K,
                                  n,
                                  p,
                                  class = "archetypes",
                                  fitted_dim = NULL,
                                  residual_dim = NULL) {
    expect_s3_class(fit, class)
    expect_matrix_dim(coordinates(fit), K, p)
    expect_matrix_like_dim(fit[["coefficients"]], K, n)
    expect_matrix_like_dim(fit[["compositions"]], n, K)
    expect_row_stochastic(fit[["coefficients"]])
    expect_row_stochastic(fit[["compositions"]])
    expect_true(is_all_finite(fit[["loss"]][["loss"]]))

    if (!is.null(fitted_dim)) {
        expect_equal(dim(fitted(fit)), fitted_dim)
    }
    if (!is.null(residual_dim)) {
        expect_equal(dim(residuals(fit)), residual_dim)
    }
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
    rownames(A) <- paste0("A", 1:3)
    B <- rbind(
        c(1, 0, 0, 0),
        c(0, 1, 0, 0),
        c(0, 0, 1, 0)
    )
    rownames(B) <- rownames(A)
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
    colnames(S) <- rownames(A)
    loss <- data.frame(loss = c(2, 1, 0.5), r2 = c(0, 0.5, 0.75), k_S = 1, k_A = 1)

    archetypes(
        A = A,
        coefficients = B,
        compositions = S,
        loss = loss,
        data = X,
        feature_map = .aa_identity_feature_map(A)
    )
}

local_test_pdf <- function(envir = parent.frame()) {
    path <- tempfile(fileext = ".pdf")
    grDevices::pdf(path)
    withr::defer({
        if (names(grDevices::dev.cur()) != "null device") {
            grDevices::dev.off()
        }
        unlink(path)
    }, envir = envir)
    invisible(path)
}

directional_matrix <- function(n = 90L) {
    theta <- seq(0.05, pi / 2 - 0.05, length.out = n)
    X <- cbind(cos(theta), sin(theta), 0.35 * sin(2 * theta))
    X <- X / sqrt(rowSums(X * X))
    scale <- 0.5 + seq_len(n) / n
    X * scale
}

sparse_test_matrix <- function() {
    Matrix::Matrix(
        matrix(
            c(
                1, 0, 2, 0,
                0, 3, 0, 1,
                4, 0, 0, 0,
                0, 5, 1, 0,
                2, 0, 0, 3,
                0, 1, 0, 0
            ),
            nrow = 6L,
            byrow = TRUE,
            dimnames = list(paste0("x", 1:6), paste0("v", 1:4))
        ),
        sparse = TRUE
    )
}

missing_test_matrix <- function() {
    matrix(
        c(
            1, NA, 2, 0,
            0, 3, NA, 1,
            4, 0, 0, NA,
            NA, 5, 1, 0,
            2, NA, 0, 3,
            0, 1, NA, 0
        ),
        nrow = 6L,
        byrow = TRUE,
        dimnames = list(paste0("x", 1:6), paste0("v", 1:4))
    )
}
