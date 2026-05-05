#' Perform Archetypal Analysis using NNLS
#'
#' @param data data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param init function, method string, or numeric coordinate matrix to initialize
#'   archetypes (default: `"furthest_sum"`). When a matrix is supplied it must
#'   have dimension `K x ncol(data)`. Rows outside the convex hull of `data` are
#'   projected into it with a warning; row names, when present, are used as
#'   archetype names.
#' @param init_args list of additional arguments for the initialization function
#' @param weights optional vector of sample weights (default: NULL)
#' @param robust whether to use Tukey bisquare row reweighting (default: FALSE)
#' @param tukey_c tuning constant for Tukey bisquare weights (default: 4.685)
#' @param sd_threshold threshold for feature standard deviation to filter
#'   low-variance features (default: 1e-6)
#' @param max_iter maximum number of iterations (default: 100)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-6)
#' @param tol_r2 convergence tolerance based on R^2 (default: 0.9999)
#' @param max_kappa maximum condition number for archetypes (default: 1000)
#' @param eps small positive number to ensure numerical stability
#'   (default: 0 for sparse input 1e-8 for dense)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param ols_solver method for solving the OLS problem min_A X = SA (default: "qr")
#' @param bigM large constant to enforce simplex constraint (default: 200)
#'
#' @returns An object of class \code{\link{archetypes}}
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "YAAAP"))
#' archetypes_nnls(as.matrix(toy), K = 3)
#'
#' @references Alcacer, A., Epifanio, I., Mair, S., & Mørup, M. (2025).
#' A Survey on Archetypal Analysis. *arXiv preprint arXiv:2504.12392*.
#' \url{https://arxiv.org/abs/2504.12392}
#'
#' @export
archetypes_nnls <- function(data,
                            K,
                            init = "furthest_sum",
                            init_args = list(),
                            weights = NULL,
                            robust = FALSE,
                            tukey_c = 4.685,
                            sd_threshold = 1e-6,
                            max_iter = 100L,
                            tol = 1e-6,
                            tol_r2 = 0.9999,
                            max_kappa = 1000,
                            eps = ifelse(inherits(data, "sparseMatrix"), 0, 1e-8),
                            verbose = FALSE,
                            # NNLS specific
                            ols_solver = c("qr", "ginv", "BFGS"),
                            bigM = 200) {
    .aa_run(
        call = match.call(),
        data = data,
        K = K,
        method = "nnls",
        init = init,
        init_args = init_args,
        weights = weights,
        robust = robust,
        tukey_c = tukey_c,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        method_args = list(
            ols_solver = ols_solver,
            bigM = bigM
        )
    )
}

.aa_fit_nnls <- function(X,
                         weight_fun,
                         max_iter,
                         tol,
                         tol_r2,
                         max_kappa,
                         eps,
                         verbose,
                         A,
                         B,
                         S,
                         loss,
                         ols_solver) {
    # Nomenclature following arXiv:2504.12392v1:
    #   X ~ SA (N x M) Data Matrix
    #   A = BX (K x M) Archetypes
    #   B = (K x N) Archetypes Coefficients (base transform, C in the paper)
    #   S = (N x K) Archetypes Scores (new coordinates)
    A0 <- A
    nnls_svd_kappa_threshold <- 500
    loss_terms <- .aa_loss_terms(
        X,
        A,
        S,
        weight_fun = weight_fun,
        return_S_terms = max_kappa > 1
    )
    row_weights <- loss_terms[["row_weights"]]
    x2 <- loss_terms[["x2"]]  # cache to avoid recomputing in every iteration
    loss[["k_A"]][1L] <- kappa(A, exact = TRUE)
    loss <- .aa_update_loss(
        loss,
        1L,
        loss_terms,
        verbose = verbose,
        max_kappa = max_kappa
    )
    use_svd_for_S <- loss[["k_A"]][1L] > nnls_svd_kappa_threshold
    converged <- FALSE


    # Main Loop  --------------------------------------------------------------

    if (verbose) message("Starting main loop...")
    for (i in seq_len(max_iter)) {
        check_kappa <- i %% 10L == 0L  # Check kappa every 10 iterations
        # Step
        S <- fit_nnls(X, t(A), eps = eps, use_svd = use_svd_for_S) # Project X to A-simplex
        A <- fit_ols(S, X, method = ols_solver, row_weights = row_weights)
        B <- fit_nnls(A, t(X), eps = eps, use_svd = FALSE) # Project A to X-simplex
        A <- B %*% X
        loss_terms <- .aa_loss_terms(
            X,
            A,
            S,
            weight_fun = weight_fun,
            return_S_terms = check_kappa && max_kappa > 1,
            x2 = x2
        )
        row_weights <- loss_terms[["row_weights"]]

        loss <- .aa_update_loss(
            loss,
            i + 1L,
            loss_terms,
            verbose = verbose,
            max_kappa = max_kappa
        )
        if (check_kappa) {
            k_A <- loss[["k_A"]][i + 1L]
            # TODO: exact kappa already computes the SVD? maybe we should resuse it.
            if (is.na(k_A))
                loss[["k_A"]][i + 1L] <- k_A <- kappa(A, exact = TRUE)
            use_svd_for_S <- k_A > nnls_svd_kappa_threshold
        }

        # Check convergence
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }

    list(
        A0 = A0,
        A = A,
        B = B,
        S = S,
        delta = 0,
        i = i,
        loss = loss,
        converged = converged
    )
}

# Fit Non-negative Least Squares for every row of Y
#
# This function solves the problem: $\min_{B} ||Y - BX||_2 s.t. B >= 0$
#
# @param Y data matrix (rows = samples, columns = dimensions)
# @param X data matrix (rows = samples, columns = dimensions)
# @param eps small positive number to ensure numerical stability (default: 1e-8)
# @param project function to project the results onto the simplex (default: `proj_l1`)
# @param use_svd logical, whether to use SVD for dimensionality reduction (default: FALSE)
fit_nnls <- function(Y, X, eps = 1e-8, project = proj_l1, use_svd = FALSE) {
    # min ||Y - Beta %*% X||_2 s.t. Beta >= eps

    if (use_svd) {
        s <- svd(X)
        # TODO: choose rank based on explained variance
        X <- with(s, t(v * d))
        Y <- Y %*% s$u
    }
    # TODO: parallelize
    Beta <- matrix(eps, nrow = nrow(Y), ncol = ncol(X))
    for (i in seq_len(nrow(Y)))
        Beta[i, ] <- stats::coef(nnls::nnls(X, Y[i, ]))
    Beta <- project(Beta, eps = eps)
    Beta
}

# Fit Ordinary Least Squares (OLS) for every column of X
#
# This function solves the problem: $\min_{A} ||X - SA||_F$
#
# @param S data matrix (rows = samples, columns = archetypes)
# @param X data matrix (rows = samples, columns = dimensions)
# @param method method to use for solving the OLS problem (default: "qr")
# @param a0 initial guess for the coefficients (optional)
# @param row_weights optional vector of row weights
# @param ... additional arguments passed to the solver
fit_ols <- function(S, X, method, a0 = NULL, row_weights = NULL, ...) {
    # Solve min ||X - S %*% A||_F
    if (!is.null(row_weights)) {
        sqrt_weights <- sqrt(row_weights)
        S <- S * sqrt_weights
        X <- X * sqrt_weights
    }

    if (tolower(method) == "qr") return(qr.solve(S, X))
    if (tolower(method) == "ginv") return(MASS::ginv(S) %*% X)
    # TODO: test if computing Gram matrix is faster than using `qr.solve` or `ginv`
    # TODO: if "qr" return StS from Gram matrix as it's part of computation
    # A = solve(t(S) %*% S) %*% t(S) %*% X

    # method is one of the `optim` methods (if not an error will be thrown)
    M <- ncol(X)
    K <- ncol(S)
    if (is.null(a0))  # random initial guess
        a0 <- apply(X, 2L, function(x) stats::rnorm(K, mean(x), stats::sd(x)))

    a0 <- as.vector(a0)
    stopifnot(length(a0) == K * M)

    # Squared Frobenius objective for smooth BFGS optimization.
    fn <- function(a) {
        A <- matrix(a, nrow = K, ncol = M)
        R <- X - S %*% A
        sum(R * R)
    }

    # d/dA ||X - S %*% A||_F^2 = -2 * t(S) %*% (X - S %*% A)
    gr <- function(a) {
        A <- matrix(a, nrow = K, ncol = M)
        R <- X - S %*% A
        as.vector(-2 * crossprod(S, R))
    }

    res <- stats::optim(a0, fn, gr, method = method)
    matrix(res$par, nrow = K, ncol = M)
}
