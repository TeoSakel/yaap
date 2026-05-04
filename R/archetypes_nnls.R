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

    # Input Checks  -----------------------------------------------------------

    # Generic Checks
    .aa_check_inputs(
        data = data,
        K = K,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps
    )
    # NNLS specific checks
    ols_solver <- match.arg(ols_solver)
    stopifnot("robust must be TRUE or FALSE" =
                  is.logical(robust) && length(robust) == 1L && !is.na(robust))
    stopifnot("tukey_c must be positive" =
                  length(tukey_c) == 1L && is.finite(tukey_c) && tukey_c > 0)
    stopifnot("`bigM` must be greater than 0" = bigM > 0)
    weight_fun <- if (robust) {
        function(row_rss) .aa_bisquare_weights(row_rss, c = tukey_c)
    } else {
        .aa_unit_weights
    }

    # Edge Case checks
    out <- .aa_checks_edge_cases(data, K, verbose)  # edge cases
    if (!is.null(out)) return(out)  # return early if edge case

    # Prepossessing Data  -----------------------------------------------------

    cl <- match.call()
    pre <- .aa_preprocess(data, sd_threshold, weights, verbose, bigM = bigM)
    X <- pre[["X"]]                    # preprocessed data
    # N <- nrow(X)
    undo_scale <- pre[["undo_scale"]]  # function to undo preprocessing
    init <- .aa_preprocess_init(init, data, X)
    rm(pre)


    # Initialize Variables  ---------------------------------------------------

    # Nomenclature following arXiv:2504.12392v1:
    #   X ~ SA (N x M) Data Matrix
    #   A = BX (K x M) Archetypes
    #   B = (K x N) Archetypes Coefficients (base transform, C in the paper)
    #   S = (N x K) Archetypes Scores (new coordinates)

    init_vars <- .aa_init_vars(X, K, init, init_args, eps, max_iter, verbose)
    A <- A0 <- init_vars[["A"]]
    B <- init_vars[["B"]]
    S <- init_vars[["S"]]
    loss <- init_vars[["loss"]]
    rm(init_vars)
    loss_terms <- .aa_loss_terms(
        X,
        A,
        S,
        weight_fun = weight_fun,
        return_S_terms = FALSE
    )
    row_weights <- loss_terms[["row_weights"]]
    loss <- .aa_update_loss(
        loss,
        1L,
        loss_terms,
        verbose = verbose,
        max_kappa = max_kappa
    )
    converged <- FALSE


    # Main Loop  --------------------------------------------------------------

    if (verbose) message("Starting main loop...")
    # TODO: use heuristic to choose between fit_nnls and fit_nnls_svd
    for (i in seq_len(max_iter)) {
        # Step
        S <- fit_nnls(X, t(A), eps = eps) # Project X to A-simplex
        A <- fit_ols(S, X, method = ols_solver, row_weights = row_weights)
        B <- fit_nnls(A, t(X), eps = eps) # Project A to X-simplex
        A <- B %*% X
        loss_terms <- .aa_loss_terms(
            X,
            A,
            S,
            weight_fun = weight_fun,
            return_S_terms = (i + 1L) %% 10L == 0L && max_kappa > 1
        )
        row_weights <- loss_terms[["row_weights"]]

        # Check convergence
        loss <- .aa_update_loss(
            loss,
            i + 1L,
            loss_terms,
            verbose = verbose,
            max_kappa = max_kappa
        )
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }

    # Prepare Output  ---------------------------------------------------------

    .aa_prepare_output(
        call = cl,
        data = data,
        A0 = A0,
        X = X,
        A = A,
        B = B,
        S = S,
        i = i,
        loss = loss,
        converged = converged,
        undo_scale = undo_scale,
        max_iter = max_iter,
        verbose = verbose
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
