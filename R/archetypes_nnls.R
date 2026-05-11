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
#' @param scale behaves like the `scale.` argument in base R `scale` function.
#'   with the additional option to specify a custom positive-semidefine matrix
#'   for metric embedding (default: TRUE, i.e. column-wise unit variance)
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
#' @param bigM large constant to enforce simplex constraint, or `NULL` to set it automatically.
#' @param max_no_update maximum consecutive iterations without improvement before
#'   considering NNLS stalled (default: 5)
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
                            scale = TRUE,
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
                            bigM = NULL,
                            max_no_update = 5L) {
    .aa_run_aa_default(
        call = match.call(),
        data = data,
        K = K,
        method = "nnls",
        init = init,
        init_args = init_args,
        weights = weights,
        scale = scale,
        robust = robust,
        tukey_c = tukey_c,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        ols_solver = ols_solver,
        bigM = bigM,
        max_no_update = max_no_update
    )
}

.aa_nnls_loss_terms <- function(X, A, S, xss = NULL, ...) {
    # computes RSS using the trace formulation used in the pgd version
    # eventhough this is not the bottleneck in the NNLS version

    iM <- attr(X, "bigM")
    if (!is.null(iM)) {
        X <- X[, -iM, drop = FALSE]
        A <- A[, -iM, drop = FALSE]
    }

    if (is.null(xss)) xss <- norm(X, "F")^2  # computed once
    AAt <- tcrossprod(A)     # (K x K) archetype Gram matrix
    XAt <- tcrossprod(X, A)  # (N x K) projection of data onto archetypes
    StS <- crossprod(S)      # (K x K) archetype score Gram matrix
    rss <- xss - 2 * sum(S * XAt) + sum(StS * AAt)

    list(
        rss = rss,
        xss = xss,
        StS = StS,
        S_weighted = S,
        AAt = AAt,
        XAt = XAt,
        A = A
    )
}

.aa_nnls_weighted_loss_terms <- function(X, A, S, weight_fun, row_xss = NULL, ...) {
    iM <- attr(X, "bigM")
    if (!is.null(iM)) {
        X <- X[, -iM, drop = FALSE]
        A <- A[, -iM, drop = FALSE]
    }

    if (is.null(row_xss)) row_xss <- rowSums(X * X)
    AAt <- tcrossprod(A)
    XAt <- tcrossprod(X, A)

    row_rss <- .aa_trace_row_rss(row_xss, S, XAt, AAt)
    row_weights <- weight_fun(row_rss)
    .aa_check_row_weights(row_weights, nrow(X))
    if (.aa_trivial_row_weights(row_weights))
        row_weights <- NULL

    S_weighted <- .aa_weight_rows(S, row_weights)
    list(
        rss = sum(.aa_weight_rows(row_rss, row_weights)),
        xss = sum(.aa_weight_rows(row_xss, row_weights)),
        row_xss = row_xss,
        row_rss = row_rss,
        row_weights = row_weights,
        StS = crossprod(S_weighted, S),
        S_weighted = S_weighted,
        AAt = AAt,
        XAt = XAt,
        A = A
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
                         loss_fun,
                         A,
                         B,
                         S,
                         loss,
                         ols_solver,
                         max_no_update) {
    # Nomenclature following arXiv:2504.12392v1:
    #   X ~ SA (N x M) Data Matrix
    #   A = BX (K x M) Archetypes
    #   B = (K x N) Archetypes Coefficients (base transform, C in the paper)
    #   S = (N x K) Archetypes Scores (new coordinates)
    #   rss = ||X - SA||^2 = ||X||^2 - 2*tr(SAXt) + tr(StS AAt) Residual Sum of Squares
    A0 <- A
    Xt <- t(X)  # compute Xt once to reuse for B update
    nnls_svd_kappa_threshold <- 500
    loss_terms <- loss_fun(X, A, S, weight_fun = weight_fun)
    xss <- loss_terms[["xss"]]
    row_xss <- loss_terms[["row_xss"]]
    row_weights <- loss_terms[["row_weights"]]
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
    no_update <- 0L
    max_simplex_error <- 0
    project <- proj_l1
    best_loss <- loss_terms[["rss"]]
    best_args <- list(A = A, B = B, S = S)

    # edge case: if max_iter = 0 return initial solution without any updates
    if (max_iter == 0L) {
        return(list(
            A0 = A0,
            A = A,
            B = B,
            S = S,
            delta = 0,
            i = 0L,
            loss = loss,
            converged = TRUE
        ))
    }


    # Main Loop  --------------------------------------------------------------

    if (verbose) message("Starting main loop...")
    for (i in seq_len(max_iter)) {
        j <- i + 1L  # loss row to update
        check_kappa <- i %% 10L == 0L  # Check kappa every 10 iterations

        # S update
        S_raw <- fit_nnls(X, t(A), use_svd = use_svd_for_S) # Project X to A-simplex
        max_simplex_error <- max(max_simplex_error, max(abs(rowSums(S_raw) - 1)))
        S <- project(S_raw, eps = eps)
        # A update
        A <- fit_ols(S, X, method = ols_solver, row_weights = row_weights)
        # B update
        B_raw <- fit_nnls(A, Xt, use_svd = FALSE) # Project A to X-simplex
        max_simplex_error <- max(max_simplex_error, max(abs(rowSums(B_raw) - 1)))
        B <- project(B_raw, eps = eps)
        # Final A update to ensure A = BX
        A <- B %*% X

        # Update objective
        loss_terms <- loss_fun(
            X,
            A,
            S,
            weight_fun = weight_fun,
            xss = xss,
            row_xss = row_xss
        )
        row_weights <- loss_terms[["row_weights"]]
        loss <- .aa_update_loss(
            loss,
            j,
            loss_terms,
            verbose = verbose,
            max_kappa = max_kappa
        )

        if (loss_terms[["rss"]] < best_loss) {
            best_loss <- loss_terms[["rss"]]
            best_args <- list(A = A, B = B, S = S)
            no_update <- 0L
        } else {
            no_update <- no_update + 1L

            if (verbose) {
                fmt <- paste(
                    "Iteration %d: NNLS update did not improve best loss;",
                    "continuing from current iterate (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                fmt <- paste(
                    "NNLS stalled: no loss improvement after",
                    "%d consecutive candidate updates"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }
        }

        # Check if A is ill-conditioned and if we should switch to SVD for S update
        if (check_kappa) {
            k_A <- loss[["k_A"]][j]
            # TODO: exact kappa already computes the SVD? maybe we should resuse it.
            if (is.na(k_A))
                loss[["k_A"]][j] <- k_A <- kappa(A, exact = TRUE)
            use_svd_for_S <- k_A > nnls_svd_kappa_threshold
        }

        # Check convergence
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }
    if (max_simplex_error > 0.05) {
        fmt <- paste(
            "Raw NNLS coefficients were not close to simplex constraints;",
            "maximum absolute row-sum error was %.3g. Consider increasing `bigM`."
        )
        warning(sprintf(fmt, max_simplex_error), call. = FALSE)
    }

    c(best_args, list(delta = 0, A0 = A0, i = i, loss = loss, converged = converged))
}

# Fit Non-negative Least Squares for every row of Y
#
# This function solves the problem: $\min_{B} ||Y - BX||_2 s.t. B >= 0$
#
# @param Y data matrix (rows = samples, columns = dimensions)
# @param X data matrix (rows = samples, columns = dimensions)
# @param use_svd logical, whether to use SVD for dimensionality reduction (default: FALSE)
fit_nnls <- function(Y, X, use_svd = FALSE) {
    # min ||Y - Beta %*% X||_2 s.t. Beta >= 0

    if (use_svd) {
        s <- svd(X)
        # TODO: choose rank based on explained variance
        X <- with(s, t(v * d))
        Y <- Y %*% s$u
    }
    # TODO: parallelize
    Beta <- matrix(0, nrow = nrow(Y), ncol = ncol(X))
    for (i in seq_len(nrow(Y)))
        Beta[i, ] <- stats::coef(nnls::nnls(X, Y[i, ]))
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
    if (!.aa_trivial_row_weights(row_weights)) {
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
        norm(R, "F")^2
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
