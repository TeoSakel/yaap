#' Archetypes Analysis using Projected Gradient Descent
#'
#' @param data data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param init function, method string, or numeric coordinate matrix to initialize
#'   archetypes (default: `"furthest_sum"`). When a matrix is supplied it must
#'   have dimension `K x ncol(data)`. Rows outside the `delta`-relaxed convex
#'   hull of `data` are projected into it with a warning; row names, when
#'   present, are used as archetype names.
#' @param init_args list of additional arguments for the initialization function
#' @param weights optional vector of sample weights (default: NULL)
#' @param scale scaling or metric embedding used before fitting. `TRUE` applies
#'   the default z-score preprocessing, `FALSE` leaves columns on their original
#'   scale, a positive numeric vector divides columns by user-supplied scale
#'   factors, and a symmetric positive-definite matrix applies the corresponding
#'   feature metric embedding in the original data column space.
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
#' @param missing whether to fit the missing-data PGD objective. When `TRUE`,
#'   only observed entries are optimized; dense `NA` values are treated as
#'   missing and sparse structural zeros are treated as missing.
#' @param delta maximum allowed relaxation of archetypes convexity constraint (default: 0)
#' @param pseudo_pgd whether to use pseudo projected gradient descent (default: TRUE)
#' @param step_size initial step size for the gradient descent (default: 1.0)
#' @param max_iter_optimizer maximum iterations for the line search optimizer (default: 10)
#' @param step_shrinkage factor to reduce step size if no improvement during line search (default: 0.5)
#' @param max_no_update maximum consecutive outer iterations with no accepted
#'   line-search update before PGD stops as stalled (default: 5)
#'
#' @returns An object of class \code{\link{archetypes}}
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "YAAAP"))
#' archetypes_pgd(as.matrix(toy), K = 3)
#'
#' @references
#' Mørup, M., & Hansen, L. K. (2012).
#' Archetypal analysis for machine learning and data mining.
#' *Neurocomputing*, 80, 54-63. \url{https://dx.doi.org/10.1016/j.neucom.2011.06.033}
#'
#' @export
archetypes_pgd <- function(data,
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
                           missing = any(is.na(data)),
                           # PGD specific
                           delta = 0,
                           pseudo_pgd = TRUE,
                           step_size = 1.0,
                           max_iter_optimizer = 10L,
                           step_shrinkage = 0.5,
                           max_no_update = 5L) {
    .aa_run_aa_default(
        call = match.call(),
        data = data,
        K = K,
        method = "pgd",
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
        missing = missing,
        delta = delta,
        pseudo_pgd = pseudo_pgd,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = max_no_update
    )
}

.aa_fit_pgd_missing <- function(X,
                                M,
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
                                delta,
                                pseudo_pgd,
                                step_size,
                                max_iter_optimizer,
                                step_shrinkage,
                                max_no_update) {
    if (!is.null(weight_fun))
        stop("`robust = TRUE` is not supported with `missing = TRUE`.", call. = FALSE)

    if (pseudo_pgd) {
        project <- proj_l1
    } else {
        project <- proj_simplex
    }

    update_alpha <- delta > 0
    a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
    a_hi <- 1 + delta
    clip <- function(a) pmax(pmin(a, a_hi), a_lo)
    denom_eps <- ifelse(eps > 0, eps, 1e-8)

    A0 <- A
    a <- rowSums(B)
    slack_tol <- 1e-6
    if (any(a < a_lo - slack_tol) || any(a > a_hi + slack_tol)) {
        fmt <- "Initialize B marginals are outside the specified delta range [%.3f, %.3f]"
        stop(sprintf(fmt, a_lo, a_hi))
    }
    H <- B
    B <- B / a

    A <- .aa_pgd_missing_A(H, X, M, denom_eps)
    R <- .aa_pgd_missing_resid(M, S, A, X)  # residuals matrix
    xss <- norm(X, "F")^2  # X sum of squares; invariant for fixed X
    rss <- norm(R, "F")^2  # residual sum of squares -> objective to minimize
    loss_terms <- list(
        rss = rss,
        xss = xss,
        StS = crossprod(S),
        A = A
    )
    loss <- .aa_update_loss(
        loss,
        1L,
        loss_terms,
        verbose = verbose,
        max_kappa = max_kappa
    )
    converged <- FALSE

    step_S <- step_size
    step_B <- step_size
    step_a <- step_size
    no_update <- 0L

    if (verbose) message("Starting missing-data optimization loop...")
    for (i in seq_len(max_iter)) {
        accepted_update <- FALSE

        # Update S
        grad <- .aa_pgd_missing_grad_S(S, R, A, pseudo_pgd)
        for (k in seq_len(max_iter_optimizer)) {
            S_new <- project(S - step_S * grad, eps = eps)
            R_new <- .aa_pgd_missing_resid(M, S_new, A, X)
            rss_new <- norm(R_new, "F")^2
            if (rss_new < rss) {
                S <- S_new
                R <- R_new
                rss <- rss_new
                step_S <- step_S / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }

        # Update B
        grad_H <- .aa_pgd_missing_grad_H(H, X, M, S, A, R, denom_eps)
        grad <- if (pseudo_pgd) {
            grad_H - rowSums(B * grad_H)
        } else {
            grad_H
        }
        grad_B_for_alpha <- grad_H
        grad <- a * grad
        for (k in seq_len(max_iter_optimizer)) {
            B_new <- project(B - step_B * grad, eps = eps)
            H_new <- a * B_new
            A_new <- .aa_pgd_missing_A(H_new, X, M, denom_eps)
            R_new <- .aa_pgd_missing_resid(M, S, A_new, X)
            rss_new <- norm(R_new, "F")^2
            if (rss_new < rss) {
                B <- B_new
                H <- H_new
                A <- A_new
                R <- R_new
                rss <- rss_new
                step_B <- step_B / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_B <- step_B * step_shrinkage
        }

        # Update a:
        if (update_alpha) {
            grad <- grad_alpha(B, grad_B_for_alpha)
            for (k in seq_len(max_iter_optimizer)) {
                a_new <- clip(a - step_a * grad)
                H_new <- a_new * B
                A_new <- .aa_pgd_missing_A(H_new, X, M, denom_eps)
                R_new <- .aa_pgd_missing_resid(M, S, A_new, X)
                rss_new <- norm(R_new, "F")^2
                if (rss_new < rss) {
                    a <- a_new
                    H <- H_new
                    A <- A_new
                    R <- R_new
                    rss <- rss_new
                    step_a <- step_a / step_shrinkage
                    accepted_update <- TRUE
                    break
                }
                step_a <- step_a * step_shrinkage
            }
        }

        # Update loss terms
        loss_terms <- list(
            rss = rss,
            xss = xss,
            StS = if (i %% 10 == 0 && max_kappa > 1) crossprod(S) else NULL,
            A = A
        )
        loss <- .aa_update_loss(
            loss,
            i + 1L,
            loss_terms,
            verbose = verbose,
            max_kappa = max_kappa
        )

        # Check update acceptance and convergence
        if (!accepted_update) {
            no_update <- no_update + 1L
            step_S <- step_S * step_shrinkage
            step_B <- step_B * step_shrinkage
            if (update_alpha)
                step_a <- step_a * step_shrinkage

            if (verbose) {
                fmt <- paste(
                    "Iteration %d: no line-search update accepted;",
                    "shrinking steps (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                fmt <- paste(
                    "Missing-data PGD stalled: no line-search update accepted after",
                    "%d consecutive step shrinkages"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }
            next
        }

        no_update <- 0L
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }

    list(
        A0 = A0,
        A = A,
        B = H,
        S = S,
        delta = delta,
        i = i,
        loss = loss,
        converged = converged
    )
}

.aa_fit_pgd <- function(X,
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
                        delta,
                        pseudo_pgd,
                        step_size,
                        max_iter_optimizer,
                        step_shrinkage,
                        max_no_update) {
    # Nomenclature following www.doi.org/10.1016/j.neucom.2011.06.033:
    #   X ~ SA (N x M) Data Matrix
    #   A = aBX (K x M) Archetypes
    #   B = (K x N) Archetypes Coefficients (base transform, C in the paper)
    #   a = (K x 1) Archetypes scaling (allows them to lay outside convex hull)
    #               In theory they are stored in diagonal matrix.
    #   S = (N x K) Archetypes Scores (new coordinates)

    # PGD specific preparations -----------------------------------------------

    # define gradient & projection functions
    if (pseudo_pgd) {
        grad_S  <- grad_S_l1
        grad_B  <- grad_B_l1
        project <- proj_l1
    } else {
        grad_S  <- grad_S_simplex
        grad_B  <- grad_B_simplex
        project <- proj_simplex
    }

    # Setup alpha updates
    update_alpha <- delta > 0
    a_lo <- ifelse(eps > 0, eps, 1e-8) # to avoid division by 0
    a_lo <- max(1 - delta, a_lo) # lower bound for a
    a_hi <- 1 + delta # upper bound for a
    clip <- function(a) pmax(pmin(a, a_hi), a_lo) # clip a to [a_lo, a_hi]

    A0 <- A
    a  <- rowSums(B)

    slack_tol <- 1e-6
    if (any(a < a_lo - slack_tol) || any(a > a_hi + slack_tol)) {
        fmt <- "Initialize B marginals are outside the specified delta range [%.3f, %.3f]"
        stop(sprintf(fmt, a_lo, a_hi))
    }
    B  <- B / a # normalize B to row-stochastic

    # Compute auxiliary variables
    A     <- A0                # (K x M) = Archetypes = a*B %*% X
    loss_terms <- .aa_loss_terms(
        X, A, S,
        weight_fun = weight_fun,
        return_S_terms = TRUE
    )
    XAt <- loss_terms[["XAt"]] # (N x K)
    AAt <- loss_terms[["AAt"]] # (K x K)
    StS <- loss_terms[["StS"]] # (K x K)
    StX <- loss_terms[["StX"]] # (K x M)
    xss <- loss_terms[["xss"]] # scalar
    rss <- loss_terms[["rss"]] # scalar
    XXt <- tcrossprod(X)       # (N x N)
    row_xss <- loss_terms[["row_xss"]]  # (N x 1) row-wise sum of squares; invariant for fixed X
    row_weights <- loss_terms[["row_weights"]]
    S_weighted <- loss_terms[["S_weighted"]]
    StXXt <- crossprod(S_weighted, XXt) # (K x N)

    loss <- .aa_update_loss(
        loss,
        1L,
        loss_terms,
        verbose = verbose,
        max_kappa = max_kappa
    )
    converged <- FALSE

    ## Step sizes
    step_S <- step_size # mu_S in the paper
    step_B <- step_size # mu_C in the paper
    step_a <- step_size # mu_a in the paper
    no_update <- 0L

    # Main optimization loop  -------------------------------------------------

    if (verbose) message("Starting optimization loop...")
    for (i in seq_len(max_iter)) {
        accepted_update <- FALSE

        ## Update S
        grad <- grad_S(S, AAt, XAt, row_weights = row_weights) # (N x K) - diag(a) is absorbed in A
        for (k in seq_len(max_iter_optimizer)) { # line search
            S_new   <- project(S - step_S * grad, eps = eps)
            S_new_weighted <- .aa_weight_rows(S_new, row_weights)
            StS_new <- crossprod(S_new_weighted, S_new)
            StX_new <- crossprod(S_new_weighted, X)
            rss_new <- xss - 2 * sum(A * StX_new) + sum(StS_new * AAt)
            if (rss_new < rss) { # update variables
                S      <- S_new
                StS    <- StS_new
                StX    <- StX_new
                rss    <- rss_new
                StXXt  <- crossprod(S_new_weighted, XXt)
                step_S <- step_S / step_shrinkage # leave room for shrinkage
                accepted_update <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }

        ## Update B
        grad <- a * grad_B(B, A, X, StS, StXXt) # (K x N)
        for (k in seq_len(max_iter_optimizer)) { # line search
            B_new <- project(B - step_B * grad, eps = eps)
            A_new <- (a * B_new) %*% X
            AAt_new <- tcrossprod(A_new)
            rss_new <- xss - 2 * sum(A_new * StX) + sum(StS * AAt_new)
            if (rss_new < rss) {
                B   <- B_new
                A   <- A_new
                rss <- rss_new
                AAt <- AAt_new
                XAt <- tcrossprod(X, A)
                step_B <- step_B / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_B <- step_B * step_shrinkage
        }

        ## Update a:
        if (update_alpha) {
            grad <- grad_alpha(B, grad / a)
            for (k in seq_len(max_iter_optimizer)) {
                a_new    <- clip(a - step_a * grad)
                a_update <- a_new / a # to undo aA multiplication
                A_new    <- a_update * A
                AAt_new  <- AAt * tcrossprod(a_update) # fast outer product
                rss_new  <- xss - 2 * sum(A_new * StX) + sum(StS * AAt_new)
                if (rss_new < rss) { # update variables
                    a   <- a_new
                    A   <- A_new
                    rss <- rss_new
                    AAt <- AAt_new
                    XAt <- sweep(XAt, 2, a_update, "*")
                    step_a <- step_a / step_shrinkage # leave room for shrinkage
                    accepted_update <- TRUE
                    break
                }
                step_a <- step_a * step_shrinkage
            }
        }

        loss_terms <- .aa_loss_terms(X, A, S, weight_fun = weight_fun, return_S_terms = TRUE,
                                     xss = xss, rss = rss, row_xss = row_xss,
                                     row_weights = row_weights,
                                     StS = StS, StX = StX, AAt = AAt, XAt = XAt)
        row_weights <- loss_terms[["row_weights"]]
        AAt   <- loss_terms[["AAt"]]
        StS   <- loss_terms[["StS"]]
        StX   <- loss_terms[["StX"]]
        xss   <- loss_terms[["xss"]]
        rss   <- loss_terms[["rss"]]
        S_weighted <- loss_terms[["S_weighted"]]
        StXXt <- crossprod(S_weighted, XXt)

        # Check convergence
        loss <- .aa_update_loss(
            loss,
            i + 1L,
            loss_terms,
            verbose = verbose,
            max_kappa = max_kappa
        )
        if (!accepted_update) {
            no_update <- no_update + 1L
            step_S <- step_S * step_shrinkage
            step_B <- step_B * step_shrinkage
            if (update_alpha)
                step_a <- step_a * step_shrinkage

            if (verbose) {
                fmt <- paste(
                    "Iteration %d: no line-search update accepted;",
                    "shrinking steps (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                fmt <- paste(
                    "PGD stalled: no line-search update accepted after",
                    "%d consecutive step shrinkages"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }
            next
        }

        no_update <- 0L
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }


    list(
        A0 = A0,
        A = A,
        B = a * B,
        S = S,
        delta = delta,
        i = i,
        loss = loss,
        converged = converged
    )
}

# Gradient functions for RSS-normed archetypes analysis
# ||X - SBX||^2 = tr(XXt - 2SBXXt + SBXXBtSt)
grad_S_simplex <- function(S, AAt, XAt, row_weights = NULL) {
    # 2 * (SBXXtBt - XXtBt) = 2 * (SAAt - XAt)
    .aa_weight_rows(S %*% AAt - XAt, row_weights)
}

grad_B_simplex <- function(B, A, X, StS, StXXt) {
    # 2 * (StSBXXt - StXXt) = 2 * (StSAXt - StXXt)
    tcrossprod(StS %*% A, X) - StXXt
}

grad_S_l1 <- function(S, ...) {
    grad <- grad_S_simplex(S, ...) # (N x K)
    row_dot <- rowSums(S * grad)
    grad - row_dot
}

grad_B_l1 <- function(B, ...) {
    grad <- grad_B_simplex(B, ...) # (K x N)
    row_dot <- rowSums(B * grad)
    grad - row_dot
}

grad_alpha <- function(B, grad_aB) {
    # dR/da = dR/d(aB) * d(aB)/da = dR/d(aB) * B
    rowSums(grad_aB * B)
}

.aa_pgd_missing_A <- function(H, X, M, eps) {
    numerator <- H %*% X
    denominator <- H %*% M + eps
    as.matrix(numerator / denominator)
}

.aa_pgd_missing_resid <- function(M, S, A, X) {
    if (!inherits(M, "sparseMatrix")) {
        Xhat <- S %*% A
        Xhat[!M] <- 0
        return(Xhat - X)
    }

    entries <- Matrix::summary(M)  # i, j, x for nonzero entries of M
    if (length(entries[["i"]]) == 0L)
        return(Matrix::sparseMatrix(dims = dim(M), dimnames = dimnames(M)))
    vals <- rowSums(S[entries[["i"]], , drop = FALSE] * t(A[, entries[["j"]], drop = FALSE]))
    Xhat <- Matrix::sparseMatrix(
        i = entries[["i"]],
        j = entries[["j"]],
        x = vals,
        dims = dim(M),
        dimnames = dimnames(M)
    )
    Xhat - X
}

.aa_pgd_missing_grad_S <- function(S, R, A, pseudo_pgd) {
    grad <- R %*% t(A)
    if (!pseudo_pgd)
        return(grad)
    grad - rowSums(S * grad)
}

.aa_pgd_missing_grad_H <- function(H, X, M, S, A, R, eps) {
    grad <- crossprod(R, S) / (t(H %*% M) + eps)
    t(X %*% grad - M %*% (grad * t(A)))
}
