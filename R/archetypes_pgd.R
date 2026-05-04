#' Archetypes Analysis using Projected Gradient Descent
#'
#' @param data data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param init function or method string to initialize archetypes
#'   (default: `"furthest_sum"`)
#' @param init_args list of additional arguments for the initialization function
#' @param weights optional vector of sample weights (default: NULL)
#' @param sd_threshold threshold for feature standard deviation to filter low-variance features (default: 1e-4)
#' @param max_iter maximum number of iterations (default: 500)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-6)
#' @param tol_r2 convergence tolerance based on R^2 (default: 0.9999)
#' @param max_kappa maximum condition number for archetypes (default: 1000)
#' @param eps small positive number to ensure numerical stability
#'   (default: 0 for sparse input 1e-8 for dense)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param delta maximum allowed relaxation of archetypes convexity constraint (default: 0)
#' @param pseudo_pgd whether to use pseudo projected gradient descent (default: TRUE)
#' @param step_size initial step size for the gradient descent (default: 1.0)
#' @param max_iter_optimizer maximum iterations for the line search optimizer (default: 10)
#' @param step_shrinkage factor to reduce step size if no improvement during line search (default: 0.5)
#'
#' @returns An object of class \code{\link{archetypes}}
#'
#' @examples
#' data(toy)
#' archetypes_pgd(toy, K = 3)
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
                           sd_threshold = 1e-4,
                           max_iter = 500L,
                           tol = 1e-6,
                           tol_r2 = 0.9999,
                           max_kappa = 1000,
                           eps = ifelse(is(data, "sparseMatrix"), 0, 1e-8),
                           verbose = FALSE,
                           # PGD specific
                           delta = 0,
                           pseudo_pgd = TRUE,
                           step_size = 1.0,
                           max_iter_optimizer = 10L,
                           step_shrinkage = 0.5) {

    # Input Check -------------------------------------------------------------

    # Generic Checks
    .aa_check_inputs(data=data, K=K, tol=tol, tol_r2=tol_r2, max_kappa=max_kappa, eps=eps)
    # PGD specific checks
    stopifnot("step_size must be positive" = step_size > 0)
    stopifnot("step_shrinkage must be between (0, 1)" = step_shrinkage > 0 && step_shrinkage < 1)
    stopifnot("delta must be single non-negative number" = length(delta) == 1 && delta >= 0)

    # Edge Case checks
    out <- .aa_checks_edge_cases(data, K, verbose)  # edge cases
    if (!is.null(out)) return(out)  # return early if edge case

    # Prepossessing Data  -----------------------------------------------------

    # Nomenclature following www.doi.org/10.1016/j.neucom.2011.06.033:
    #   X ~ SA (N x M) Data Matrix
    #   A = aBX (K x M) Archetypes
    #   B = (K x N) Archetypes Coefficients (base transform, C in the paper)
    #   a = (K x 1) Archetypes scaling (allows them to lay outside convex hull)
    #               In theory they are stored in diagonal matrix.
    #   S = (N x K) Archetypes Scores (new coordinates)

    cl <- match.call()
    pre <- .aa_preprocess(data, sd_threshold, weights, verbose)
    X <- pre[["X"]]                    # preprocessed data
    # N <- nrow(X)
    undo_scale <- pre[["undo_scale"]]  # function to undo preprocessing
    xss <- pre[["xss"]]                # total sum of squares
    rm(pre)

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

    # Initialization  ---------------------------------------------------------

    init_vars <- .aa_init_vars(X, K, init, init_args, eps, max_iter, verbose)
    A0 <- init_vars[["A"]]
    B  <- init_vars[["B"]]
    S  <- init_vars[["S"]]
    a  <- rowSums(B)

    if (any(a < a_lo) || any(a > a_hi)) {
        fmt <- "Initialize B marginals are outside the specified delta range [%.3f, %.3f]"
        stop(sprintf(fmt, a_lo, a_hi))
    }
    B  <- B / a # normalize B to row-stochastic

    # Compute auxiliary variables
    A     <- A0                # (K x M) = Archetypes = a*B %*% X
    AAt   <- tcrossprod(A)     # (K x K)
    XAt   <- tcrossprod(X, A)  # (N x K)
    StS   <- crossprod(S)      # (K x K)
    StX   <- crossprod(S, X)   # (K x M)
    XXt   <- tcrossprod(X)     # (N x N)
    StXXt <- crossprod(S, XXt) # (K x N)

    # Loss
    loss <- init_vars[["loss"]]
    ### RSS: ||X - SA||^2 = ||X||^2 - 2 * tr(X^T SA) + tr(SA A^T S^T)
    ### Reorganize the terms inside trace to minimize memory footprint
    rss <- xss - 2 * sum(A * StX) + sum(StS * AAt)
    loss <- .aa_update_loss(
        loss,
        1L,
        verbose = verbose,
        max_kappa = max_kappa,
        rss = rss,
        xss = xss,
        StS = StS,
        AAt = AAt,
        A = A
    )
    converged <- FALSE

    ## Step sizes
    step_S <- step_size # mu_S in the paper
    step_B <- step_size # mu_C in the paper
    step_a <- step_size # mu_a in the paper

    rm(init_vars)

    # Main optimization loop  -------------------------------------------------

    if (verbose) message("Starting optimization loop...")
    for (i in seq_len(max_iter)) {
        ## Update S
        grad <- grad_S(S, AAt, XAt) # (N x K) - diag(a) is absorbed in A
        for (k in seq_len(max_iter_optimizer)) { # line search
            S_new   <- project(S - step_S * grad, eps = eps)
            StS_new <- crossprod(S_new)
            StX_new <- crossprod(S_new, X)
            rss_new <- xss - 2 * sum(A * StX_new) + sum(StS_new * AAt)
            if (rss_new < rss) { # update variables
                S      <- S_new
                StS    <- StS_new
                StX    <- StX_new
                rss    <- rss_new
                StXXt  <- crossprod(S, XXt)
                step_S <- step_S / step_shrinkage # leave room for shrinkage
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
                    break
                }
                step_a <- step_a * step_shrinkage
            }
        }


        # Check convergence
        loss <- .aa_update_loss(
            loss,
            i + 1L,
            verbose = verbose,
            max_kappa = max_kappa,
            rss = rss,
            xss = xss,
            StS = StS,
            AAt = AAt,
            A = A
        )
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }


    # Prepare Output  ---------------------------------------------------------

    out <- .aa_prepare_output(
        call = cl,
        data = data,
        X = X,
        A0 = A0,
        A = A,
        B = a*B,
        S = S,
        delta = delta,
        i = i,
        loss = loss,
        converged = converged,
        undo_scale = undo_scale,
        max_iter = max_iter,
        verbose = verbose
    )
    return(out)
}

# Gradient functions for RSS-normed archetypes analysis
# ||X - SBX||^2 = tr(XXt - 2SBXXt + SBXXBtSt)
grad_S_simplex <- function(S, AAt, XAt) {
    # 2 * (SBXXtBt - XXtBt) = 2 * (SAAt - XAt)
    return(S %*% AAt - XAt)
}

grad_B_simplex <- function(B, A, X, StS, StXXt) {
    # 2 * (StSBXXt - StXXt) = 2 * (StSAXt - StXXt)
    return(tcrossprod(StS %*% A, X) - StXXt)
}

grad_S_l1 <- function(S, ...) {
    grad <- grad_S_simplex(S, ...) # (N x K)
    row_dot <- rowSums(S * grad)
    return(grad - row_dot)
}

grad_B_l1 <- function(B, ...) {
    grad <- grad_B_simplex(B, ...) # (K x N)
    row_dot <- rowSums(B * grad)
    return(grad - row_dot)
}

grad_alpha <- function(B, grad_aB) {
    # dR/da = dR/d(aB) * d(aB)/da = dR/d(aB) * B
    return(rowSums(grad_aB * B))
}
