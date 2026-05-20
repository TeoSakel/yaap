#' Archetypes Analysis using Projected Gradient Descent
#'
#' Fits an archetypal analysis (AA) model by minimising the reconstruction error
#' \eqn{\|X - SA\|_F^2}, where \eqn{S} is the composition matrix and \eqn{A} is
#' the archetype coordinates via projected gradient descent (PGD).
#'
#' @param x data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param init initialization method; see [run_aa()] for available options.
#' @param init_args list of additional arguments for the initialization function.
#' @param weights optional vector of sample weights (default: NULL)
#' @param scale scaling or metric embedding used before fitting. `FALSE`
#'   (default) leaves columns on their original scale, `TRUE` applies
#'   z-score preprocessing, a positive numeric vector divides columns by
#'   user-supplied scale factors, and a symmetric positive-definite matrix
#'   applies the corresponding feature metric embedding in the original data
#'   column space.
#' @param robust robust row reweighting selector. Use `FALSE` for ordinary
#'   squared error, `TRUE` for `"psi.bisquare"`, a MASS psi function name,
#'   or a custom psi function. See [MASS::rlm()] for psi details; `method =
#'   "MM"` is not supported because it is not applicable to AA.
#' @param robust_args list of tuning arguments passed to the robust psi function.
#' @param sd_threshold threshold for feature standard deviation to filter
#'   low-variance features (default: 1e-6)
#' @param max_iter maximum number of iterations (default: 100)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-6)
#' @param tol_r2 convergence tolerance based on R^2 (default: 0.9999)
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
#' @details
#' ## Model family
#'
#' Three arguments change *which* objective is optimised and therefore the
#' statistical family of the fitted model:
#'
#' * **`delta`** relaxes the convexity constraint on the archetypes.  With the
#'   default `delta = 0`, every archetype is a convex combination of observed
#'   data points (i.e., it lies inside or on the boundary of the data convex
#'   hull).  Setting `delta > 0` allows archetypes to step outside the hull by
#'   up to a factor \eqn{1 + \delta} and inside by up to \eqn{1 - \delta} —
#'   useful when the true extremes are absent from the data due to truncation or
#'   censoring.  Use with care: large values remove the interpretability
#'   guarantee that archetypes are representable as extreme prototypes, and the
#'   uniqueness result of Theorem 1 in Mørup & Hansen (2012) no longer holds.
#'
#' * **`robust`** switches the loss from ordinary squared error to an
#'   iteratively re-weighted version based on MASS-style psi row weights.
#'   Use `TRUE` for `"psi.bisquare"`, or pass `"psi.huber"`, `"psi.hampel"`,
#'   or a custom psi function with tuning values in `robust_args`. See
#'   [MASS::rlm()] for psi details. `method = "MM"` from `MASS::rlm()` is not
#'   supported because MM estimation is not directly applicable to the AA
#'   objective. Robust fitting is incompatible with `missing = TRUE`.
#'
#' * **`missing`** activates the missing-data objective described in Section 3.5
#'   of Mørup & Hansen (2012).  When `TRUE`, only observed entries contribute to
#'   the loss and the archetype profiles are normalised entry-wise so that each
#'   archetype is a weighted average of the *observed* values of the data points
#'   that define it.  For dense matrices, `NA` marks missing entries; for sparse
#'   \code{sparseMatrix} inputs, structural zeros are treated as missing.  The
#'   default is `any(is.na(x))`, so missing-data mode is activated automatically
#'   when the input contains `NA`s.  For kernel and probabilistic variants of AA
#'   see \code{vignette("non-gaussian-aa", package = "yaap")}.
#'
#' ## Algorithm mechanics and tuning
#'
#' The solver alternates gradient steps for \eqn{S} and \eqn{B} with an inner
#' line search of up to `max_iter_optimizer` steps (default 10) that shrinks the
#' step size by `step_shrinkage` (default 0.5) until the objective decreases.
#' The outer loop repeats for at most `max_iter` iterations (default 100).
#'
#' When `pseudo_pgd = TRUE` (the default), the gradients are projected onto the
#' tangent space of the \eqn{\ell_1} unit sphere before the simplex projection
#' step.  This is the "pseudo-PGD" variant of Mørup & Hansen (2012), which
#' removes the radial gradient component and thereby prevents the step from
#' increasing the \eqn{\ell_1} norm before projection.  In practice it leads to
#' larger effective steps and faster convergence; set `pseudo_pgd = FALSE` to
#' use the plain projected gradient if needed.
#'
#' `step_size` sets the initial step size for both \eqn{S} and \eqn{B} updates
#' (and for \eqn{\alpha} when `delta > 0`).  If no step in the inner line
#' search is accepted for `max_no_update` consecutive outer iterations (default
#' 5), the solver declares the iterate stalled and stops early.  Step sizes are
#' also reduced after a failed outer iteration, so in practice the solver
#' self-tunes step sizes down as it approaches a local minimum; there is usually
#' no need to change `step_size` unless the default is far from the optimal
#' scale for a particular dataset.
#'
#' Convergence is declared when the relative decrease in RSS between successive
#' accepted updates falls below `tol` (default 1e-6) *or* \eqn{R^2} exceeds
#' `tol_r2` (default 0.9999).
#'
#' ## Preprocessing
#'
#' Before fitting, features with standard deviation below `sd_threshold`
#' (default 1e-6) are dropped to avoid ill-conditioning; their values are
#' restored as column means in the output.  The `scale` argument controls the
#' feature metric: `FALSE` (default) leaves columns on their original scale, `TRUE`
#' applies column-wise z-scoring, a positive numeric vector divides
#' each column by the supplied factor, and a symmetric positive-definite matrix
#' defines a full quadratic feature metric — see
#' \code{vignette("non-gaussian-aa", package = "yaap")} for the metric-AA use
#' case.  Sample-level `weights` re-weight rows in the loss but do not affect
#' the preprocessing step.
#'
#' @returns An object of class \code{archetypes}
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' archetypes_pgd(as.matrix(toy), K = 3)
#'
#' @references
#' Mørup, M., & Hansen, L. K. (2012).
#' Archetypal analysis for machine learning and data mining.
#' *Neurocomputing*, 80, 54-63. \doi{10.1016/j.neucom.2011.06.033}
#'
#' @export
archetypes_pgd <- function(x,
                           K,
                           init = "furthest_sum",
                           init_args = list(),
                           weights = NULL,
                           scale = FALSE,
                           robust = FALSE,
                           robust_args = list(),
                           sd_threshold = 1e-6,
                           max_iter = 100L,
                           tol = 1e-6,
                           tol_r2 = 0.9999,
                           eps = ifelse(inherits(x, "sparseMatrix"), 0, 1e-8),
                           verbose = FALSE,
                           missing = any(is.na(x)),
                           # PGD specific
                           delta = 0,
                           pseudo_pgd = TRUE,
                           step_size = 1.0,
                           max_iter_optimizer = 10L,
                           step_shrinkage = 0.5,
                           max_no_update = 5L) {
    .aa_fit_engine(
        call = match.call(),
        x = x,
        K = K,
        method = "pgd",
        init = init,
        init_args = init_args,
        weights = weights,
        scale = scale,
        robust = robust,
        robust_args = robust_args,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
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

.aa_pgd_block <- function(ctx,
                          delta = 0,
                          pseudo_pgd = TRUE,
                          step_size = 1.0,
                          max_iter_optimizer = 10L,
                          step_shrinkage = 0.5,
                          max_no_update = 5L) {
    loss_fun <- if (!identical(ctx[["robust"]], FALSE)) .aa_pgd_weighted_loss else .aa_pgd_loss
    fit_fun <- if (ctx[["missing"]]) .aa_fit_pgd_missing else .aa_fit_pgd

    list(
        check = function(ctx) {
            .aa_euclidean_check(ctx)
            .aa_check_projected_gradient_controls(
                step_size = step_size,
                max_iter_optimizer = max_iter_optimizer,
                step_shrinkage = step_shrinkage,
                max_no_update = max_no_update
            )
            if (!is_non_negative(delta)) {
                stop("`delta` must be a single non-negative number.", call. = FALSE)
            }
            if (!is_logical(pseudo_pgd)) {
                stop("`pseudo_pgd` must be TRUE or FALSE.", call. = FALSE)
            }
            invisible(TRUE)
        },
        preprocess = function(ctx) .aa_euclidean_preprocess(ctx, bigM = 0),
        edge_case = .aa_euclidean_edge_case,
        init = function(ctx, prep) .aa_euclidean_init(ctx, prep, delta = delta),
        fit = function(ctx, prep, init_vars) {
            common_args <- list(
                X = prep[["X"]],
                weight_fun = .aa_weight_fun(ctx[["robust"]], ctx[["robust_args"]]),
                max_iter = ctx[["max_iter"]],
                tol = ctx[["tol"]],
                tol_r2 = ctx[["tol_r2"]],
                eps = ctx[["eps"]],
                verbose = ctx[["verbose"]]
            )
            if (ctx[["missing"]]) {
                common_args[["M"]] <- prep[["M"]]
            } else {
                common_args[["loss_fun"]] <- loss_fun
            }
            do.call(
                fit_fun,
                c(
                    common_args,
                    init_vars,
                    list(
                        delta = delta,
                        pseudo_pgd = pseudo_pgd,
                        step_size = step_size,
                        max_iter_optimizer = as.integer(max_iter_optimizer),
                        step_shrinkage = step_shrinkage,
                        max_no_update = as.integer(max_no_update)
                    )
                )
            )
        },
        final_loss = .aa_final_loss,
        prepare_output = function(ctx, prep, fit) {
            .aa_euclidean_output(ctx, prep, fit, fit_info = list(method = "pgd"))
        }
    )
}

.aa_fit_pgd_missing <- function(X,
                                M,
                                weight_fun,
                                max_iter,
                                tol,
                                tol_r2,
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
    if (!is.null(weight_fun)) {
        stop("`robust = TRUE` is not supported with `missing = TRUE`.", call. = FALSE)
    }

    if (pseudo_pgd) {
        project <- proj_l1
        grad_S <- grad_S_matrix_l1
        grad_B <- grad_B_matrix_l1
    } else {
        project <- proj_simplex
        grad_S <- grad_S_matrix
        grad_B <- grad_B_matrix
    }

    update_alpha <- delta > 0
    a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
    a_hi <- 1 + delta
    clip <- function(a) pmax(pmin(a, a_hi), a_lo)
    denom_eps <- ifelse(eps > 0, eps, 1e-8)

    A0 <- A
    aB <- B
    a <- rowSums(aB)
    slack_tol <- 1e-6
    if (any(a < a_lo - slack_tol) || any(a > a_hi + slack_tol)) {
        fmt <- "Initialize B marginals are outside the specified delta range [%.3f, %.3f]"
        stop(sprintf(fmt, a_lo, a_hi))
    }
    B <- aB / a

    A <- .aa_pgd_missing_A(aB, X, M, denom_eps)
    R <- .aa_pgd_missing_resid(M, S, A, X) # residuals matrix
    xss <- norm(X, "F")^2 # X sum of squares; invariant for fixed X
    rss <- norm(R, "F")^2 # residual sum of squares -> objective to minimize
    loss_terms <- list(rss = rss, xss = xss, A = A)
    loss[["loss"]][1L] <- rss
    loss[["r2"]][1L] <- 1 - rss / xss

    # Edge case: if max_iter = 0 return initial fit
    if (max_iter == 0) {
        return(list(
            A0 = A0,
            A = A,
            B = aB,
            S = S,
            delta = delta,
            i = 0L,
            loss = loss,
            converged = TRUE
        ))
    }

    converged <- FALSE

    step_S <- step_size
    step_B <- step_size
    step_a <- step_size
    no_update <- 0L

    if (verbose) message("Starting missing-data optimization loop...")
    for (i in seq_len(max_iter)) {
        accepted_update <- FALSE

        # Update S
        grad <- grad_S(S, R, A)
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
        grad_aB <- grad_aB_matrix(aB, X, M, S, A, R, denom_eps)
        grad <- grad_B(grad_aB, a, B)
        for (k in seq_len(max_iter_optimizer)) {
            B_new <- project(B - step_B * grad, eps = eps)
            aB_new <- a * B_new
            A_new <- .aa_pgd_missing_A(aB_new, X, M, denom_eps)
            R_new <- .aa_pgd_missing_resid(M, S, A_new, X)
            rss_new <- norm(R_new, "F")^2
            if (rss_new < rss) {
                B <- B_new
                aB <- aB_new
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
            grad <- grad_alpha(B, grad_aB)
            for (k in seq_len(max_iter_optimizer)) {
                a_new <- clip(a - step_a * grad)
                aB_new <- a_new * B
                A_new <- .aa_pgd_missing_A(aB_new, X, M, denom_eps)
                R_new <- .aa_pgd_missing_resid(M, S, A_new, X)
                rss_new <- norm(R_new, "F")^2
                if (rss_new < rss) {
                    a <- a_new
                    aB <- aB_new
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

        # Check progress
        if (!accepted_update) {
            no_update <- no_update + 1L

            # Copy previous loss entries to keep history.
            loss[["loss"]][i + 1L] <- loss_terms[["rss"]]
            loss[["r2"]][i + 1L] <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

            if (no_update >= max_no_update) {
                fmt <- paste(
                    "Missing-data PGD stalled: no line-search update accepted after",
                    "%d consecutive step shrinkages"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }

            # Refine step sizes for next iteration.
            step_S <- step_S * step_shrinkage
            step_B <- step_B * step_shrinkage
            if (update_alpha) {
                step_a <- step_a * step_shrinkage
            }

            if (verbose) {
                fmt <- paste(
                    "Iteration %d: no line-search update accepted;",
                    "shrinking steps (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
        } else {
            loss_terms <- list(
                rss = rss,
                xss = xss,
                A = A
            )
            loss[["loss"]][i + 1L] <- loss_terms[["rss"]]
            loss[["r2"]][i + 1L] <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

            no_update <- 0L
            converged <- .aa_check_convergence(loss, i, tol, tol_r2, verbose)
            if (converged) break
        }
    }

    list(
        A0 = A0,
        A = A,
        B = aB,
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
                        eps,
                        verbose,
                        loss_fun,
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
    #   rss = ||X - SA||^2 = ||X||^2 - 2*tr(SAXt) + tr(StS AAt) Residual Sum of Squares

    # PGD specific preparations -----------------------------------------------

    # define gradient & projection functions
    if (pseudo_pgd) {
        grad_S <- grad_S_l1
        grad_aB <- grad_aB_l1
        project <- proj_l1
    } else {
        grad_S <- grad_S_trace
        grad_aB <- grad_aB_trace
        project <- proj_simplex
    }

    # Setup alpha updates
    update_alpha <- delta > 0
    a_lo <- ifelse(eps > 0, eps, 1e-8) # to avoid division by 0
    a_lo <- max(1 - delta, a_lo) # lower bound for a
    a_hi <- 1 + delta # upper bound for a
    clip <- function(a) pmax(pmin(a, a_hi), a_lo) # clip a to [a_lo, a_hi]
    a <- rowSums(B)
    slack_tol <- 1e-6
    if (any(a < a_lo - slack_tol) || any(a > a_hi + slack_tol)) {
        fmt <- "Initialize B marginals are outside the specified delta range [%.3f, %.3f]"
        stop(sprintf(fmt, a_lo, a_hi))
    }
    B <- B / a # normalize B to row-stochastic

    # Initialize Variables and loss terms
    A0 <- A # keep initial archetypes for output
    loss_terms <- loss_fun(
        X, A, S,
        weight_fun = weight_fun
    )
    loss[["loss"]][1L] <- loss_terms[["rss"]]
    loss[["r2"]][1L] <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

    # Edge case: if max_iter = 0 return initial fit
    if (max_iter == 0) {
        return(list(
            A0 = A0,
            A = A,
            B = a * B,
            S = S,
            delta = delta,
            i = 0L,
            loss = loss,
            converged = TRUE
        ))
    }

    # Auxiliary variables for efficient loss and gradient updates
    XAt <- loss_terms[["XAt"]] # (N x K)
    AAt <- loss_terms[["AAt"]] # (K x K)
    StS <- loss_terms[["StS"]] # (K x K)
    StX <- loss_terms[["StX"]] # (K x M)
    xss <- loss_terms[["xss"]] # scalar = ||X||^2; invariant for fixed X
    rss <- loss_terms[["rss"]] # scalar = ||X - SA||^2
    row_xss <- loss_terms[["row_xss"]] # (N x 1) row-wise sum of squares; invariant for fixed X
    row_weights <- loss_terms[["row_weights"]]
    StXXt <- loss_terms[["StXXt"]] # (K x N)

    ## Step sizes
    step_S <- step_size # mu_S in the paper
    step_B <- step_size # mu_C in the paper
    step_a <- step_size # mu_a in the paper
    no_update <- 0L
    converged <- FALSE

    # Main optimization loop  -------------------------------------------------

    if (verbose) message("Starting optimization loop...")
    for (i in seq_len(max_iter)) {
        accepted_update <- FALSE

        ## Update S
        grad <- grad_S(S, AAt, XAt, row_weights = row_weights) # (N x K) - diag(a) is absorbed in A
        for (k in seq_len(max_iter_optimizer)) { # line search
            S_new <- project(S - step_S * grad, eps = eps)
            S_new_weighted <- .aa_weight_rows(S_new, row_weights)
            StS_new <- crossprod(S_new_weighted, S_new)
            StX_new <- crossprod(S_new_weighted, X)
            StXXt_new <- tcrossprod(StX_new, X)
            rss_new <- xss - 2 * sum(A * StX_new) + sum(StS_new * AAt)
            if (rss_new < rss) { # update variables
                S <- S_new
                StS <- StS_new
                StX <- StX_new
                rss <- rss_new
                StXXt <- StXXt_new
                step_S <- step_S / step_shrinkage # leave room for shrinkage
                accepted_update <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }

        ## Update B
        grad_aB_step <- grad_aB(B, A, X, StS, StXXt) # (K x N)
        grad <- a * grad_aB_step # grad_B
        for (k in seq_len(max_iter_optimizer)) { # line search
            B_new <- project(B - step_B * grad, eps = eps)
            A_new <- (a * B_new) %*% X
            AAt_new <- tcrossprod(A_new)
            rss_new <- xss - 2 * sum(A_new * StX) + sum(StS * AAt_new)
            if (rss_new < rss) {
                B <- B_new
                A <- A_new
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
            grad <- grad_alpha(B, grad_aB_step) # grad_a
            for (k in seq_len(max_iter_optimizer)) {
                a_new <- clip(a - step_a * grad)
                a_update <- a_new / a # to undo aA multiplication
                A_new <- a_update * A
                AAt_new <- AAt * tcrossprod(a_update) # fast outer product
                rss_new <- xss - 2 * sum(A_new * StX) + sum(StS * AAt_new)
                if (rss_new < rss) { # update variables
                    a <- a_new
                    A <- A_new
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

        # Check progress
        if (!accepted_update) {
            no_update <- no_update + 1L

            # Copy previous loss entries to keep history
            loss[["loss"]][i + 1L] <- loss_terms[["rss"]]
            loss[["r2"]][i + 1L] <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

            if (no_update >= max_no_update) {
                fmt <- paste(
                    "PGD stalled: no line-search update accepted after",
                    "%d consecutive step shrinkages"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }

            # Refine step sizes for next iteration
            step_S <- step_S * step_shrinkage
            step_B <- step_B * step_shrinkage
            step_a <- step_a * step_shrinkage

            if (verbose) {
                fmt <- paste(
                    "Iteration %d: no line-search update accepted;",
                    "shrinking steps (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
        } else { # update accepted
            # Update loss only after accepting an update
            loss_terms <- loss_fun(
                X, A, S,
                weight_fun = weight_fun,
                xss = xss,
                rss = rss,
                row_xss = row_xss,
                row_weights = row_weights,
                StS = StS,
                StX = StX,
                StXXt = StXXt,
                AAt = AAt,
                XAt = XAt
            )
            loss[["loss"]][i + 1L] <- loss_terms[["rss"]]
            loss[["r2"]][i + 1L] <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

            # Refresh cached terms. In the weighted route, robust reweighting
            # changes S_weighted and therefore StS, StX, StXXt, xss, and rss.
            AAt <- loss_terms[["AAt"]]
            XAt <- loss_terms[["XAt"]]
            StS <- loss_terms[["StS"]]
            StX <- loss_terms[["StX"]]
            xss <- loss_terms[["xss"]]
            rss <- loss_terms[["rss"]]
            row_weights <- loss_terms[["row_weights"]]
            StXXt <- loss_terms[["StXXt"]]

            no_update <- 0L
            converged <- .aa_check_convergence(loss, i, tol, tol_r2, verbose)
            if (converged) break
        }
    }

    # Final projection to ensure predict(fit, type = "compositions") == S
    S <- fit_simplex(A, X, eps = eps)

    list(
        A0 = A0,
        A = A,
        B = a * B,
        S = S,
        delta = delta,
        i = i,
        loss = loss,
        row_weights = row_weights,
        converged = converged
    )
}

.aa_pgd_loss <- function(X, A, S, xss = NULL, rss = NULL,
                         StS = NULL, StX = NULL, AAt = NULL, XAt = NULL,
                         StXXt = NULL,
                         ...) {
    xss <- xss %||% norm(X, "F")^2
    AAt <- AAt %||% tcrossprod(A)
    XAt <- XAt %||% tcrossprod(X, A)
    StS <- StS %||% crossprod(S)
    StX <- StX %||% crossprod(S, X)
    rss <- rss %||% (xss - 2 * sum(A * StX) + sum(StS * AAt))
    StXXt <- StXXt %||% tcrossprod(StX, X)

    list(
        rss = rss,
        xss = xss,
        StS = StS,
        StX = StX,
        StXXt = StXXt,
        S_weighted = S,
        AAt = AAt,
        XAt = XAt,
        A = A
    )
}

.aa_pgd_weighted_loss <- function(X, A, S, weight_fun, row_xss = NULL,
                                  AAt = NULL, XAt = NULL, ...) {
    row_xss <- row_xss %||% rowSums(X * X)
    AAt <- AAt %||% tcrossprod(A)
    XAt <- XAt %||% tcrossprod(X, A)

    row_rss <- .aa_trace_row_rss(row_xss, S, XAt, AAt)
    row_weights <- weight_fun(row_rss)
    row_weights <- .aa_check_row_weights(row_weights, nrow(X))

    S_weighted <- .aa_weight_rows(S, row_weights)
    StS <- crossprod(S_weighted, S)
    StX <- crossprod(S_weighted, X)
    StXXt <- tcrossprod(StX, X)

    list(
        rss = sum(.aa_weight_rows(row_rss, row_weights)),
        xss = sum(.aa_weight_rows(row_xss, row_weights)),
        row_xss = row_xss,
        row_rss = row_rss,
        row_weights = row_weights,
        StS = StS,
        StX = StX,
        StXXt = StXXt,
        S_weighted = S_weighted,
        AAt = AAt,
        XAt = XAt,
        A = A
    )
}

# Gradient functions for squared-loss archetypes analysis
# ||X - SBX||^2 = tr(XXt - 2SBXXt + SBXXBtSt)
grad_S_trace <- function(S, AAt, XAt, row_weights = NULL) {
    # 2 * (SBXXtBt - XXtBt) = 2 * (SAAt - XAt)
    .aa_weight_rows(S %*% AAt - XAt, row_weights)
}

# Keep B as argument for consistency with grad_aB_trace
grad_aB_trace <- function(B, A, X, StS, StXXt) {
    # 2 * (StSBXXt - StXXt) = 2 * (StSAXt - StXXt)
    tcrossprod(StS %*% A, X) - StXXt
}

grad_alpha <- function(B, grad_aB) {
    # dR/da = dR/d(aB) * d(aB)/da = dR/d(aB) * B
    rowSums(grad_aB * B)
}

# L1 gradients remove the radial component of the grad
# (the increase of the norm induced by the grad itself)
grad_S_l1 <- function(S, ...) {
    grad <- grad_S_trace(S, ...) # (N x K)
    row_dot <- rowSums(S * grad)
    grad - row_dot
}

grad_aB_l1 <- function(B, ...) {
    grad <- grad_aB_trace(B, ...) # (K x N)
    row_dot <- rowSums(B * grad)
    grad - row_dot
}

# Gradient of functions using the whole residual matrix R = X - SA
# Keep S as argument for consistency with pseudo version
grad_S_matrix <- function(S, R, A) tcrossprod(R, A) # d(R^2)/dS = 2 * R * dR/dS = -2RAt


grad_S_matrix_l1 <- function(S, R, A) {
    grad <- grad_S_matrix(S, R, A)
    grad - rowSums(S * grad)
}

grad_aB_matrix <- function(aB, X, M, S, A, R, eps) {
    # 1) Base case (aB = B, M = 1): with R = SBX - X => dR = S dB X = dR/dB = S (.) X
    #    d(R^2)/dB = d(R^2)/dR * dR/dB = 2 St R Xt
    # 2) For a != 1 we optimize normalized archetypes A(aB) = (aBX) / (aB1 + eps),
    #    so the denominator adds a quotient-rule correction (not just a scaling by a).
    # 3) With missingness, replace aB1 by aBM and R by M o (SA - X),
    #    which masks unobserved entries and yields grad_aB = (XG - M(G o At))t,
    #    where G = Rt S / ((aBM)t + eps).
    grad <- crossprod(R, S) / (t(aB %*% M) + eps)
    t(X %*% grad - M %*% (grad * t(A)))
}

grad_B_matrix <- function(grad_aB, a, B) {
    # d(R^2)/d(B) = d(R^2)/d(aB) * daB/dB = grad_aB * a
    a * grad_aB
}

grad_B_matrix_l1 <- function(grad_aB, a, B) {
    # Remove radial component before simplex projection in the pseudo/L1 variant.
    a * (grad_aB - rowSums(B * grad_aB))
}

.aa_pgd_missing_A <- function(aB, X, M, eps) {
    numerator <- aB %*% X
    denominator <- aB %*% M + eps
    as.matrix(numerator / denominator)
}

.aa_pgd_missing_resid <- function(M, S, A, X) {
    if (!inherits(M, "sparseMatrix")) {
        Xhat <- S %*% A
        Xhat[!M] <- 0
        return(Xhat - X)
    }

    entries <- Matrix::summary(M) # i, j, x for nonzero entries of M
    if (length(entries[["i"]]) == 0L) {
        return(Matrix::sparseMatrix(dims = dim(M), dimnames = dimnames(M)))
    }
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
