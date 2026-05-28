#' Archetypal Analysis using Frank-Wolfe Updates
#'
#' Fits an archetypal analysis (AA) model using greedy Frank-Wolfe steps.
#'
#' @param x data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param init initialization method; see [run_aa()] for available options.
#' @param init_args list of additional arguments for the initialization function.
#' @param weights optional vector of sample weights (default: NULL)
#' @param robust FALSE for ordinary squared error, TRUE for "psi.bisquare",
#'   a MASS psi function name, or a custom psi function.
#' @param robust_args list of tuning arguments passed to the robust psi function.
#' @param scale scaling or metric embedding used before fitting; see
#'   [archetypes_pgd()] for details (default: FALSE).
#' @param sd_threshold threshold for feature standard deviation to filter
#'   low-variance features (default: 1e-6)
#' @param max_iter maximum number of outer iterations (default: 100)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-4)
#' @param tol_r2 convergence tolerance based on R\eqn{^2} (default: 0.9999)
#' @param eps small positive number to ensure numerical stability
#'   (default: 0 for sparse input, 1e-8 for dense)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param missing whether to fit the missing-data observed-entry objective.
#'   The default activates missing-data mode when x contains NAs.
#' @param nrep number of random restarts; the best fit is returned (default: 1)
#' @param max_iter_optimizer number of Frank-Wolfe inner steps used for each
#'   alternating update of the composition and coefficient matrices (default: 10)
#' @param fw_step_offset denominator offset in the open-loop Frank-Wolfe step size.
#'   Smaller values give more aggressive vertex moves; must be at least 1 (default: 2).
#' @param max_no_update maximum consecutive outer iterations with no accepted
#'   Frank-Wolfe update before stopping as stalled (default: 5)
#'
#' @details
#' Like [archetypes_pgd()], this method minimises the ordinary squared-error AA
#' objective \eqn{\|X - SBX\|_F^2}. The difference is that each update moves
#' directly toward a simplex vertex selected by the current gradient, rather
#' than taking a gradient step and projecting the result back to the simplex.
#' This projection-free update is the rapid AA strategy proposed by Bauckhage
#' et al. (2015), motivated by the fact that the composition and coefficient
#' columns each lie in a standard simplex.
#'
#' The method is often attractive when a fast, sparse approximation is preferred
#' over the smoother projected-gradient trajectory. The paper also frames AA as
#' an autoencoder because the factorisation maps \eqn{X} through archetypal
#' coefficients and back to a reconstruction \eqn{XBA}.
#'
#' `max_iter_optimizer` controls how many Frank-Wolfe vertex moves are attempted
#' for each alternating composition or coefficient update. Each inner move uses
#' the open-loop step size \eqn{2 / (j + c)}, where \eqn{j} is the inner
#' iteration and \eqn{c} is `fw_step_offset`. Smaller offsets make the first
#' moves more aggressive and closer to a pure vertex replacement; larger offsets
#' preserve more of the current iterate and make the alternating updates more
#' conservative.
#'
#' `method = "fw"` supports ordinary Euclidean fitting with optional scaling,
#' sample weights, robust row weights, random restarts, and sparse matrix
#' input. With `missing = TRUE`, FW minimizes reconstruction error over only
#' observed entries. Robust fitting remains incompatible with missing-data fits.
#'
#' @returns An object of class \code{archetypes}.
#'
#' @seealso [run_aa()] for the common entry point and full parameter documentation.
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' archetypes_fw(as.matrix(toy), K = 3)
#'
#' @references
#' Bauckhage, C., Kersting, K., Hoppe, F., & Thurau, C. (2015).
#' Archetypal analysis as an autoencoder. *Workshop New Challenges in Neural
#' Computation*, 8-15.
#'
#' @export
archetypes_fw <- function(x,
                          K,
                          init = "furthest_sum",
                          init_args = list(),
                          weights = NULL,
                          robust = FALSE,
                          robust_args = list(),
                          scale = FALSE,
                          sd_threshold = 1e-6,
                          max_iter = 100L,
                          tol = 1e-4,
                          tol_r2 = 0.9999,
                          eps = ifelse(inherits(x, "sparseMatrix"), 0, 1e-8),
                          verbose = FALSE,
                          missing = any(is.na(x)),
                          nrep = 1L,
                          max_iter_optimizer = 10L,
                          fw_step_offset = 2,
                          max_no_update = 5L) {
    .aa_fit_engine(
        call = match.call(),
        x = x,
        K = K,
        method = "fw",
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
        nrep = nrep,
        max_iter_optimizer = max_iter_optimizer,
        fw_step_offset = fw_step_offset,
        max_no_update = max_no_update
    )
}

.aa_fw_block <- function(ctx,
                         max_iter_optimizer = 10L,
                         fw_step_offset = 2,
                         max_no_update = 5L) {
    loss_fun <- if (!identical(ctx[["robust"]], FALSE)) .aa_pgd_weighted_loss else .aa_pgd_loss
    fit_fun <- if (ctx[["missing"]]) .aa_fit_fw_missing else .aa_fit_fw

    list(
        check = function(ctx) {
            .aa_euclidean_check(ctx)
            if (!is_count(max_iter_optimizer)) {
                stop("`max_iter_optimizer` must be a positive integer.", call. = FALSE)
            }
            if (!(is_number(fw_step_offset) && fw_step_offset >= 1)) {
                stop("`fw_step_offset` must be a single number greater than or equal to 1.", call. = FALSE)
            }
            .aa_check_max_no_update(max_no_update)
            invisible(TRUE)
        },
        preprocess = function(ctx) .aa_euclidean_preprocess(ctx, bigM = 0),
        edge_case = .aa_euclidean_edge_case,
        init = function(ctx, prep) .aa_euclidean_init(ctx, prep, delta = 0),
        fit = function(ctx, prep, init_vars) {
            common_args <- list(
                X = prep[["X"]],
                max_iter = ctx[["max_iter"]],
                tol = ctx[["tol"]],
                tol_r2 = ctx[["tol_r2"]],
                eps = ctx[["eps"]],
                verbose = ctx[["verbose"]],
                A = init_vars[["A"]],
                B = init_vars[["B"]],
                S = init_vars[["S"]],
                loss = init_vars[["loss"]],
                max_iter_optimizer = as.integer(max_iter_optimizer),
                fw_step_offset = fw_step_offset,
                max_no_update = as.integer(max_no_update)
            )
            if (ctx[["missing"]]) {
                common_args[["M"]] <- prep[["M"]]
            } else {
                common_args[["loss_fun"]] <- loss_fun
                common_args[["weight_fun"]] <- .aa_weight_fun(ctx[["robust"]], ctx[["robust_args"]])
            }
            do.call(fit_fun, common_args)
        },
        final_loss = .aa_final_loss,
        prepare_output = function(ctx, prep, fit) {
            .aa_euclidean_output(ctx, prep, fit, fit_info = list(method = "fw"))
        }
    )
}

.aa_fit_fw_missing <- function(X,
                               M,
                               max_iter,
                               tol,
                               tol_r2,
                               eps,
                               verbose,
                               A,
                               B,
                               S,
                               loss,
                               max_iter_optimizer,
                               fw_step_offset,
                               max_no_update) {
    # Missing FW optimizes only observed entries; B stays row-stochastic.
    denom_eps <- ifelse(eps > 0, eps, 1e-8)
    A0 <- A
    A <- .aa_pgd_missing_A(B, X, M, denom_eps)
    R <- .aa_pgd_missing_resid(M, S, A, X)
    xss <- norm(X, "F")^2
    rss <- norm(R, "F")^2
    loss_terms <- list(rss = rss, xss = xss, A = A)
    loss[["loss"]][1L] <- rss
    loss[["r2"]][1L] <- 1 - rss / xss
    loss[["rloss"]] <- rep(NA_real_, length(loss[["loss"]]))
    loss[["rloss"]][1L] <- rss

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

    best_loss_terms <- loss_terms
    best_args <- list(A = A, B = B, S = S)
    no_update <- 0L
    converged <- FALSE

    if (verbose) message("Starting missing-data Frank-Wolfe optimization loop...")
    for (i in seq_len(max_iter)) {
        i_loss <- i + 1L

        for (j in seq_len(max_iter_optimizer)) {
            grad <- grad_S_matrix(S, R, A)
            S <- .aa_fw_row_step(S, grad, step = 2 / (j + fw_step_offset))
            R <- .aa_pgd_missing_resid(M, S, A, X)
        }

        for (j in seq_len(max_iter_optimizer)) {
            grad <- grad_aB_matrix(B, X, M, S, A, R, denom_eps)
            B <- .aa_fw_row_step(B, grad, step = 2 / (j + fw_step_offset))
            A <- .aa_pgd_missing_A(B, X, M, denom_eps)
            R <- .aa_pgd_missing_resid(M, S, A, X)
        }

        rss <- norm(R, "F")^2
        loss_terms <- list(rss = rss, xss = xss, A = A)

        if (loss_terms[["rss"]] < best_loss_terms[["rss"]]) {
            best_loss_terms <- loss_terms
            best_args <- list(A = A, B = B, S = S)
            no_update <- 0L
        } else {
            no_update <- no_update + 1L
            if (verbose) {
                fmt <- paste(
                    "Iteration %d: missing-data Frank-Wolfe candidate did not improve best loss;",
                    "continuing from current iterate (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                loss[["loss"]][i_loss] <- best_loss_terms[["rss"]]
                loss[["r2"]][i_loss] <- 1 - best_loss_terms[["rss"]] / best_loss_terms[["xss"]]
                loss[["rloss"]][i_loss] <- loss_terms[["rss"]]
                fmt <- paste(
                    "Missing-data Frank-Wolfe stalled: no loss improvement after",
                    "%d consecutive candidate updates"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }
        }

        loss[["loss"]][i_loss] <- best_loss_terms[["rss"]]
        loss[["r2"]][i_loss] <- 1 - best_loss_terms[["rss"]] / best_loss_terms[["xss"]]
        loss[["rloss"]][i_loss] <- loss_terms[["rss"]]

        converged <- .aa_check_convergence(loss, i, tol, tol_r2, verbose)
        if (converged) break
    }

    c(
        best_args,
        list(
            A0 = A0,
            delta = 0,
            i = i,
            loss = loss,
            converged = converged
        )
    )
}

.aa_fit_fw <- function(X,
                       max_iter,
                       tol,
                       tol_r2,
                       eps,
                       verbose,
                       A,
                       B,
                       S,
                       loss,
                       loss_fun,
                       weight_fun,
                       max_iter_optimizer,
                       fw_step_offset,
                       max_no_update) {
    # Paper X is M x N. yaap X is N x M, S is paper A^T, and B is paper B^T.
    A0 <- A
    loss_terms <- loss_fun(X, A, S, weight_fun = weight_fun)
    loss[["loss"]][1L] <- loss_terms[["rss"]]
    loss[["r2"]][1L] <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]
    loss[["rloss"]] <- rep(NA_real_, length(loss[["loss"]]))
    loss[["rloss"]][1L] <- loss_terms[["rss"]]

    if (max_iter == 0L) {
        return(list(
            A0 = A0,
            A = A,
            B = B,
            S = S,
            delta = 0,
            i = 0L,
            loss = loss,
            converged = TRUE,
            row_weights = loss_terms[["row_weights"]]
        ))
    }

    xss <- loss_terms[["xss"]]
    row_xss <- loss_terms[["row_xss"]]
    row_weights <- loss_terms[["row_weights"]]
    AAt <- loss_terms[["AAt"]]
    XAt <- loss_terms[["XAt"]]
    StS <- loss_terms[["StS"]]
    StX <- loss_terms[["StX"]]
    StXXt <- loss_terms[["StXXt"]]
    best_loss_terms <- loss_terms
    best_args <- list(A = A, B = B, S = S)
    no_update <- 0L
    converged <- FALSE

    if (verbose) message("Starting Frank-Wolfe optimization loop...")
    for (i in seq_len(max_iter)) {
        i_loss <- i + 1L

        # Paper Algorithm 1: update compositions while archetypes are fixed.
        for (j in seq_len(max_iter_optimizer)) {
            grad <- grad_S_trace(S, AAt, XAt, row_weights = row_weights)
            S <- .aa_fw_row_step(S, grad, step = 2 / (j + fw_step_offset))
        }
        S_weighted <- .aa_weight_rows(S, row_weights)
        StS <- crossprod(S_weighted, S)
        StX <- crossprod(S_weighted, X)
        StXXt <- tcrossprod(StX, X)

        # Paper Algorithm 2: update archetype coefficients with S fixed.
        for (j in seq_len(max_iter_optimizer)) {
            grad <- grad_aB_trace(B, A, X, StS, StXXt)
            B <- .aa_fw_row_step(B, grad, step = 2 / (j + fw_step_offset))
            A <- B %*% X
        }
        AAt <- tcrossprod(A)
        XAt <- tcrossprod(X, A)

        loss_terms <- loss_fun(
            X, A, S,
            weight_fun = weight_fun,
            xss = xss,
            row_xss = row_xss,
            AAt = AAt,
            XAt = XAt,
            StS = StS,
            StX = StX,
            StXXt = StXXt
        )
        row_weights <- loss_terms[["row_weights"]]

        if (loss_terms[["rss"]] < best_loss_terms[["rss"]]) {
            best_loss_terms <- loss_terms
            best_args <- list(A = A, B = B, S = S)
            no_update <- 0L
        } else {
            no_update <- no_update + 1L
            if (verbose) {
                fmt <- paste(
                    "Iteration %d: Frank-Wolfe candidate did not improve best loss;",
                    "continuing from current iterate (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                loss[["loss"]][i_loss] <- best_loss_terms[["rss"]]
                loss[["r2"]][i_loss] <- 1 - best_loss_terms[["rss"]] / best_loss_terms[["xss"]]
                loss[["rloss"]][i_loss] <- loss_terms[["rss"]]
                fmt <- paste(
                    "Frank-Wolfe stalled: no loss improvement after",
                    "%d consecutive candidate updates"
                )
                warning(sprintf(fmt, max_no_update), call. = FALSE)
                break
            }
        }

        loss[["loss"]][i_loss] <- best_loss_terms[["rss"]]
        loss[["r2"]][i_loss] <- 1 - best_loss_terms[["rss"]] / best_loss_terms[["xss"]]
        loss[["rloss"]][i_loss] <- loss_terms[["rss"]]

        converged <- .aa_check_convergence(loss, i, tol, tol_r2, verbose)
        if (converged) break
    }

    best_args[["S"]] <- fit_simplex(best_args[["A"]], X, eps = eps)
    c(
        best_args,
        list(
            A0 = A0,
            delta = 0,
            i = i,
            loss = loss,
            converged = converged,
            row_weights = best_loss_terms[["row_weights"]]
        )
    )
}

.aa_fw_row_step <- function(W, gradient, step) {
    # Move each simplex row toward the vertex with the smallest gradient entry.
    gradient <- as.matrix(gradient)
    vertex <- max.col(-gradient, ties.method = "first")  # for reproducibility, ties should be rare
    vertex <- cbind(seq_len(nrow(W)), vertex)
    # W = (1 - step) * W + step * e_vertex, where e_vertex = onehot encoding of vertex.
    W <- (1 - step) * W
    W[vertex] <- W[vertex] + step
    W
}
