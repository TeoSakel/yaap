#' Archetypal Analysis using Non-Negative Least Squares
#'
#' Fits archetypal analysis (AA) by alternating between a non-negative
#' least-squares (NNLS) step for convex-mixtures ($S$, $B$) and an ordinary
#' least-squares (OLS) step for archetype coordinates.
#'
#' @param x data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param init initialization method; see [run_aa()] for available options.
#' @param init_args list of additional arguments for the initialization function.
#' @param weights optional vector of sample weights (default: NULL)
#' @param scale scaling or metric embedding used before fitting; see
#'   [archetypes_pgd()] for details (default: TRUE).
#' @param robust whether to use Tukey bisquare row reweighting (default: FALSE)
#' @param tukey_c tuning constant for Tukey bisquare weights (default: 4.685)
#' @param sd_threshold threshold for feature standard deviation to filter
#'   low-variance features (default: 1e-6).
#' @param max_iter maximum number of iterations (default: 100)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-6)
#' @param tol_r2 convergence tolerance based on R\eqn{^2} (default: 0.9999)
#' @param max_kappa maximum condition number for archetypes (default: 1000)
#' @param eps small positive number to ensure numerical stability
#'   (default: 0 for sparse input 1e-8 for dense)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param bigM large constant to enforce simplex constraint,
#'   or `NULL` to set it automatically.
#' @param max_no_update maximum consecutive iterations without improvement
#'   before considering NNLS stalled (default: 5)
#'
#' @details
#' ## NNLS solver
#'
#' Standard AA alternates between updating the composition matrix S (holding
#' archetypes fixed) and the archetype coordinates A (holding compositions fixed).
#' The NNLS solver reformulates each S-update as a non-negative least-squares
#' problem with a simplex-constraint row appended (the "big-M" trick), making it
#' well-suited for settings where NNLS is more numerically stable than projected
#' gradient descent, such as sparse inputs. The archetype composition matrix B
#' is also updated via NNLS.
#'
#' ## OLS and big-M
#'
#' The A-update solves `min_A ||X - SA||_F^2` as an unconstrained least-squares
#' problem via QR factorisation.
#'
#' `bigM` is the weight of the row appended to enforce the simplex (sum-to-one) constraint.
#' When `NULL` (default) it is set automatically based heuristics. If it fails
#' to enforce the simplex constraint, warning messages will ask you to increase it.
#' if the solver is slow or numerically unstable try decreasing it.
#'
#' @returns An object of class \code{archetypes}.
#'
#' @seealso [run_aa()] for the common entry point and full parameter documentation.
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' archetypes_nnls(as.matrix(toy), K = 3)
#'
#' @references Alcacer, A., Epifanio, I., Mair, S., & Mørup, M. (2025).
#' A Survey on Archetypal Analysis. *arXiv preprint arXiv:2504.12392*.
#' \url{https://arxiv.org/abs/2504.12392}
#'
#' @export
archetypes_nnls <- function(x,
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
                            eps = ifelse(inherits(x, "sparseMatrix"), 0, 1e-8),
                            verbose = FALSE,
                            max_no_update = 5L,
                            # NNLS specific
                            bigM = NULL) {
    .aa_fit_engine(
        call = match.call(),
        x = x,
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
        bigM = bigM,
        max_no_update = max_no_update
    )
}

.aa_nnls_block <- function(ctx,
                           bigM = NULL,
                           max_no_update = 5L) {
    list(
        check = function(ctx) {
            .aa_euclidean_check(ctx)
            if (!is.null(bigM)) {
                stopifnot("`bigM` must be NULL or a positive number" = is_positive(bigM))
            }
            .aa_check_max_no_update(max_no_update)
            invisible(TRUE)
        },
        preprocess = function(ctx) .aa_euclidean_preprocess(ctx, bigM = bigM),
        edge_case = .aa_euclidean_edge_case,
        init = function(ctx, prep) .aa_euclidean_init(ctx, prep, delta = 0),
        fit = function(ctx, prep, init_vars) {
            do.call(
                .aa_fit_nnls,
                c(
                    list(
                        X = prep[["X"]],
                        weight_fun = .aa_weight_fun(ctx[["robust"]], ctx[["tukey_c"]]),
                        max_iter = ctx[["max_iter"]],
                        tol = ctx[["tol"]],
                        tol_r2 = ctx[["tol_r2"]],
                        max_kappa = ctx[["max_kappa"]],
                        eps = ctx[["eps"]],
                        verbose = ctx[["verbose"]],
                        loss_fun = if (ctx[["robust"]])
                            .aa_nnls_weighted_loss_terms
                        else
                            .aa_nnls_loss_terms
                    ),
                    init_vars,
                    list(
                        max_no_update = as.integer(max_no_update)
                    )
                )
            )
        },
        final_loss = .aa_final_loss,
        prepare_output = function(ctx, prep, fit) {
            .aa_euclidean_output(ctx, prep, fit, fit_info = list(method = "nnls"))
        }
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

    xss <- xss %||% norm(X, "F")^2  # computed once
    AAt <- tcrossprod(A)     # (K x K)
    XAt <- tcrossprod(X, A)  # (N x K)
    StS <- crossprod(S)      # (K x K)
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

    row_xss <- row_xss %||% rowSums(X * X)
    AAt <- tcrossprod(A)
    XAt <- tcrossprod(X, A)

    row_rss <- .aa_trace_row_rss(row_xss, S, XAt, AAt)
    row_weights <- weight_fun(row_rss)
    row_weights <- .aa_check_row_weights(row_weights, nrow(X))

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
                         max_no_update) {
    # Nomenclature following arXiv:2504.12392v1:
    #   X ~ SA (N x M) Data Matrix
    #   A = BX (K x M) Archetypes
    #   B = (K x N) Archetypes Coefficients (base transform, C in the paper)
    #   S = (N x K) Archetypes Scores (new coordinates)
    #   rss = ||X - SA||^2 = ||X||^2 - 2*tr(SAXt) + tr(StS AAt) Residual Sum of Squares
    A0 <- A
    Xt <- t(X)  # compute Xt once to reuse for B update
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
    converged <- FALSE
    no_update <- 0L
    max_simplex_error <- 0
    project <- proj_l1
    best_loss_terms <- loss_terms
    best_args <- list(A = A, B = B, S = S)
    last_update_accepted <- TRUE

    # edge case: if max_iter = 0 return initial solution without any updates
    if (max_iter == 0L) {
        return(list(
            # (kappa(A) check removed)
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
        # S update
        S_raw <- fit_nnls(X, t(A), use_svd = FALSE) # Project X to A-simplex
        max_simplex_error <- max(max_simplex_error, max(abs(rowSums(S_raw) - 1)))
        S <- project(S_raw, eps = eps)
        # A update: X = SA
        if (is.null(row_weights)) {
            A <- qr.solve(S, X)
        } else {
            sqrt_weights <- sqrt(row_weights)
            A <- qr.solve(S * sqrt_weights, X * sqrt_weights)
        }
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

        if (loss_terms[["rss"]] < best_loss_terms[["rss"]]) {
            best_loss_terms <- loss_terms
            best_args <- list(A = A, B = B, S = S)
            no_update <- 0L
            last_update_accepted <- TRUE
        } else {
            no_update <- no_update + 1L
            last_update_accepted <- FALSE

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

        # Check convergence
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }
    if (!last_update_accepted) {
        loss <- .aa_update_loss(
            loss,
            i + 1L,
            best_loss_terms,
            verbose = FALSE,
            max_kappa = max_kappa
        )
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
