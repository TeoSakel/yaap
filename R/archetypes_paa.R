#' Probabilistic Archetypal Analysis using Projected Gradient Descent
#'
#' Fits probabilistic archetypal analysis (PAA) by minimising the
#' log-likelihood of the data of several exponential-family model via
#' projected gradient descent (PGD).
#'
#' @param x data matrix (rows = samples, columns = dimensions).
#' @param K number of archetypes.
#' @param family observation family. One of `"gaussian"`, `"binomial"`,
#'   `"poisson"`, or `"multinomial"` (default: `"gaussian"`).
#' @param init initialization method; see [run_aa()] for available options.
#' @param init_args list of additional arguments for the initialization function.
#' @param max_iter maximum number of outer iterations (default: 100)
#' @param tol convergence tolerance based on objective loss (default: 1e-4)
#' @param tol_r2 convergence tolerance based on R\eqn{^2} (default: 0.9999)
#' @param eps small positive number for numerical stability (default: 1e-8)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param step_size initial line-search step size (default: 1.0)
#' @param max_iter_optimizer maximum line-search iterations per update (default: 10)
#' @param step_shrinkage factor used to shrink rejected line-search steps (default: 0.5)
#' @param max_no_update maximum consecutive outer iterations with no accepted
#'   line-search update before stopping as stalled (default: 5)
#'
#' @details
#' ## Observation families
#'
#' PAA models each observation as drawn from an exponential-family distribution
#' whose natural parameter is a convex combination of K archetypal profiles.
#' Before optimisation begins, each data point is mapped to a fixed profile by
#' computing its per-sample MLE parameter under the chosen family (e.g. the raw
#' values for `"gaussian"`, high probability for `"binomial"`, etc.). The
#' archetypes are then found as convex combinations of these fixed profiles.
#' Built-in families and their data requirements:
#'
#' \describe{
#'   \item{`gaussian`}{Squared reconstruction error; data can be any real matrix.}
#'   \item{`binomial`}{Binary cross-entropy; each entry must lie in \[0, 1\].}
#'   \item{`poisson`}{Poisson log-likelihood; all entries must be non-negative.}
#'   \item{`multinomial`}{Multinomial log-likelihood; entries must be non-negative
#'                        with positive row sums, treated as count vectors.}
#' }
#'
#' Note: `gaussian` PAA is equivalent to standard AA; use [archetypes_pgd()]
#' instead for a more efficient solver with support for missing data, sample
#' weights, and metric scaling.
#'
#' @returns An object of class \code{archetypes}.
#'
#' @seealso [run_aa()] for the common entry point and full parameter documentation.
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' archetypes_paa(as.matrix(toy), K = 3)
#'
#' @references
#' Seth, S., & Eugster, M. J. A. (2014). Probabilistic archetypal analysis.
#' arXiv:1312.7604.
#'
#' @export
archetypes_paa <- function(x,
                           K,
                           family = c("gaussian", "binomial", "poisson", "multinomial"),
                           init = "furthest_sum",
                           init_args = list(),
                           max_iter = 100L,
                           tol = 1e-4,
                           tol_r2 = 0.9999,
                           eps = 1e-8,
                           verbose = FALSE,
                           step_size = 1.0,
                           max_iter_optimizer = 10L,
                           step_shrinkage = 0.5,
                           max_no_update = 5L) {
    family <- .aa_paa_normalize_family(family)
    .aa_fit_engine(
        call = match.call(),
        x = x,
        K = K,
        method = "paa",
        family = family,
        init = init,
        init_args = init_args,
        weights = NULL,
        scale = FALSE,
        robust = FALSE,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        eps = eps,
        verbose = verbose,
        missing = FALSE,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = max_no_update
    )
}

.aa_paa_normalize_family <- function(family) {
    match.arg(family, c("gaussian", "binomial", "poisson", "multinomial"))
}

.aa_paa_block <- function(ctx,
                          step_size = 1.0,
                          max_iter_optimizer = 10L,
                          step_shrinkage = 0.5,
                          max_no_update = 5L) {
    family <- .aa_paa_normalize_family(ctx[["family"]])
    list(
        check = function(ctx) {
            if (ctx[["missing"]]) {
                stop("`missing = TRUE` is not supported for `method = 'paa'`.", call. = FALSE)
            }
            if (!identical(ctx[["robust"]], FALSE)) {
                stop("`robust` is not supported for `method = 'paa'`.", call. = FALSE)
            }
            if (!is.null(ctx[["weights"]])) {
                stop("`weights` are not supported for `method = 'paa'`.", call. = FALSE)
            }
            if (!identical(ctx[["scale"]], FALSE)) {
                warning("`scale` is ignored for `method = 'paa'`.", call. = FALSE)
            }
            .aa_check_fit_controls(ctx)
            if (!(is_number(ctx[["eps"]]) && ctx[["eps"]] > 0)) {
                stop("`eps` must be positive for `method = \"paa\"`.", call. = FALSE)
            }
            .aa_check_projected_gradient_controls(
                step_size = step_size,
                max_iter_optimizer = max_iter_optimizer,
                step_shrinkage = step_shrinkage,
                max_no_update = max_no_update
            )
        },
        preprocess = function(ctx) {
            spec <- .aa_paa_family(family, eps = ctx[["eps"]])
            prep <- spec[["prepare"]](ctx[["x"]])
            prep[["data_X"]] <- prep[["X"]]
            prep[["X"]] <- prep[["P"]]
            prep[["spec"]] <- spec
            prep[["family"]] <- family
            prep
        },
        edge_case = function(ctx, prep) NULL,
        init = function(ctx, prep) .aa_euclidean_init(ctx, prep, delta = 0),
        fit = function(ctx, prep, init_vars) {
            .aa_fit_paa_pgd(
                X = prep[["data_X"]],
                P = prep[["P"]],
                spec = prep[["spec"]],
                max_iter = ctx[["max_iter"]],
                tol = ctx[["tol"]],
                tol_r2 = ctx[["tol_r2"]],
                eps = ctx[["eps"]],
                verbose = ctx[["verbose"]],
                A0 = init_vars[["A"]],
                B = init_vars[["B"]],
                S = init_vars[["S"]],
                loss = init_vars[["loss"]],
                step_size = step_size,
                max_iter_optimizer = max_iter_optimizer,
                step_shrinkage = step_shrinkage,
                max_no_update = as.integer(max_no_update)
            )
        },
        final_loss = .aa_final_loss,
        prepare_output = function(ctx, prep, fit) {
            .aa_euclidean_output(
                ctx,
                prep,
                fit,
                fit_info = list(method = "paa", scaling = "none")
            )
        }
    )
}

.aa_paa_family <- function(family, eps) {
    family <- .aa_paa_normalize_family(family)
    clip_prob <- function(Y) pmin(pmax(Y, eps), 1 - eps)
    clip_pos <- function(Y) pmax(Y, eps)
    as_numeric_matrix <- function(data) {
        as.matrix(data)
    }
    normalize_rows <- function(X) {
        X <- pmax(X, eps)
        X / rowSums(X)
    }

    switch(family,
        gaussian = list(
            family = family,
            prepare = function(data) {
                X <- as_numeric_matrix(data)
                list(X = X, P = X)
            },
            objective = function(X, Y) norm(X - Y, "F")^2,
            gradient = function(X, Y) 2 * (Y - X)
        ),
        binomial = list(
            family = family,
            prepare = function(data) {
                X <- as_numeric_matrix(data)
                if (!all(X >= 0 & X <= 1)) {
                    stop("Binomial data must lie in [0, 1].", call. = FALSE)
                }
                P <- clip_prob(X)
                list(X = X, P = P)
            },
            objective = function(X, Y) {
                Y <- clip_prob(Y)
                -sum(X * log(Y) + (1 - X) * log(1 - Y))
            },
            gradient = function(X, Y) {
                Y <- clip_prob(Y)
                (1 - X) / (1 - Y) - X / Y
            }
        ),
        poisson = list(
            family = family,
            prepare = function(data) {
                X <- as_numeric_matrix(data)
                if (!all(X >= 0)) {
                    stop("Poisson data must be non-negative.", call. = FALSE)
                }
                list(X = X, P = clip_pos(X))
            },
            objective = function(X, Y) {
                Y <- clip_pos(Y)
                sum(-X * log(Y) + Y)
            },
            gradient = function(X, Y) {
                Y <- clip_pos(Y)
                1 - X / Y
            }
        ),
        multinomial = list(
            family = family,
            prepare = function(data) {
                X <- as_numeric_matrix(data)
                if (!all(X >= 0)) {
                    stop("Multinomial data must be non-negative.", call. = FALSE)
                }
                if (!all(rowSums(X) > 0)) {
                    stop("Multinomial rows must have positive totals.", call. = FALSE)
                }
                list(X = X, P = normalize_rows(X))
            },
            objective = function(X, Y) {
                Y <- normalize_rows(Y)
                -sum(X * log(Y))
            },
            gradient = function(X, Y) {
                Y <- normalize_rows(Y)
                -X / Y
            }
        )
    )
}

.aa_fit_paa_pgd <- function(X,
                            P,
                            spec,
                            max_iter,
                            tol,
                            tol_r2,
                            eps,
                            verbose,
                            A0,
                            B,
                            S,
                            loss,
                            step_size,
                            max_iter_optimizer,
                            step_shrinkage,
                            max_no_update) {
    # Nomenclature follows Seth and Eugster (2014), transposed to yaap's
    # row-oriented convention:
    #   X: observed data, rows are observations x_n.
    #   P: paper Θ^T, the fixed MLE parameter profiles derived from X.
    #   B: paper W^T, archetype generator weights over observed profiles.
    #   S: paper H^T, sample compositions over archetypes.
    #   A = B %*% P: paper Z^T, archetypal parameter profiles.
    #   Y = S %*% A: fitted observation parameters Θ W H.
    project <- proj_l1
    pseudo_grad <- function(Z, grad) grad - rowSums(Z * grad)

    A <- B %*% P
    Y <- S %*% A
    obj <- spec[["objective"]](X, Y)
    loss_ref <- max(obj, .Machine$double.eps)
    loss[["loss"]][1L] <- obj
    loss[["r2"]][1L] <- 1 - obj / loss_ref
    converged <- FALSE
    no_update <- 0L
    step_S <- step_size
    step_B <- step_size

    if (max_iter == 0L) {
        return(list(
            A0 = A0, A = A, B = B, S = S, i = 0L,
            loss = loss, converged = TRUE
        ))
    }

    if (verbose) message("Starting PAA optimization loop...")
    for (i in seq_len(max_iter)) {
        j <- i + 1L
        accepted_update <- FALSE

        # S-update: hold archetypal profiles A fixed and fit each observation's
        # simplex composition under the selected exponential-family likelihood.
        grad_Y <- spec[["gradient"]](X, Y)
        grad_S <- pseudo_grad(S, tcrossprod(grad_Y, A))
        for (k in seq_len(max_iter_optimizer)) {
            S_new <- project(S - step_S * grad_S, eps = eps)
            Y_new <- S_new %*% A
            obj_new <- spec[["objective"]](X, Y_new)
            if (obj_new < obj) {
                S <- S_new
                Y <- Y_new
                obj <- obj_new
                step_S <- step_S / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }

        # B-update: hold compositions S fixed and move archetypes within the
        # convex hull of the fixed parameter profiles P.
        if (accepted_update) {
            grad_Y <- spec[["gradient"]](X, Y)
        }
        grad_A <- crossprod(S, grad_Y)
        grad_B <- pseudo_grad(B, tcrossprod(grad_A, P))
        for (k in seq_len(max_iter_optimizer)) {
            B_new <- project(B - step_B * grad_B, eps = eps)
            A_new <- B_new %*% P
            Y_new <- S %*% A_new
            obj_new <- spec[["objective"]](X, Y_new)
            if (obj_new < obj) {
                B <- B_new
                A <- A_new
                Y <- Y_new
                obj <- obj_new
                step_B <- step_B / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_B <- step_B * step_shrinkage
        }

        # TODO: consider adding slack parametrers for poisson and potentially binomial
        if (!accepted_update) {
            no_update <- no_update + 1L
            for (nm in names(loss)) {
                loss[[nm]][j] <- loss[[nm]][i]
            }
            if (verbose) {
                fmt <- paste(
                    "Iteration %d: PAA candidate did not improve loss;",
                    "keeping previous iterate (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                warning(sprintf(
                    "PAA PGD stalled: no loss improvement after %d consecutive updates",
                    max_no_update
                ), call. = FALSE)
                break
            }
            next
        }

        no_update <- 0L
        loss[["loss"]][j] <- obj
        loss[["r2"]][j] <- 1 - obj / loss_ref
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, verbose)
        if (converged) break
    }

    list(A0 = A0, A = A, B = B, S = S, i = i, loss = loss, converged = converged)
}

.aa_paa_predict_S <- function(object,
                              newdata,
                              max_iter = 100L,
                              tol = 1e-4,
                              eps = 1e-8,
                              step_size = 1.0,
                              max_iter_optimizer = 10L,
                              step_shrinkage = 0.5,
                              ...) {
    family <- object[["family"]] %||% "gaussian"
    spec <- .aa_paa_family(family, eps = eps)
    prep <- spec[["prepare"]](newdata)
    X <- prep[["X"]]
    P <- prep[["P"]]
    A <- .aa_input_coordinates_matrix(object)
    if (ncol(P) != ncol(A)) {
        fmt <- "`newdata` has %d columns but `object$coordinates` has %d columns"
        stop(sprintf(fmt, ncol(P), ncol(A)), call. = FALSE)
    }
    S <- fit_simplex(A, P, eps = eps)
    Y <- S %*% A
    obj <- spec[["objective"]](X, Y)
    step_S <- step_size
    converged <- FALSE
    eps2 <- 1e-16 #  small constant to prevent division by zero in convergence check
    for (i in seq_len(max_iter)) {
        grad_Y <- spec[["gradient"]](X, Y)
        grad_S <- tcrossprod(grad_Y, A)
        grad_S <- grad_S - rowSums(grad_S * S)
        accepted <- FALSE
        for (k in seq_len(max_iter_optimizer)) {
            S_new <- proj_l1(S - step_S * grad_S, eps = eps)
            Y_new <- S_new %*% A
            obj_new <- spec[["objective"]](X, Y_new)
            if (obj_new < obj) {
                S <- S_new
                Y <- Y_new
                obj <- obj_new
                step_S <- step_S / step_shrinkage
                accepted <- TRUE
                converged <- (obj - obj_new) < tol * max(obj, eps2)
                break
            }
            step_S <- step_S * step_shrinkage
        }
        if (!accepted || converged) break
    }
    colnames(S) <- rownames(A)
    rownames(S) <- rownames(X)
    S
}
