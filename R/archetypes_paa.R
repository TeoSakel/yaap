#' Probabilistic Archetypal Analysis using Projected Gradient Descent
#'
#' Fits probabilistic archetypal analysis (PAA) with one projected-gradient
#' optimizer shared across likelihood families. Archetypes are convex
#' combinations of family-specific parameter profiles derived from the data.
#'
#' @param data Data matrix with rows as samples.
#' @param K Number of archetypes.
#' @param family Observation family. One of `"gaussian"`, `"bernoulli"`,
#'   `"poisson"`, or `"multinomial"`.
#' @param init Initialization method, function, or coordinate matrix.
#' @param init_args List of additional arguments for the initialization method.
#' @param max_iter Maximum number of outer iterations.
#' @param tol Convergence tolerance based on objective loss.
#' @param tol_r2 Convergence tolerance based on relative objective improvement.
#' @param max_kappa Maximum condition number warning threshold.
#' @param eps Small positive number for numerical stability.
#' @param verbose Whether to print progress messages.
#' @param step_size Initial line-search step size.
#' @param max_iter_optimizer Maximum line-search iterations per update.
#' @param step_shrinkage Factor used to shrink rejected line-search steps.
#' @param max_no_update Maximum consecutive outer iterations with no accepted
#'   line-search update before stopping as stalled.
#'
#' @returns An object of class \code{\link{archetypes}}.
#'
#' @references
#' Seth, S., & Eugster, M. J. A. (2014). Probabilistic archetypal analysis.
#' arXiv:1312.7604.
#'
#' @export
archetypes_paa <- function(data,
                           K,
                           family = c("gaussian", "bernoulli", "poisson", "multinomial"),
                           init = "furthest_sum",
                           init_args = list(),
                           max_iter = 100L,
                           tol = 1e-6,
                           tol_r2 = 0.9999,
                           max_kappa = 1000,
                           eps = 1e-8,
                           verbose = FALSE,
                           step_size = 1.0,
                           max_iter_optimizer = 10L,
                           step_shrinkage = 0.5,
                           max_no_update = 5L) {
    call <- match.call()
    family <- match.arg(family)
    .aa_check_paa_inputs(
        data = data,
        K = K,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = max_no_update
    )

    spec <- .aa_paa_family(family, eps = eps)
    prep <- spec[["prepare"]](data)
    X <- prep[["X"]]
    P <- prep[["P"]]

    init_vars <- .aa_init_vars(
        X = P,
        K = K,
        init = init,
        init_args = init_args,
        eps = eps,
        max_iter = max_iter,
        verbose = verbose,
        delta = 0
    )

    fit <- .aa_fit_paa_pgd(
        X = X,
        P = P,
        spec = spec,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        A0 = init_vars[["A"]],
        B = init_vars[["B"]],
        S = init_vars[["S"]],
        loss = init_vars[["loss"]],
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = as.integer(max_no_update)
    )

    archetype_names <- rownames(fit[["B"]])
    if (is.null(archetype_names))
        archetype_names <- paste0("A", seq_len(K))
    rownames(fit[["A"]]) <- rownames(fit[["A0"]]) <- rownames(fit[["B"]]) <-
        colnames(fit[["S"]]) <- archetype_names
    colnames(fit[["B"]]) <- rownames(fit[["S"]]) <- rownames(P)
    colnames(fit[["A"]]) <- colnames(fit[["A0"]]) <- colnames(P)

    j <- fit[["i"]] + 1L
    loss <- as.data.frame(fit[["loss"]])[seq_len(j), , drop = FALSE]
    rownames(loss) <- NULL

    if (!fit[["converged"]])
        warning(sprintf("Algorithm did not converge after %d iterations", max_iter), call. = FALSE)
    if (verbose) {
        fmt <- ifelse(fit[["converged"]], "Converged after %d iterations:",
                      "Final iteration %d:")
        fmt <- paste(fmt, "loss = %.4g, R2 = %.3f")
        message(sprintf(fmt, fit[["i"]], loss[j, "loss"], loss[j, "r2"]))
    }

    archetypes(
        call = call,
        data = as.matrix(data),
        init = fit[["A0"]],
        coordinates = fit[["A"]],
        coefficients = fit[["B"]],
        compositions = fit[["S"]],
        loss = loss,
        converged = fit[["converged"]],
        family = family
    )
}

.aa_check_paa_inputs <- function(data,
                                 K,
                                 max_iter,
                                 tol,
                                 tol_r2,
                                 max_kappa,
                                 eps,
                                 step_size,
                                 max_iter_optimizer,
                                 step_shrinkage,
                                 max_no_update) {
    stopifnot("data must be a matrix-like object" =
                  is.matrix(data) || inherits(data, "data.frame") ||
                  inherits(data, "sparseMatrix"))
    stopifnot("K must be an integer" = K == as.integer(K))
    stopifnot("K must be an integer greater or equal to 1" = K >= 1L)
    stopifnot("K cannot be greater than number of samples" = K <= nrow(data))
    stopifnot("max_iter must be a non-negative integer" =
                  max_iter == as.integer(max_iter) && max_iter >= 0L)
    stopifnot("tol must be positive" = tol > 0)
    stopifnot("tol_r2 must be between (0, 1)" = tol_r2 >= 0 && tol_r2 <= 1)
    stopifnot("max_kappa must be >=1" = max_kappa >= 1)
    stopifnot("eps must be positive" = eps > 0)
    stopifnot("step_size must be positive" = step_size > 0)
    stopifnot("max_iter_optimizer must be a positive integer" =
                  max_iter_optimizer == as.integer(max_iter_optimizer) &&
                  max_iter_optimizer >= 1L)
    stopifnot("step_shrinkage must be between (0, 1)" =
                  step_shrinkage > 0 && step_shrinkage < 1)
    stopifnot("max_no_update must be a positive integer" =
                  max_no_update == as.integer(max_no_update) && max_no_update >= 1L)
    invisible(TRUE)
}

.aa_paa_family <- function(family, eps) {
    clip_prob <- function(Y) pmin(pmax(Y, eps), 1 - eps)
    clip_pos <- function(Y) pmax(Y, eps)
    as_numeric_matrix <- function(data) {
        X <- as.matrix(data)
        stopifnot("data must be numeric" = is.numeric(X))
        stopifnot("data contains missing values" = !any(is.na(X)))
        stopifnot("data contains non-finite values" = all(is.finite(X)))
        X
    }
    normalize_rows <- function(X) {
        X <- pmax(X, eps)
        X / rowSums(X)
    }

    switch(
        family,
        gaussian = list(
            family = family,
            prepare = function(data) {
                X <- as_numeric_matrix(data)
                list(X = X, P = X)
            },
            objective = function(X, Y) norm(X - Y, "F")^2,
            gradient = function(X, Y) 2 * (Y - X)
        ),
        bernoulli = list(
            family = family,
            prepare = function(data) {
                X <- as_numeric_matrix(data)
                stopifnot("Bernoulli data must lie in [0, 1]" =
                              all(X >= 0 & X <= 1))
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
                stopifnot("Poisson data must be non-negative" = all(X >= 0))
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
                stopifnot("Multinomial data must be non-negative" = all(X >= 0))
                stopifnot("Multinomial rows must have positive totals" = all(rowSums(X) > 0))
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
                            max_kappa,
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
    # Nomenclature follows Seth and Eugster (2014), transposed to YAAAP's
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
    loss_terms <- list(rss = obj, xss = loss_ref, StS = crossprod(S), A = A)
    loss <- .aa_update_loss(loss, 1L, loss_terms, verbose = verbose, max_kappa = max_kappa)
    converged <- FALSE
    no_update <- 0L
    step_S <- step_size
    step_B <- step_size

    if (max_iter == 0L) {
        return(list(A0 = A0, A = A, B = B, S = S, i = 0L,
                    loss = loss, converged = TRUE))
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
        if (accepted_update)
            grad_Y <- spec[["gradient"]](X, Y)
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

        # TODO: consider adding slack parametrers for poisson and potentially bernoulli
        if (!accepted_update) {
            no_update <- no_update + 1L
            for (nm in names(loss))
                loss[[nm]][j] <- loss[[nm]][i]
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
        loss_terms <- list(
            rss = obj,
            xss = loss_ref,
            StS = if (i %% 10L == 0L && max_kappa > 1) crossprod(S) else NULL,
            A = A
        )
        loss <- .aa_update_loss(loss, j, loss_terms, verbose = verbose, max_kappa = max_kappa)
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }

    list(A0 = A0, A = A, B = B, S = S, i = i, loss = loss, converged = converged)
}

.aa_paa_predict_S <- function(object,
                              newdata,
                              max_iter = 100L,
                              tol = 1e-6,
                              eps = 1e-8,
                              step_size = 1.0,
                              max_iter_optimizer = 10L,
                              step_shrinkage = 0.5,
                              ...) {
    family <- object[["family"]]
    if (is.null(family))
        family <- "gaussian"
    spec <- .aa_paa_family(family, eps = eps)
    prep <- spec[["prepare"]](newdata)
    X <- prep[["X"]]
    P <- prep[["P"]]
    A <- object[["coordinates"]]
    if (ncol(P) != ncol(A)) {
        fmt <- "`newdata` has %d columns but `object$coordinates` has %d columns"
        stop(sprintf(fmt, ncol(P), ncol(A)), call. = FALSE)
    }
    S <- fit_simplex(A, P, eps = eps)
    Y <- S %*% A
    obj <- spec[["objective"]](X, Y)
    step_S <- step_size
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
                if (abs(obj - obj_new) < tol * max(obj, .Machine$double.eps))
                    accepted <- FALSE
                obj <- obj_new
                step_S <- step_S / step_shrinkage
                accepted <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }
        if (!accepted) break
    }
    colnames(S) <- rownames(A)
    rownames(S) <- rownames(X)
    S
}
