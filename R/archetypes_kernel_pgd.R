#' Kernel Archetypal Analysis using Projected Gradient Descent
#'
#' Fits archetypal analysis (AA) in a reproducing kernel Hilbert space (RKHS)
#' using projected gradient descent (PGD) to minimise the kernel reconstruction
#' error \eqn{\|K - SA\|_F^2}, where \eqn{K} is the kernel Gram matrix, \eqn{S}
#' is the composition matrix, and \eqn{A} is the archetype representation in the
#' RKHS.
#'
#' @param x data matrix (rows = samples) or an `N x N` Gram matrix when
#'   `kernel = "precomputed"`.
#' @param K number of archetypes.
#' @param kernel kernel specification. One of `"linear"`, `"rbf"`,
#'   `"laplace"`, `"polynomial"`, `"precomputed"`, or a function returning an
#'   `N x N` Gram matrix.
#' @param kernel_args list of arguments passed to the kernel.
#' @param data optional original data matrix attached to precomputed-kernel fits
#'   so coordinates() can return an input-space proxy.
#' @param init initialization method; see [run_aa()] for available options.
#' @param init_args list of additional arguments for the initialization function.
#' @param robust robust row reweighting selector. Use `FALSE` for ordinary
#'   squared error, `TRUE` for `"psi.bisquare"`, a MASS psi function name,
#'   or a custom psi function. See [MASS::rlm()] for psi details; `method =
#'   "MM"` is not supported because it is not applicable to AA.
#' @param robust_args list of tuning arguments passed to the robust psi function.
#' @param max_iter maximum number of outer iterations (default: 100)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-6)
#' @param tol_r2 convergence tolerance based on R\eqn{^2} (default: 0.9999)
#' @param eps small positive number to ensure numerical stability (default: 1e-8)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param delta maximum allowed relaxation of archetype convexity constraint (default: 0)
#' @param pseudo_pgd whether to use pseudo projected gradient descent (default: TRUE)
#' @param step_size initial line-search step size (default: 1.0)
#' @param max_iter_optimizer maximum line-search iterations per update (default: 10)
#' @param step_shrinkage factor used to shrink rejected line-search steps (default: 0.5)
#' @param max_no_update maximum consecutive outer iterations with no accepted
#'   line-search update before stopping as stalled (default: 5)
#'
#' @details
#' ## Kernels and the Gram matrix
#'
#' A kernel function \eqn{k(x_i, x_j)} measures pairwise similarity between
#' samples in some transformed feature space. The \eqn{N \times N} Gram matrix
#' with entries \eqn{G_{ij} = k(x_i, x_j)} is the only input the solver needs.
#' Built-in kernels and their `kernel_args`:
#'
#' \describe{
#'   \item{`linear`}{`tcrossprod(X)`.}
#'   \item{`rbf`}{Gaussian RBF: `exp(-0.5 * dist(X)^2 / sigma^2)`.}
#'   \item{`laplace`}{Laplace: `exp(-dist(X, method = "manhattan") / sigma)`.}
#'   \item{`polynomial`}{`(gamma * tcrossprod(X) + coef0)^degree`.}
#' }
#'
#' Kernel parameters passed via `kernel_args`:
#' - `sigma`: bandwidth for `rbf` and `laplace`; auto-selected when `NULL`.
#' - `gamma`: scale for `polynomial` (default `1/ncol(x)`).
#' - `degree`: polynomial degree (default 3).
#' - `coef0`: polynomial intercept (default 1).
#'
#' Note: a "linear" kernel is equivalent to standard pgd AA on the original
#' data but less efficient.
#'
#' To use a precomputed Gram matrix, compute it externally, pass it as `x`,
#' and set `kernel = "precomputed"`. Supply the original data as the `data`
#' argument if you want coordinates() to return an input-space proxy.
#'
#' Kernel initializers that can be expressed through row selections or
#' coefficient matrices are supported: `"random"`, `"furthest_first"`,
#' `"kmeans_pp"`, `"furthest_sum"`, and `"dirichlet"`. Candidate batching can
#' be requested through `init_args = list(batch_size = ..., batch_type = ...)`;
#' distal batching uses distances implied by the Gram matrix.
#'
#' ## Archetype coordinates
#'
#' Because the space where similarities are measured is defined by the kernel,
#' archetype coordinates are not directly available as ordinary input-space
#' coordinates. The coordinates() accessor returns a proxy when the original
#' input data are stored: A = B X, using the fitted generator weights B.
#' For linear kernels this equals the usual Euclidean representation. For
#' nonlinear kernels it is only a plotting and inspection proxy; distances in
#' this coordinate matrix do not carry kernel-metric meaning.
#'
#' When the original data are unavailable, as with some precomputed Gram
#' matrices, coordinates() returns NULL. For an alternative feature-space view,
#' plot.kernel_archetypes() with projection = "pca" projects archetype
#' compositions onto kernel PCA components of the data.
#'
#' @returns An object of class \code{kernel_archetypes}, extending
#'   \code{archetypes}. Kernel methods expose Hilbert-space residual norms,
#'   input-space coordinate proxies when data are stored, and explicitly reject
#'   out-of-sample \code{predict()} until cross-Gram evaluation is part of the API.
#'
#' @seealso [run_aa()] for the common entry point and full parameter documentation.
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' archetypes_kernel_pgd(as.matrix(toy), K = 3)
#'
#' @references
#' Mørup, M., & Hansen, L. K. (2012).
#' Archetypal analysis for machine learning and data mining.
#' *Neurocomputing*, 80, 54-63. \doi{10.1016/j.neucom.2011.06.033}
#'
#' @export
archetypes_kernel_pgd <- function(x,
                                  K,
                                  kernel = c(
                                      "rbf", "laplace", "polynomial",
                                      "linear", "precomputed"
                                  ),
                                  kernel_args = list(),
                                  data = NULL,
                                  init = "furthest_sum",
                                  init_args = list(),
                                  robust = FALSE,
                                  robust_args = list(),
                                  max_iter = 100L,
                                  tol = 1e-6,
                                  tol_r2 = 0.9999,
                                  eps = 1e-8,
                                  verbose = FALSE,
                                  delta = 0,
                                  pseudo_pgd = TRUE,
                                  step_size = 1.0,
                                  max_iter_optimizer = 10L,
                                  step_shrinkage = 0.5,
                                  max_no_update = 5L) {
    if (!is.function(kernel)) {
        kernel <- match.arg(kernel)
    }
    .aa_fit_engine(
        call = match.call(),
        x = x,
        K = K,
        method = "kernel",
        init = init,
        init_args = init_args,
        weights = NULL,
        scale = FALSE,
        robust = robust,
        robust_args = robust_args,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        eps = eps,
        verbose = verbose,
        missing = FALSE,
        data = data,
        kernel = kernel,
        kernel_args = kernel_args,
        delta = delta,
        pseudo_pgd = pseudo_pgd,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = max_no_update
    )
}

.aa_kernel_block <- function(ctx,
                             kernel = c("linear", "rbf", "laplace", "polynomial", "precomputed"),
                             kernel_args = list(),
                             delta = 0,
                             pseudo_pgd = TRUE,
                             step_size = 1.0,
                             max_iter_optimizer = 10L,
                             step_shrinkage = 0.5,
                             max_no_update = 5L) {
    if (!is.function(kernel)) {
        kernel <- match.arg(kernel)
    }

    list(
        check = function(ctx) {
            if (ctx[["missing"]]) {
                stop("`missing = TRUE` is only supported for `method = 'pgd'`.", call. = FALSE)
            }
            if (!is.null(ctx[["weights"]])) {
                stop("`weights` are not supported for `method = 'kernel'`.", call. = FALSE)
            }
            if (!identical(ctx[["scale"]], FALSE)) {
                warning("`scale` is ignored for `method = 'kernel'`.", call. = FALSE)
            }
            if (!is.list(kernel_args)) {
                stop("`kernel_args` must be a list.", call. = FALSE)
            }
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
        preprocess = function(ctx) {
            prep <- .aa_kernel_prepare(
                x = ctx[["x"]],
                kernel = kernel,
                kernel_args = kernel_args,
                data = ctx[["data"]]
            )
            .aa_check_kernel_inputs(
                ctx = ctx,
                prep = prep,
                check_psd = identical(prep[["kernel"]], "precomputed")
            )
            prep
        },
        edge_case = function(ctx, prep) NULL,
        init = function(ctx, prep) {
            .aa_kernel_init_vars(
                G = prep[["gram"]],
                K = ctx[["K"]],
                init = ctx[["init"]],
                init_args = ctx[["init_args"]],
                eps = ctx[["eps"]],
                max_iter = ctx[["max_iter"]],
                verbose = ctx[["verbose"]],
                delta = delta
            )
        },
        fit = function(ctx, prep, init_vars) {
            .aa_fit_kernel_pgd(
                G = prep[["gram"]],
                weight_fun = .aa_weight_fun(ctx[["robust"]], ctx[["robust_args"]]),
                max_iter = ctx[["max_iter"]],
                tol = ctx[["tol"]],
                tol_r2 = ctx[["tol_r2"]],
                eps = ctx[["eps"]],
                verbose = ctx[["verbose"]],
                B = init_vars[["B"]],
                S = init_vars[["S"]],
                loss = init_vars[["loss"]],
                delta = delta,
                pseudo_pgd = pseudo_pgd,
                step_size = step_size,
                max_iter_optimizer = as.integer(max_iter_optimizer),
                step_shrinkage = step_shrinkage,
                max_no_update = as.integer(max_no_update)
            )
        },
        final_loss = .aa_final_loss,
        prepare_output = function(ctx, prep, fit) {
            .aa_prepare_kernel_output(
                call = ctx[["call"]],
                data = prep[["data"]],
                gram = prep[["gram"]],
                kernel = prep[["kernel"]],
                kernel_args = prep[["kernel_args"]],
                init = fit[["init"]],
                B = fit[["B"]],
                S = fit[["S"]],
                delta = fit[["delta"]],
                i = fit[["i"]],
                loss = fit[["loss"]],
                converged = fit[["converged"]],
                weights = fit[["row_weights"]],
                max_iter = ctx[["max_iter"]],
                verbose = ctx[["verbose"]],
                row_names = rownames(prep[["gram"]]),
                fit_info = list(
                    method = "kernel",
                    kernel = if (is.function(kernel)) "function" else kernel,
                    robust = !identical(ctx[["robust"]], FALSE),
                    robust_psi = .aa_robust_label(ctx[["robust"]]) %||% NA_character_,
                    robust_args = I(list(ctx[["robust_args"]]))
                )
            )
        }
    )
}

.aa_fit_kernel_pgd <- function(G,
                               weight_fun,
                               max_iter,
                               tol,
                               tol_r2,
                               eps,
                               verbose,
                               B,
                               S,
                               loss,
                               delta,
                               pseudo_pgd,
                               step_size,
                               max_iter_optimizer,
                               step_shrinkage,
                               max_no_update) {
    if (pseudo_pgd) {
        grad_S <- grad_S_l1
        grad_kernel_B <- grad_kernel_B_l1
        project <- proj_l1
    } else {
        grad_S <- grad_S_trace
        grad_kernel_B <- function(B, grad) grad
        project <- proj_simplex
    }

    update_alpha <- delta > 0
    a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
    a_hi <- 1 + delta
    clip <- function(a) pmax(pmin(a, a_hi), a_lo)

    init <- aB <- B # actual archetype coefficients used for loss and gradient computations
    a <- rowSums(aB)
    B <- aB / a # projected archetype coefficients used for updates
    slack_tol <- 1e-6
    if (any(a < a_lo - slack_tol) || any(a > a_hi + slack_tol)) {
        fmt <- "Initialize B marginals are outside the specified delta range [%.3f, %.3f]"
        stop(sprintf(fmt, a_lo, a_hi))
    }

    diagG <- diag(G)
    AAt <- aB %*% G %*% t(aB)
    XAt <- G %*% t(aB)
    row_rss <- pmax(diagG - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt)), 0)
    row_weights <- if (is.null(weight_fun)) NULL else weight_fun(row_rss)
    row_weights <- .aa_check_row_weights(row_weights, nrow(G))
    xss <- sum(.aa_weight_rows(diagG, row_weights))
    rss <- sum(.aa_weight_rows(row_rss, row_weights))
    S_weighted <- .aa_weight_rows(S, row_weights)
    StS <- crossprod(S_weighted, S)
    StG <- crossprod(S_weighted, G)

    loss[["loss"]][1L] <- rss
    loss[["r2"]][1L] <- 1 - rss / xss
    converged <- FALSE

    step_S <- step_size
    step_B <- step_size
    step_a <- step_size
    no_update <- 0L

    if (verbose) message("Starting kernel optimization loop...")
    for (i in seq_len(max_iter)) {
        accepted_update <- FALSE

        grad <- grad_S(S, AAt, XAt, row_weights = row_weights)
        for (k in seq_len(max_iter_optimizer)) {
            S_new <- project(S - step_S * grad, eps = eps)
            S_new_weighted <- .aa_weight_rows(S_new, row_weights)
            StS_new <- crossprod(S_new_weighted, S_new)
            rss_new <- xss - 2 * sum(S_new_weighted * XAt) + sum(StS_new * AAt)
            if (rss_new < rss) {
                S <- S_new
                S_weighted <- S_new_weighted
                StS <- StS_new
                StG <- crossprod(S_weighted, G)
                rss <- rss_new
                step_S <- step_S / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }

        grad_aB <- StS %*% aB %*% G - StG
        grad <- a * grad_kernel_B(B, grad_aB) # grad_B
        for (k in seq_len(max_iter_optimizer)) {
            B_new <- project(B - step_B * grad, eps = eps)
            aB_new <- a * B_new
            XAt_new <- G %*% t(aB_new)
            AAt_new <- aB_new %*% G %*% t(aB_new)
            rss_new <- xss - 2 * sum(S_weighted * XAt_new) + sum(StS * AAt_new)
            if (rss_new < rss) {
                B <- B_new
                aB <- aB_new
                XAt <- XAt_new
                AAt <- AAt_new
                rss <- rss_new
                step_B <- step_B / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_B <- step_B * step_shrinkage
        }

        if (update_alpha) {
            grad <- grad_alpha(B, grad_aB)
            for (k in seq_len(max_iter_optimizer)) {
                a_new <- clip(a - step_a * grad)
                a_update <- a_new / a
                aB_new <- a_update * aB
                XAt_new <- sweep(XAt, 2L, a_update, "*")
                AAt_new <- AAt * tcrossprod(a_update)
                rss_new <- xss - 2 * sum(S_weighted * XAt_new) + sum(StS * AAt_new)
                if (rss_new < rss) {
                    a <- a_new
                    aB <- aB_new
                    XAt <- XAt_new
                    AAt <- AAt_new
                    rss <- rss_new
                    step_a <- step_a / step_shrinkage
                    accepted_update <- TRUE
                    break
                }
                step_a <- step_a * step_shrinkage
            }
        }

        # Update loss
        row_rss <- pmax(diagG - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt)), 0)
        if (!is.null(weight_fun)) {
            row_weights <- weight_fun(row_rss)
            row_weights <- .aa_check_row_weights(row_weights, nrow(G))
            xss <- sum(.aa_weight_rows(diagG, row_weights))
            rss <- sum(.aa_weight_rows(row_rss, row_weights))
            S_weighted <- .aa_weight_rows(S, row_weights)
            StS <- crossprod(S_weighted, S)
            StG <- crossprod(S_weighted, G)
        }
        loss[["loss"]][i + 1L] <- rss
        loss[["r2"]][i + 1L] <- 1 - rss / xss

        # Check for accepted update; if none, shrink steps and check for stall
        if (!accepted_update) {
            no_update <- no_update + 1L
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
            if (no_update >= max_no_update) {
                warning(sprintf(
                    "Kernel PGD stalled: no line-search update accepted after %d consecutive step shrinkages",
                    max_no_update
                ), call. = FALSE)
                break
            }
            next
        }

        no_update <- 0L
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, verbose)
        if (converged) break
    }

    list(
        init = init,
        B = aB,
        S = S,
        delta = delta,
        i = i,
        loss = loss,
        row_weights = row_weights,
        converged = converged
    )
}

# TODO: extend to use kernlab or kerntools methods
.aa_kernel_prepare <- function(x, kernel, kernel_args, data = NULL) {
    if (!is.list(kernel_args)) {
        stop("`kernel_args` must be a list.", call. = FALSE)
    }

    if (identical(kernel, "precomputed")) {
        G <- as.matrix(x)
        attached_data <- if (is.data.frame(data)) data.matrix(data) else data
        return(list(
            gram = G,
            data = attached_data,
            kernel = "precomputed",
            kernel_args = list()
        ))
    }

    if (!is.null(data)) {
        stop("`data` is only used with `kernel = 'precomputed'`.", call. = FALSE)
    }
    stored_data <- if (is.data.frame(x)) data.matrix(x) else x
    if (is.function(kernel)) {
        G <- do.call(kernel, c(list(stored_data), kernel_args))
        return(list(gram = as.matrix(G), data = stored_data, kernel = kernel, kernel_args = kernel_args))
    }

    if (!is.character(kernel) || length(kernel) != 1L) {
        stop("`kernel` must be a single string or a function.", call. = FALSE)
    }
    kernel <- match.arg(kernel, c("linear", "rbf", "laplace", "polynomial", "precomputed"))
    X <- if (inherits(stored_data, "sparseMatrix") && kernel %in% c("linear", "polynomial")) {
        stored_data
    } else {
        as.matrix(stored_data)
    }
    G <- do.call(.aa_builtin_kernel, c(list(X = X, kernel = kernel), kernel_args))
    list(gram = G, data = stored_data, kernel = kernel, kernel_args = kernel_args)
}

.aa_builtin_kernel <- function(X,
                               kernel,
                               sigma = NULL,
                               gamma = NULL,
                               degree = 3,
                               coef0 = 1) {
    # Catch any mis-specification of kernel arguments that would be ignored by the kernel function
    if (!is.null(gamma) && kernel %in% c("rbf", "laplace")) {
        stop(
            sprintf(
                "%s kernels use `sigma`; `gamma` is only used for polynomial kernels",
                ifelse(kernel == "rbf", "RBF", "Laplace")
            ),
            call. = FALSE
        )
    } else if (!is.null(gamma) && kernel != "polynomial") {
        warning(
            "`gamma` is only used for the polynomial kernel; ignoring for other kernels",
            call. = FALSE
        )
    } else if (!is.null(sigma) && !(kernel %in% c("rbf", "laplace"))) {
        warning(
            "`sigma` is only used for RBF and Laplace kernels; ignoring for other kernels",
            call. = FALSE
        )
    }

    switch(kernel,
        linear = tcrossprod(X),
        rbf = .aa_rbf_kernel(X, sigma = sigma),
        laplace = .aa_laplace_kernel(X, sigma = sigma),
        polynomial = .aa_polynomial_kernel(X, gamma, degree, coef0)
    )
}

.aa_auto_rbf_sigma <- function(D) {
    d_upper <- if (is.matrix(D)) D[upper.tri(D)] else as.vector(D)
    d_90 <- d_upper[d_upper <= stats::quantile(d_upper, 0.9)]
    sigma <- .aa_otsu_threshold(d_90) * 0.75
    if (!is.finite(sigma) || sigma <= 0) 1 else sigma
}

.aa_rbf_kernel <- function(X, sigma = NULL) {
    D <- stats::dist(X)
    sigma <- sigma %||% .aa_auto_rbf_sigma(D)
    if (!is_positive(sigma)) {
        stop("`sigma` must be a positive finite number.", call. = FALSE)
    }
    G <- exp(-0.5 * (D / sigma)^2)
    # convert to full symmetric matrix at the end to avoid redundant computations
    G <- as.matrix(G)
    diag(G) <- 1
    G
}

.aa_laplace_kernel <- function(X, sigma = NULL) {
    D <- stats::dist(X, method = "manhattan")
    sigma <- sigma %||% .aa_auto_rbf_sigma(D)
    if (!is_positive(sigma)) {
        stop("`sigma` must be a positive finite number.", call. = FALSE)
    }
    G <- exp(-D / sigma)
    # convert to full symmetric matrix at the end to avoid redundant computations
    G <- as.matrix(G)
    diag(G) <- 1
    G
}

.aa_polynomial_kernel <- function(X, gamma, degree, coef0) {
    if (!is_positive(degree)) {
        stop("`degree` must be a positive finite number.", call. = FALSE)
    }
    if (!is_number(coef0)) {
        stop("`coef0` must be a finite number.", call. = FALSE)
    }
    gamma <- gamma %||% (1 / ncol(X))
    if (!is_positive(gamma)) {
        stop("`gamma` must be a positive finite number.", call. = FALSE)
    }
    (gamma * tcrossprod(X) + coef0)^degree
}

.aa_check_kernel_inputs <- function(ctx, prep, check_psd = FALSE) {
    G <- prep[["gram"]]
    data <- prep[["data"]]
    if (nrow(G) != ncol(G)) {
        stop("Gram matrix must be square.", call. = FALSE)
    }
    if (!is.numeric(G)) {
        stop("Gram matrix must be numeric.", call. = FALSE)
    }
    if (any(is.na(G)) || !all(is.finite(G))) {
        stop("Gram matrix contains missing or non-finite values.", call. = FALSE)
    }
    if (!isSymmetric(G, tol = 1e-8)) {
        stop("Gram matrix must be symmetric.", call. = FALSE)
    }
    if (check_psd) {
        eigenvalues <- eigen(G, symmetric = TRUE, only.values = TRUE)[["values"]]
        tol_psd <- 1e-8 * max(1, max(abs(eigenvalues)))
        if (min(eigenvalues) < -tol_psd) {
            stop("Gram matrix must be positive semidefinite.", call. = FALSE)
        }
    }
    .aa_check_fit_controls(ctx, n = nrow(G))
    if (!is.null(data)) {
        if (nrow(data) != nrow(G)) {
            stop("`data` rows must match Gram matrix dimensions.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

.aa_kernel_init_vars <- function(G, K, init, init_args, eps, max_iter, verbose, delta = 0) {
    if (verbose) message("Initializing kernel archetypes...")
    L <- max_iter + 1L

    if (is_tabular(init)) {
        B <- as.matrix(init)
        if (nrow(B) != K) {
            fmt <- "nrow(init) = %d does not match K (%d)"
            stop(sprintf(fmt, nrow(B), K))
        }
        if (ncol(B) != nrow(G)) {
            fmt <- "ncol(init) = %d does not match number of samples (%d)"
            stop(sprintf(fmt, ncol(B), nrow(G)))
        }
        if (!all(B >= 0)) {
            stop("Initial coefficient matrix must be non-negative.", call. = FALSE)
        }
        a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
        a_hi <- 1 + delta
        a <- rowSums(B)
        if (!all(a >= a_lo - 1e-8 & a <= a_hi + 1e-8)) {
            stop("Initial coefficient row sums are outside the allowed delta range.", call. = FALSE)
        }
        nm <- .aa_init_names(B)
        rownames(B) <- nm
        colnames(B) <- rownames(G)
        S <- .aa_kernel_init_S(G, B, eps)
        return(list(
            init = B,
            B = B,
            S = S,
            loss = list(loss = rep(NA_real_, L), r2 = rep(NA_real_, L))
        ))
    }

    if (is_single_string(init) && init == "dirichlet") {
        B <- do.call(
            .aa_kernel_dirichlet,
            args = c(list(G = G, K = K), init_args)
        )
        nm <- .aa_init_names(B)
        rownames(B) <- nm
        colnames(B) <- rownames(G)
        S <- .aa_kernel_init_S(G, B, eps)
        return(list(
            init = B,
            B = B,
            S = S,
            loss = list(loss = rep(NA_real_, L), r2 = rep(NA_real_, L))
        ))
    }

    # TODO: disallow "row selection" init
    if (is_single_string(init) && init %in% c("random", "furthest_first", "kmeans_pp", "furthest_sum")) {
        method <- match.arg(init, c("random", "furthest_first", "kmeans_pp", "furthest_sum"))
        ind <- do.call(
            .aa_kernel_init_indices,
            args = c(list(G = G, K = K, method = method), init_args)
        )
    } else if (is.numeric(init) || is.character(init) || is.logical(init)) {
        ind <- .aa_normalize_row_indices(init, nrow(G), rownames(G))
    } else {
        stop("`init` must be a method string, row indices/names, or a coefficient matrix",
            call. = FALSE
        )
    }
    if (length(ind) != K) {
        fmt <- "length(init) = %d does not match K (%d)"
        stop(sprintf(fmt, length(ind), K))
    }

    B <- onehot(ind, sparse = FALSE, nc = nrow(G))
    nm <- names(ind) %||% paste0("A", seq_along(ind))
    rownames(B) <- nm
    colnames(B) <- rownames(G)
    S <- .aa_kernel_init_S(G, B, eps)
    list(
        init = B,
        B = B,
        S = S,
        loss = list(loss = rep(NA_real_, L), r2 = rep(NA_real_, L))
    )
}

.aa_kernel_init_indices <- function(G,
                                    K,
                                    method = "furthest_sum",
                                    batch_size = NULL,
                                    batch_type = c("distal", "uniform"),
                                    batch_replace = FALSE,
                                    ...) {
    batch_type <- match.arg(batch_type)
    if (!is_logical(batch_replace)) {
        stop("`batch_replace` must be TRUE or FALSE.", call. = FALSE)
    }
    batch_size <- .aa_validate_batch_size(
        batch_size,
        n = nrow(G),
        K = K,
        replace = batch_replace
    )
    distances <- .aa_kernel_dist2(G)
    center_dists <- .aa_kernel_dist2_center(G)
    candidates <- .aa_sample(
        center_dists,
        size = batch_size,
        type = batch_type,
        replace = batch_replace
    )
    candidate_distances <- distances[candidates, candidates, drop = FALSE]
    candidate_center_dists <- center_dists[candidates]
    G_candidates <- G[candidates, candidates, drop = FALSE]

    ind <- switch(method,
        random = sample(seq_along(candidates), K, replace = FALSE),
        furthest_first = furthest_first(
            G_candidates,
            K,
            distances = candidate_distances,
            center_dists = candidate_center_dists
        ),
        kmeans_pp = kmeans_pp(
            G_candidates,
            K,
            distances = candidate_distances,
            center_dists = candidate_center_dists
        ),
        furthest_sum = furthest_sum(G_candidates, K, distances = candidate_distances)
    )
    candidates[ind]
}

.aa_kernel_dirichlet <- function(G,
                                 K,
                                 alpha = 1,
                                 batch_size = NULL,
                                 batch_type = c("distal", "uniform"),
                                 batch_replace = FALSE,
                                 ...) {
    if (!is_positive(alpha)) {
        stop("`alpha` must be a single positive number.", call. = FALSE)
    }
    batch_type <- match.arg(batch_type)
    if (!is_logical(batch_replace)) {
        stop("`batch_replace` must be TRUE or FALSE.", call. = FALSE)
    }
    batch_size <- .aa_validate_batch_size(
        batch_size,
        n = nrow(G),
        K = K,
        replace = batch_replace
    )
    center_dists <- .aa_kernel_dist2_center(G)
    candidates <- .aa_sample(
        center_dists,
        size = batch_size,
        type = batch_type,
        replace = batch_replace
    )
    B <- matrix(0, nrow = K, ncol = nrow(G))
    B[, candidates] <- stats::rgamma(K * length(candidates), shape = alpha)
    B <- B / rowSums(B)
    rownames(B) <- paste0("A", seq_len(K))
    colnames(B) <- rownames(G)
    B
}

.aa_kernel_dist2 <- function(G, ind = NULL) {
    d <- diag(G)
    ind <- ind %||% seq_len(ncol(G))
    pmax(outer(d, d[ind], "+") - 2 * G[, ind, drop = FALSE], 0)
}

.aa_kernel_dist2_center <- function(G) {
    N <- nrow(G)
    row_mean <- rowMeans(G)
    center_norm <- sum(G) / (N * N)
    pmax(diag(G) - 2 * row_mean + center_norm, 0)
}

.aa_kernel_init_S <- function(G, B, eps = 0) {
    AAt <- B %*% G %*% t(B)
    XAt <- G %*% t(B)
    D2 <- pmax(outer(diag(G), diag(AAt), "+") - 2 * XAt, 0)
    S <- proj_l1(1 / D2, eps = eps)
    S[is.nan(S)] <- 1
    colnames(S) <- rownames(B)
    rownames(S) <- rownames(G)
    S
}

.aa_prepare_kernel_output <- function(call, data, gram, kernel, kernel_args, init,
                                      B, S, delta, i, loss, converged,
                                      weights = NULL,
                                      max_iter, verbose, row_names = NULL,
                                      fit_info = list()) {
    j <- i + 1L
    loss <- as.data.frame(loss)[seq_len(j), , drop = FALSE]
    rownames(loss) <- NULL

    archetype_names <- rownames(B) %||% paste0("A", seq_len(nrow(B)))
    sample_names <- row_names %||% paste0("x", seq_len(ncol(B)))

    rownames(B) <- colnames(S) <- archetype_names
    colnames(B) <- rownames(S) <- rownames(gram) <- colnames(gram) <- sample_names
    if (!is.null(init)) {
        rownames(init) <- archetype_names
    }

    if (!converged) {
        warning(sprintf("Algorithm did not converge after %d iterations", max_iter))
    }
    if (verbose) {
        fmt <- ifelse(converged, "Converged after %d iterations:", "Final iteration %d:")
        fmt <- paste(fmt, "loss = %.4g, R2 = %.3f")
        message(sprintf(fmt, i, loss[j, "loss"], loss[j, "r2"]))
    }

    .aa_new_kernel_archetypes(
        coefficients = B,
        compositions = S,
        gram = gram,
        slack = delta,
        loss = loss,
        converged = converged,
        call = call,
        data = data,
        init = init,
        kernel = kernel,
        kernel_args = kernel_args,
        weights = weights,
        fit_info = fit_info %|p|% list(
            family = "gaussian",
            robust = fit_info[["robust"]] %||% !is.null(weights),
            robust_psi = fit_info[["robust_psi"]] %||% NA_character_,
            robust_args = fit_info[["robust_args"]] %||% I(list(list())),
            missing = FALSE,
            delta = delta,
            init = "kernel coefficients",
            scaling = "none"
        )
    )
}

grad_kernel_B_l1 <- function(B, grad) {
    row_dot <- rowSums(B * grad)
    grad - row_dot
}
