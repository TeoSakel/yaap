#' Kernel Archetypal Analysis using Projected Gradient Descent
#'
#' Fits kernel archetypal analysis from a precomputed Gram matrix or from a
#' kernel computed on `data`. For nonlinear kernels, `coordinates_proxy` is the
#' input-space convex combination `coefficients %*% data`, matching the plotting
#' convention used by Mørup and Hansen for Figure 5. It is a visualization
#' proxy, not the exact Hilbert-space archetype.
#'
#' @param data Optional data matrix with rows as samples.
#' @param K Number of archetypes.
#' @param gram Optional precomputed `N x N` Gram matrix.
#' @param kernel Optional kernel specification. One of `"linear"`, `"rbf"`,
#'   `"polynomial"`, or a function returning an `N x N` Gram matrix.
#' @param kernel_args List of arguments passed to the kernel.
#' @param init Initialization method, row indices/names, or a `K x N`
#'   coefficient matrix.
#' @param init_args List of additional arguments for the initialization method.
#' @param robust Whether to use Tukey bisquare row reweighting.
#' @param tukey_c Tuning constant for Tukey bisquare weights.
#' @param max_iter Maximum number of outer iterations.
#' @param tol Convergence tolerance based on residual sum of squares.
#' @param tol_r2 Convergence tolerance based on R-squared.
#' @param max_kappa Maximum condition number warning threshold.
#' @param eps Small positive number to ensure numerical stability.
#' @param verbose Whether to print progress messages.
#' @param delta Maximum allowed relaxation of archetype convexity constraint.
#' @param pseudo_pgd Whether to use pseudo projected gradient descent.
#' @param step_size Initial line-search step size.
#' @param max_iter_optimizer Maximum line-search iterations per update.
#' @param step_shrinkage Factor used to shrink rejected line-search steps.
#' @param max_no_update Maximum consecutive outer iterations with no accepted
#'   line-search update before stopping as stalled.
#'
#' @returns An object of class `kernel_archetypes`.
#'
#' @export
archetypes_kernel_pgd <- function(data = NULL,
                                  K,
                                  gram = NULL,
                                  kernel = NULL,
                                  kernel_args = list(),
                                  init = "furthest_sum",
                                  init_args = list(),
                                  robust = FALSE,
                                  tukey_c = 4.685,
                                  max_iter = 100L,
                                  tol = 1e-6,
                                  tol_r2 = 0.9999,
                                  max_kappa = 1000,
                                  eps = 1e-8,
                                  verbose = FALSE,
                                  delta = 0,
                                  pseudo_pgd = TRUE,
                                  step_size = 1.0,
                                  max_iter_optimizer = 10L,
                                  step_shrinkage = 0.5,
                                  max_no_update = 5L) {
    call <- match.call()
    kernel_spec <- .aa_kernel_prepare(data, gram, kernel, kernel_args)
    G <- kernel_spec[["gram"]]
    N <- nrow(G)

    .aa_check_kernel_inputs(
        G = G,
        data = data,
        K = K,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        robust = robust,
        tukey_c = tukey_c
    )
    pgd_args <- .aa_pgd_method_args(
        delta = delta,
        pseudo_pgd = pseudo_pgd,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = max_no_update
    )

    weight_fun <- if (robust) {
        function(row_rss) .aa_bisquare_weights(row_rss, c = tukey_c)
    } else {
        NULL
    }

    init_vars <- .aa_kernel_init_vars(
        G = G,
        K = K,
        init = init,
        init_args = init_args,
        eps = eps,
        max_iter = max_iter,
        verbose = verbose,
        delta = pgd_args[["delta"]]
    )

    fit <- .aa_fit_kernel_pgd(
        G = G,
        weight_fun = weight_fun,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        B = init_vars[["B"]],
        S = init_vars[["S"]],
        loss = init_vars[["loss"]],
        delta = pgd_args[["delta"]],
        pseudo_pgd = pgd_args[["pseudo_pgd"]],
        step_size = pgd_args[["step_size"]],
        max_iter_optimizer = pgd_args[["max_iter_optimizer"]],
        step_shrinkage = pgd_args[["step_shrinkage"]],
        max_no_update = pgd_args[["max_no_update"]]
    )

    .aa_prepare_kernel_output(
        call = call,
        data = data,
        gram = G,
        kernel = kernel_spec[["kernel"]],
        kernel_args = kernel_spec[["kernel_args"]],
        init = init_vars[["init"]],
        B = fit[["B"]],
        S = fit[["S"]],
        delta = fit[["delta"]],
        i = fit[["i"]],
        loss = fit[["loss"]],
        converged = fit[["converged"]],
        max_iter = max_iter,
        verbose = verbose,
        row_names = rownames(G)
    )
}

.aa_fit_kernel_pgd <- function(G,
                               weight_fun,
                               max_iter,
                               tol,
                               tol_r2,
                               max_kappa,
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
        grad_S        <- grad_S_l1
        grad_kernel_B <- grad_kernel_B_l1
        project       <- proj_l1
    } else {
        grad_S        <- grad_S_trace
        grad_kernel_B <- function(B, grad) grad
        project       <- proj_simplex
    }

    update_alpha <- delta > 0
    a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
    a_hi <- 1 + delta
    clip <- function(a) pmax(pmin(a, a_hi), a_lo)

    init <- aB <- B  # actual archetype coefficients used for loss and gradient computations
    a <- rowSums(aB)
    B <- aB / a  # projected archetype coefficients used for updates
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
    if (!is.null(row_weights)) {
        .aa_check_row_weights(row_weights, nrow(G))
        if (.aa_trivial_row_weights(row_weights))
            row_weights <- NULL
    }
    xss <- sum(.aa_weight_rows(diagG, row_weights))
    rss <- sum(.aa_weight_rows(row_rss, row_weights))
    S_weighted <- .aa_weight_rows(S, row_weights)
    StS <- crossprod(S_weighted, S)
    StG <- crossprod(S_weighted, G)
    loss_terms <- list(rss = rss, xss = xss, StS = StS, AAt = AAt)

    loss <- .aa_update_loss(
        loss, 1L, loss_terms, verbose = verbose, max_kappa = max_kappa,
        k_A = "gram"
    )
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
        grad <- a * grad_kernel_B(B, grad_aB)  # grad_B
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
            .aa_check_row_weights(row_weights, nrow(G))
            if (.aa_trivial_row_weights(row_weights))
                row_weights <- NULL
            xss <- sum(.aa_weight_rows(diagG, row_weights))
            rss <- sum(.aa_weight_rows(row_rss, row_weights))
            S_weighted <- .aa_weight_rows(S, row_weights)
            StS <- crossprod(S_weighted, S)
            StG <- crossprod(S_weighted, G)
        }
        loss_terms <- list(rss = rss, xss = xss, StS = StS, AAt = AAt)
        loss <- .aa_update_loss(
            loss, i + 1L, loss_terms, verbose = verbose, max_kappa = max_kappa,
            k_A = "gram"
        )

        # Check for accepted update; if none, shrink steps and check for stall
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
                warning(sprintf(
                    "Kernel PGD stalled: no line-search update accepted after %d consecutive step shrinkages",
                    max_no_update
                ), call. = FALSE)
                break
            }
            next
        }

        no_update <- 0L
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, max_kappa, verbose)
        if (converged) break
    }

    list(
        init = init,
        B = aB,
        S = S,
        delta = delta,
        i = i,
        loss = loss,
        converged = converged
    )
}

.aa_kernel_prepare <- function(data, gram, kernel, kernel_args) {
    has_gram <- !is.null(gram)
    has_kernel <- !is.null(kernel)
    if (has_gram && has_kernel)
        stop("Supply exactly one of `gram` or `kernel`, not both.", call. = FALSE)
    if (!has_gram && !has_kernel)
        stop("Supply either a precomputed `gram` matrix or a `kernel`.", call. = FALSE)
    if (has_kernel && is.null(data))
        stop("`data` is required when `kernel` is supplied.", call. = FALSE)
    if (!is.list(kernel_args))
        stop("`kernel_args` must be a list.", call. = FALSE)

    if (has_gram) {
        G <- as.matrix(gram)
        return(list(gram = G, kernel = "precomputed", kernel_args = list()))
    }

    X <- as.matrix(data)
    if (is.function(kernel)) {
        G <- do.call(kernel, c(list(X), kernel_args))
        return(list(gram = as.matrix(G), kernel = kernel, kernel_args = kernel_args))
    }

    # TODO: extend to use kernlab or kerntools methods
    if (!is.character(kernel) || length(kernel) != 1L)
        stop("`kernel` must be a single string or a function.", call. = FALSE)
    kernel <- match.arg(kernel, c("linear", "rbf", "laplace", "polynomial"))
    G <- do.call(.aa_builtin_kernel, c(list(X = X, kernel = kernel), kernel_args))
    list(gram = G, kernel = kernel, kernel_args = kernel_args)
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
            sprintf("%s kernels use `sigma`; `gamma` is only used for polynomial kernels",
                    ifelse(kernel == "rbf", "RBF", "Laplace")),
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

    switch(
        kernel,
        linear = tcrossprod(X),
        rbf = rbf_kernel(X, sigma = sigma),
        laplace = laplace_kernel(X, sigma = sigma),
        polynomial = polynomial_kernel(X, gamma, degree, coef0)
    )
}

.auto_rbf_sigma <- function(D) {
    d_upper <- if (is.matrix(D)) D[upper.tri(D)] else as.vector(D)
    d_90    <- d_upper[d_upper <= stats::quantile(d_upper, 0.9)]
    sigma   <- .otsu_threshold(d_90) * 0.75
    if (!is.finite(sigma) || sigma <= 0) 1 else sigma
}

rbf_kernel <- function(X, sigma = NULL) {
    D <- stats::dist(X)
    if (is.null(sigma))
        sigma <- .auto_rbf_sigma(D)
    stopifnot("`sigma` must be a positive finite number" =
                  length(sigma) == 1L && is.finite(sigma) && sigma > 0)
    G <- exp(-0.5 * (D / sigma)^2)
    # convert to full symmetric matrix at the end to avoid redundant computations
    G <- as.matrix(G)
    diag(G) <- 1
    G
}

laplace_kernel <- function(X, sigma = NULL) {
    D <- stats::dist(X, method = "manhattan")
    if (is.null(sigma))
        sigma <- .auto_rbf_sigma(D)
    stopifnot("`sigma` must be a positive finite number" =
                  length(sigma) == 1L && is.finite(sigma) && sigma > 0)
    G <- exp(-D / sigma)
    # convert to full symmetric matrix at the end to avoid redundant computations
    G <- as.matrix(G)
    diag(G) <- 1
    G
}

polynomial_kernel <- function(X, gamma, degree, coef0) {
    stopifnot("`degree` must be a positive finite number" =
                  length(degree) == 1L && is.finite(degree) && degree > 0)
    stopifnot("`coef0` must be a finite number" = length(coef0) == 1L && is.finite(coef0))
    if (is.null(gamma))
        gamma <- 1 / ncol(X)
    stopifnot("`gamma` must be a positive finite number" =
                  length(gamma) == 1L && is.finite(gamma) && gamma > 0)
    (gamma * tcrossprod(X) + coef0)^degree
}

.aa_check_kernel_inputs <- function(G, data, K, tol, tol_r2, max_kappa,
                                    eps, robust, tukey_c) {
    stopifnot("Gram matrix must be square" = nrow(G) == ncol(G))
    stopifnot("Gram matrix must be numeric" = is.numeric(G))
    stopifnot("Gram matrix contains missing or non-finite values" =
                  !any(is.na(G)) && all(is.finite(G)))
    stopifnot("Gram matrix must be symmetric" = isSymmetric(G, tol = 1e-8))
    .aa_check_inputs(
        data = G,
        K = K,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        robust = robust,
        tukey_c = tukey_c,
        scale = FALSE
    )
    if (!is.null(data))
        stopifnot("`data` rows must match Gram matrix dimensions" = nrow(data) == nrow(G))
    invisible(TRUE)
}

.aa_kernel_init_vars <- function(G, K, init, init_args, eps, max_iter, verbose, delta = 0) {
    if (verbose) message("Initializing kernel archetypes...")
    L <- max_iter + 1L

    if (is.matrix(init) || inherits(init, "data.frame")) {
        B <- as.matrix(init)
        if (nrow(B) != K) {
            fmt <- "nrow(init) = %d does not match K (%d)"
            stop(sprintf(fmt, nrow(B), K))
        }
        if (ncol(B) != nrow(G)) {
            fmt <- "ncol(init) = %d does not match number of samples (%d)"
            stop(sprintf(fmt, ncol(B), nrow(G)))
        }
        stopifnot("init coefficient matrix must be non-negative" = all(B >= 0))
        a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
        a_hi <- 1 + delta
        a <- rowSums(B)
        stopifnot("init coefficient row sums are outside the allowed delta range" =
                      all(a >= a_lo - 1e-8 & a <= a_hi + 1e-8))
        nm <- .aa_init_names(B)
        rownames(B) <- nm
        colnames(B) <- rownames(G)
        S <- .aa_kernel_init_S(G, B, eps)
        return(list(
            init = B,
            B = B,
            S = S,
            loss = .aa_new_loss(L)
        ))
    }

    if (is.character(init) && length(init) == 1L &&
            init %in% c("uniform_archetypes", "furthest_first", "kmeans_pp", "furthest_sum")) {
        method <- match.arg(init, c("uniform_archetypes", "furthest_first", "kmeans_pp", "furthest_sum"))
        ind <- do.call(
            .aa_kernel_init_indices,
            args = c(list(G = G, K = K, method = method), init_args)
        )
    } else if (is.numeric(init) || is.character(init) || is.logical(init)) {
        ind <- .aa_normalize_row_indices(init, nrow(G), rownames(G))
    } else {
        stop("`init` must be a method string, row indices/names, or a coefficient matrix",
             call. = FALSE)
    }
    if (length(ind) != K) {
        fmt <- "length(init) = %d does not match K (%d)"
        stop(sprintf(fmt, length(ind), K))
    }

    B <- onehot(ind, sparse = FALSE, nc = nrow(G))
    nm <- names(ind)
    if (is.null(nm))
        nm <- paste0("A", seq_along(ind))
    rownames(B) <- nm
    colnames(B) <- rownames(G)
    S <- .aa_kernel_init_S(G, B, eps)
    list(
        init = B,
        B = B,
        S = S,
        loss = .aa_new_loss(L)
    )
}

.aa_kernel_init_indices <- function(G, K, method = "furthest_sum", ...) {
    distances <- .aa_kernel_dist2(G)
    center_dists <- .aa_kernel_dist2_center(G)
    switch(
        method,
        uniform_archetypes = sample(nrow(G), K, replace = FALSE),
        furthest_first = furthest_first(G, K, distances = distances, center_dists = center_dists),
        kmeans_pp = kmeans_pp(G, K, distances = distances, center_dists = center_dists),
        furthest_sum = furthest_sum(G, K, distances = distances)
    )
}

.aa_kernel_dist2 <- function(G, ind = NULL) {
    d <- diag(G)
    if (is.null(ind))
        ind <- seq_len(ncol(G))
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
                                      max_iter, verbose, row_names = NULL) {
    j <- i + 1L
    loss <- as.data.frame(loss)[seq_len(j), , drop = FALSE]
    rownames(loss) <- NULL

    archetype_names <- rownames(B)
    if (is.null(archetype_names))
        archetype_names <- paste0("A", seq_len(nrow(B)))
    sample_names <- row_names
    if (is.null(sample_names))
        sample_names <- paste0("x", seq_len(ncol(B)))

    rownames(B) <- colnames(S) <- archetype_names
    colnames(B) <- rownames(S) <- rownames(gram) <- colnames(gram) <- sample_names
    if (!is.null(init))
        rownames(init) <- archetype_names

    coordinates_proxy <- NULL
    if (!is.null(data)) {
        X <- as.matrix(data)
        coordinates_proxy <- B %*% X
        rownames(coordinates_proxy) <- archetype_names
        colnames(coordinates_proxy) <- colnames(X)
    }

    if (!converged)
        warning(sprintf("Algorithm did not converge after %d iterations", max_iter))
    if (verbose) {
        fmt <- ifelse(converged, "Converged after %d iterations:", "Final iteration %d:")
        fmt <- paste(fmt, "loss = %.4g, R2 = %.3f")
        message(sprintf(fmt, i, loss[j, "loss"], loss[j, "r2"]))
    }

    kernel_archetypes(
        coefficients = B,
        compositions = S,
        gram = gram,
        coordinates_proxy = coordinates_proxy,
        slack = delta,
        loss = loss,
        converged = converged,
        call = call,
        data = data,
        init = init,
        kernel = kernel,
        kernel_args = kernel_args
    )
}

grad_kernel_B_l1 <- function(B, grad) {
    row_dot <- rowSums(B * grad)
    grad - row_dot
}
