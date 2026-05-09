#' Directional Archetypal Analysis
#'
#' Fits a matrix-first directional archetypal analysis model. Directional AA
#' models observations by their direction rather than Euclidean magnitude:
#' archetypes are convex combinations of row-normalized input data, optionally
#' flipped onto a common hemisphere, and sample reconstructions are evaluated by
#' a Watson-style polarity-invariant angular loss.
#'
#' @param data Dense numeric data matrix (rows = samples, columns = dimensions).
#' @param K Number of archetypes.
#' @param init Initialization method. `"random"` samples exponential generator
#'   weights as in the reference implementation; any method accepted by
#'   [aa_init()] is also supported. A numeric `K x ncol(data)` matrix may be
#'   supplied as initial archetype directions.
#' @param init_args List of additional arguments for the initialization method.
#' @param weights Optional non-negative sample weights.
#' @param max_iter Maximum number of outer iterations.
#' @param tol Convergence tolerance based on directional residual loss.
#' @param tol_r2 Convergence tolerance based on directional R-squared.
#' @param max_kappa Accepted for consistency with other fitters. Directional AA
#'   does not solve by matrix inversion, so condition numbers are not computed.
#' @param eps Small positive number for numerical stability.
#' @param verbose Whether to print progress messages.
#' @param hemisphere Hemisphere handling. `"pca"` flips generator rows onto the
#'   dominant PCA hemisphere; `"none"` leaves signs unchanged.
#' @param precision Precision weighting. `"row_norm"` keeps row magnitudes in
#'   the Watson loss, matching Olsen et al.; `"unit"` gives every row equal
#'   precision.
#' @param step_size Initial projected-gradient step size.
#' @param max_iter_optimizer Maximum line-search iterations per update.
#' @param step_shrinkage Factor used to shrink rejected line-search steps.
#' @param max_no_update Maximum consecutive outer iterations with no accepted
#'   line-search update before stopping as stalled.
#'
#' @returns An object of class `directional_archetypes`, extending
#'   \code{\link{archetypes}}.
#'
#' @references
#' Olsen, A. S., Høegh, R. M. T., Hinrich, J. L., Madsen, K. H., & Mørup, M.
#' (2022). Combining electro- and magnetoencephalography data using directional
#' archetypal analysis. *Frontiers in Neuroscience*, 16, 911034.
#' \doi{10.3389/fnins.2022.911034}
#'
#' @examples
#' theta <- seq(0, pi / 2, length.out = 50)
#' X <- cbind(cos(theta), sin(theta))
#' fit <- archetypes_directional(X, K = 3, max_iter = 5)
#'
#' @export
archetypes_directional <- function(data,
                                   K,
                                   init = "random",
                                   init_args = list(),
                                   weights = NULL,
                                   max_iter = 100L,
                                   tol = 1e-6,
                                   tol_r2 = 0.9999,
                                   max_kappa = 1000,
                                   eps = 1e-8,
                                   verbose = FALSE,
                                   hemisphere = c("pca", "none"),
                                   precision = c("row_norm", "unit"),
                                   step_size = 1.0,
                                   max_iter_optimizer = 10L,
                                   step_shrinkage = 0.5,
                                   max_no_update = 5L) {
    call <- match.call()
    hemisphere <- match.arg(hemisphere)
    precision <- match.arg(precision)
    .aa_check_directional_inputs(
        data = data,
        K = K,
        weights = weights,
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

    X <- as.matrix(data)
    row_norms <- sqrt(rowSums(X * X))

    # Nomenclature follows Olsen et al. (2022), transposed to YAAAP's
    # row-oriented convention:
    #   X_loss: paper X^T, retaining row magnitude as Watson precision.
    #   X_gen: paper X_tilde^T, row-normalized generator data.
    #   X_flip: paper X_tilde_f^T, generator data flipped to one hemisphere.
    X_gen <- .aa_unit_rows(X, row_norms)
    X_loss <- if (precision == "row_norm") X else X_gen
    if (!is.null(weights)) {
        weights <- weights / mean(weights)
        X_loss <- X_loss * sqrt(weights)
    }
    flip <- .aa_directional_hemisphere(X_gen, method = hemisphere)
    X_flip <- flip[["X"]]

    init_vars <- .aa_directional_init_vars(
        X_flip = X_flip,
        X_gen = X_gen,
        K = K,
        init = init,
        init_args = init_args,
        eps = eps
    )

    fit <- .aa_fit_directional(
        X_loss = X_loss,
        X_flip = X_flip,
        A0 = init_vars[["A0"]],
        B = init_vars[["B"]],
        S = init_vars[["S"]],
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        eps = eps,
        verbose = verbose,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = as.integer(max_no_update)
    )

    directional_archetypes(
        coordinates = fit[["A"]],
        coefficients = fit[["B"]],
        compositions = fit[["S"]],
        loss = fit[["loss"]],
        converged = fit[["converged"]],
        call = call,
        data = X,
        init = init_vars[["A0"]],
        generator_data = X_flip,
        hemisphere_direction = flip[["direction"]],
        row_norms = row_norms,
        precision = precision
    )
}

.aa_check_directional_inputs <- function(data,
                                         K,
                                         weights,
                                         max_iter,
                                         tol,
                                         tol_r2,
                                         max_kappa,
                                         eps,
                                         step_size,
                                         max_iter_optimizer,
                                         step_shrinkage,
                                         max_no_update) {
    if (inherits(data, "sparseMatrix"))
        stop("Directional AA currently requires a dense numeric matrix.", call. = FALSE)
    stopifnot("data must be a matrix-like object" =
                  is.matrix(data) || inherits(data, "data.frame"))
    data <- as.matrix(data)
    stopifnot("data must be numeric" = is.numeric(data))
    stopifnot("data contains missing values" = !any(is.na(data)))
    stopifnot("data contains non-finite values" = all(is.finite(data)))
    stopifnot("K must be an integer" = K == as.integer(K))
    stopifnot("K must be an integer greater or equal to 1" = K >= 1L)
    stopifnot("K cannot be greater than number of samples" = K <= nrow(data))
    stopifnot("max_iter must be a non-negative integer" =
                  max_iter == as.integer(max_iter) && max_iter >= 0L)
    stopifnot("tol must be positive" = tol > 0)
    stopifnot("tol_r2 must be between (0, 1)" = tol_r2 >= 0 && tol_r2 <= 1)
    stopifnot("max_kappa must be >=1" = max_kappa >= 1)
    stopifnot("eps must be non-negative" = eps >= 0)
    stopifnot("step_size must be positive" = step_size > 0)
    stopifnot("max_iter_optimizer must be a positive integer" =
                  max_iter_optimizer == as.integer(max_iter_optimizer) &&
                      max_iter_optimizer >= 1L)
    stopifnot("step_shrinkage must be between (0, 1)" =
                  step_shrinkage > 0 && step_shrinkage < 1)
    stopifnot("max_no_update must be a positive integer" =
                  max_no_update == as.integer(max_no_update) && max_no_update >= 1L)
    if (!is.null(weights)) {
        stopifnot("weights must match rows in data" = length(weights) == nrow(data))
        stopifnot("weights contain NA values" = !any(is.na(weights)))
        stopifnot("weights must be non-negative" = all(weights >= 0))
        stopifnot("at least one weight must be positive" = any(weights > 0))
    }
    .aa_check_no_zero_rows(data)
    invisible(TRUE)
}

.aa_check_no_zero_rows <- function(X) {
    zero <- which(sqrt(rowSums(X * X)) <= .Machine$double.eps)
    if (length(zero) > 0L) {
        fmt <- "Directional AA cannot normalize zero-norm rows: %s"
        stop(sprintf(fmt, paste(utils::head(zero, 10L), collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}

.aa_unit_rows <- function(X, row_norms = NULL, eps = .Machine$double.eps) {
    X <- as.matrix(X)
    if (is.null(row_norms))
        row_norms <- sqrt(rowSums(X * X))
    safe <- ifelse(row_norms > eps, row_norms, 1)
    X / safe
}

.aa_align_rows <- function(Y, X) {
    sign <- ifelse(rowSums(X * Y) < 0, -1, 1)
    Y * sign
}

.aa_directional_hemisphere <- function(X_gen, method = c("pca", "none")) {
    method <- match.arg(method)
    if (method == "none")
        return(list(X = X_gen, direction = NULL))

    pc <- stats::prcomp(X_gen, center = FALSE, scale. = FALSE, rank. = 1L)
    direction <- pc[["rotation"]]  # (p x 1) unit vector normal to the PCA hyperplane
    direction <- direction / norm(direction, "F")
    projection <- as.numeric(X_gen %*% direction)
    sign <- ifelse(projection < 0, -1, 1)
    list(X = X_gen * sign, direction = as.vector(direction))
}

.aa_directional_init_vars <- function(X_flip, X_gen, K, init, init_args, eps) {
    if (is.character(init) && length(init) == 1L && identical(init, "random"))
        return(.aa_directional_random_init(X_flip, K, eps = eps))

    if (is.matrix(init) || inherits(init, "data.frame")) {
        init <- as.matrix(init)
        if (!identical(dim(init), c(K, ncol(X_flip))))
            stop("`init` matrix must have dimension `K x ncol(data)`.", call. = FALSE)
        .aa_check_no_zero_rows(init)
        B <- fit_simplex(X_flip, .aa_unit_rows(init), eps = eps)
        A <- B %*% X_flip
        S <- fit_simplex(A, X_gen, eps = eps)
        rownames(A) <- rownames(init)
        rownames(B) <- rownames(A)
        colnames(S) <- rownames(A)
        return(list(A0 = A, B = B, S = S))
    }

    init_call <- c(list(X = X_flip, K = K), init_args)
    if (is.character(init)) {
        init_call[["method"]] <- init
        init_vars <- do.call(aa_init, init_call)
    } else if (is.function(init)) {
        init_vars <- do.call(init, init_call)
    } else {
        stop("`init` must be 'random', an aa_init method, a function, or a matrix.",
             call. = FALSE)
    }
    B <- init_vars[["B"]]
    B <- proj_l1(as.matrix(B), eps = eps)
    A <- B %*% X_flip
    S <- fit_simplex(A, X_gen, eps = eps)
    list(A0 = A, B = B, S = S)
}

.aa_fit_directional <- function(X_loss,
                                X_flip,
                                A0,
                                B,
                                S,
                                max_iter,
                                tol,
                                tol_r2,
                                eps,
                                verbose,
                                step_size,
                                max_iter_optimizer,
                                step_shrinkage,
                                max_no_update) {
    # Nomeclature comparison with Olsen et al. (2022)
    #   X_loss: paper X^T, retaining row magnitude as Watson precision.
    #   X_flip: paper X_tilde_f^T, generator data flipped to one hemisphere.
    #   Y: paper X_hat^T, sample reconstructions.
    #   A: not defined in the paper = X_tilde C
    #   B: paper C^T, archetype coefficients/generators.
    #   S: paper S^T, sample compositions.
    #   intermediate terms are defined as in the paper just transposed to row orientation.

    A <- A0
    Y <- S %*% A
    xss <- norm(X_loss, "F")^2
    terms <- .aa_directional_terms(X_loss, Y, xss)
    loss <- .aa_new_loss(max_iter + 1L)
    loss <- .aa_directional_update_loss(loss, 1L, terms)
    converged <- FALSE
    no_update <- 0L
    step_S <- step_size
    step_B <- step_size

    if (verbose) message("Starting directional optimization loop...")
    if (max_iter == 0L) {
        i <- 0L
        loss <- as.data.frame(loss)[1L, , drop = FALSE]
        return(list(A0 = A0, A = A, B = B, S = S, i = i, loss = loss, converged = TRUE))
    }

    for (i in seq_len(max_iter)) {
        j <- i + 1L
        accepted_update <- FALSE

        # S-update
        grad_Y <- .aa_directional_grad_Y(X_loss, Y, terms[["z"]], terms[["q"]]) # TV in Eqs. 6-7
        grad_S <- tcrossprod(grad_Y, A)                                         # Equation 6
        grad_S <- grad_S - rowSums(grad_S * S)                                  # pseudo-grad
        for (k in seq_len(max_iter_optimizer)) {
            S_new <- proj_l1(S - step_S * grad_S, eps = eps)
            Y_new <- S_new %*% A
            terms_new <- .aa_directional_terms(X_loss, Y_new, xss)
            if (terms_new[["rss"]] < terms[["rss"]]) {
                S <- S_new
                Y <- Y_new
                terms <- terms_new
                step_S <- step_S / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }

        # B-update
        if (accepted_update) # use latest version of Y, terms
            grad_Y <- .aa_directional_grad_Y(X_loss, Y, terms[["z"]], terms[["q"]])
        grad_B <- crossprod(S, tcrossprod(grad_Y, X_flip)) # Eq 7: t(S) %*% grad_Y %*% t(X_flip)
        grad_B <- grad_B - rowSums(grad_B * B)             # pseudo-grad
        for (k in seq_len(max_iter_optimizer)) {
            B_new <- proj_l1(B - step_B * grad_B, eps = eps)
            A_new <- B_new %*% X_flip
            Y_new <- S %*% A_new
            terms_new <- .aa_directional_terms(X_loss, Y_new, xss)
            if (terms_new[["rss"]] < terms[["rss"]]) {
                B <- B_new
                A <- A_new
                Y <- Y_new
                terms <- terms_new
                step_B <- step_B / step_shrinkage
                accepted_update <- TRUE
                break
            }
            step_B <- step_B * step_shrinkage
        }

        if (!accepted_update) {
            no_update <- no_update + 1L
            for (nm in names(loss))
                loss[[nm]][j] <- loss[[nm]][i]
            if (verbose) {
                fmt <- paste(
                    "Iteration %d: directional candidate did not improve loss;",
                    "keeping previous iterate (no-update %d/%d)"
                )
                message(sprintf(fmt, i, no_update, max_no_update))
            }
            if (no_update >= max_no_update) {
                warning(sprintf(
                    "Directional AA stalled: no loss improvement after %d consecutive updates",
                    max_no_update
                ), call. = FALSE)
                break
            }
            next
        }

        no_update <- 0L
        loss <- .aa_directional_update_loss(loss, j, terms)
        converged <- .aa_check_convergence(loss, i, tol, tol_r2, Inf, verbose)
        if (converged) break
    }

    list(
        A0 = A0,
        A = A,
        B = B,
        S = S,
        i = i,
        loss = as.data.frame(loss)[seq_len(i + 1L), , drop = FALSE],
        converged = converged
    )
}

.aa_directional_terms <- function(X_loss, Y, xss = norm(X_loss, "F")^2) {
    # Computes the Watson loss between X_loss and Y (equation 3),
    # along with intermediate terms used for gradient computation  z, q, and v
    z <- rowSums(X_loss * Y)                       # (N x 1) vector of x_i^T xhat_i
    q <- pmax(rowSums(Y * Y), .Machine$double.eps) # (N x 1) vector of xhat_i^T xhat_i
    v <- z / sqrt(q)                               # (N x 1) v = z / ||xhat||
    Lw <- norm(v, "2")^2                           # eq. 4
    list(z = z, q = q, rss = xss - Lw, xss = xss)
}

.aa_directional_grad_Y <- function(X_loss, Y, z, q) {
    # Lw     = -z^2/q=z * a, from eq. 3 where a = z/q, z = XtY, q = YtY
    # da/dY  = (q dz/dY - z dq/dY) / q^2 = (X * q - 2Y z) / q^2 = (X - 2Y a) / q
    # dLw/dY = dLw/dz dz/dY + dLw/da da/dY = a dz/dY + z da/dY = a X - z * (X - 2Y a) / q
    #        = 2Y a^2 - X a
    # In the paper nomeclature: grad_Y = TV matrix in Eqs. 6-7
    alpha <- z / q
    -2 * alpha * (X_loss - Y * alpha)
}

.aa_directional_update_loss <- function(loss, i, terms) {
    loss[["loss"]][i] <- terms[["rss"]]
    loss[["r2"]][i] <- 1 - terms[["rss"]] / terms[["xss"]]
    loss
}

.aa_directional_fit_S <- function(X_loss,
                                  A,
                                  max_iter,
                                  eps,
                                  step_size,
                                  max_iter_optimizer,
                                  step_shrinkage,
                                  xss = norm(X_loss, "F")^2) {
    # Start from an ordinary simplex fit against normalized directions, then
    # refine S with projected-gradient steps on the Watson loss from Eq. 3.
    # This mirrors `predict.archetypes()`, which fits new Euclidean
    # compositions while holding coordinates fixed.
    # TODO: The fit is very crude it exits once the line search fails to find an improving step check if adequate
    X_gen <- .aa_unit_rows(X_loss)
    S <- fit_simplex(.aa_unit_rows(A), X_gen, eps = eps)
    Y <- S %*% A
    terms <- .aa_directional_terms(X_loss, Y, xss)
    step_S <- step_size
    for (i in seq_len(max_iter)) {
        grad_Y <- .aa_directional_grad_Y(X_loss, Y, terms[["z"]], terms[["q"]])
        grad_S <- grad_Y %*% t(A)
        grad_S <- grad_S - rowSums(grad_S * S)
        accepted <- FALSE
        for (k in seq_len(max_iter_optimizer)) {
            S_new <- proj_l1(S - step_S * grad_S, eps = eps)
            Y_new <- S_new %*% A
            terms_new <- .aa_directional_terms(X_loss, Y_new, xss)
            if (terms_new[["rss"]] < terms[["rss"]]) {
                S <- S_new
                Y <- Y_new
                terms <- terms_new
                step_S <- step_S / step_shrinkage
                accepted <- TRUE
                break
            }
            step_S <- step_S * step_shrinkage
        }
        if (!accepted) break
    }
    colnames(S) <- rownames(A)
    rownames(S) <- rownames(X_loss)
    S
}
