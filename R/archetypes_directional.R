#' Directional Archetypal Analysis
#'
#' Fits directional archetypal analysis (DAA) via projected gradient descent,
#' modelling observations by their orientation rather than Euclidean magnitude.
#'
#' @param x dense numeric data matrix (rows = samples, columns = dimensions).
#' @param K number of archetypes.
#' @param init initialization method; see [run_aa()] for available options.
#' @param init_args list of additional arguments for the initialization function.
#' @param weights optional non-negative sample weights (default: NULL)
#' @param max_iter maximum number of outer iterations (default: 100)
#' @param tol convergence tolerance based on directional residual loss (default: 1e-6)
#' @param tol_r2 convergence tolerance based on directional R\eqn{^2} (default: 0.9999)
#' @param max_kappa accepted for consistency with other fitters; directional AA
#'   does not solve by matrix inversion, so condition numbers are not computed
#'   (default: 1000)
#' @param eps small positive number for numerical stability (default: 1e-8)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param hemisphere hemisphere handling. `"pca"` flips generator rows onto the
#'   dominant PCA hemisphere; `"none"` leaves signs unchanged (default: `"pca"`)
#' @param precision precision weighting. `"row_norm"` keeps row magnitudes in
#'   the Watson loss, matching Olsen et al.; `"unit"` gives every row equal
#'   precision (default: `"row_norm"`)
#' @param step_size initial projected-gradient step size (default: 1.0)
#' @param max_iter_optimizer maximum line-search iterations per update (default: 10)
#' @param step_shrinkage factor used to shrink rejected line-search steps (default: 0.5)
#' @param max_no_update maximum consecutive outer iterations with no accepted
#'   line-search update before stopping as stalled (default: 5)
#'
#' @details
#' ## Directional AA
#'
#' Standard AA minimises squared Euclidean reconstruction error, which conflates
#' direction and magnitude. Directional AA row-normalises the input so archetypes
#' are convex combinations of unit vectors, and evaluates reconstructions via a
#' Watson-style angular loss that is polarity-invariant: a sample and its antipodal
#' reflection contribute equally. This makes the method suitable for data where only
#' orientation matters, such as neuroimaging source vectors or unit-sphere embeddings.
#'
#' ## Hemisphere alignment
#'
#' The Watson loss treats a unit vector and its negation as identical. Without
#' alignment, opposing signs can produce phantom "average direction" archetypes
#' near the equator.
#'
#' `hemisphere = "pca"` (default) projects each row onto the first principal
#' component and flips rows with a negative dot product, so all generators lie
#' on the same side. `hemisphere = "none"` leaves signs unchanged; use this
#' when the data are already hemispherically consistent.
#'
#' ## Precision weighting
#'
#' The Watson loss weights each sample's angular residual by the squared norm
#' of its reconstruction. `precision = "row_norm"` (default) retains the
#' original row magnitudes, giving naturally larger observations higher influence
#' matching the formulation in Olsen et al. (2022). `precision = "unit"`
#' normalises every row to unit length before computing the loss, treating all
#' samples as equally reliable.
#'
#' ## Initialization
#'
#' The default `dirichlet` initialization draws the `coefficient` matrix
#' as a Dirichlet(1, ..., 1), which is the simplex equivalent
#' of the uniform distribution. This differs from `random` (the [aa_init()]
#' method that selects K rows of the data uniformly at random), which initializes
#' `B` as a one-hot matrix so that archetypes start at actual data points
#' on the hemisphere.
#'
#' @returns An object of class \code{\link{directional_archetypes}}, extending
#'   \code{\link{archetypes}}.
#'
#' @seealso [run_aa()] for the common entry point and full parameter documentation.
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
archetypes_directional <- function(x,
                                   K,
                                   init = "dirichlet",
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
    .aa_fit_engine(
        call = match.call(),
        x = x,
        K = K,
        method = "directional",
        init = init,
        init_args = init_args,
        weights = weights,
        scale = TRUE,
        robust = FALSE,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        missing = FALSE,
        hemisphere = hemisphere,
        precision = precision,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = max_no_update
    )
}

.aa_directional_block <- function(ctx,
                                  hemisphere = c("pca", "none"),
                                  precision = c("row_norm", "unit"),
                                  step_size = 1.0,
                                  max_iter_optimizer = 10L,
                                  step_shrinkage = 0.5,
                                  max_no_update = 5L) {
    hemisphere <- match.arg(hemisphere)
    precision <- match.arg(precision)

    list(
        check = function(ctx) {
            if (ctx[["missing"]]) {
                stop("`missing = TRUE` is not supported for `method = 'directional'`.",
                    call. = FALSE
                )
            }
            if (ctx[["robust"]]) {
                stop("`robust = TRUE` is not supported for `method = 'directional'`.",
                    call. = FALSE
                )
            }
            if (!isTRUE(ctx[["scale"]])) {
                stop("`scale` is not supported for `method = 'directional'`.", call. = FALSE)
            }
            if (inherits(ctx[["x"]], "sparseMatrix")) {
                stop("Directional AA currently requires a dense numeric matrix.", call. = FALSE)
            }
            X <- as.matrix(ctx[["x"]])
            .aa_check_fit_controls(ctx, n = nrow(X))
            .aa_check_projected_gradient_controls(
                step_size = step_size,
                max_iter_optimizer = max_iter_optimizer,
                step_shrinkage = step_shrinkage,
                max_no_update = max_no_update
            )
            if (!is.null(ctx[["weights"]])) {
                stopifnot(
                    "weights must match rows in data" =
                        length(ctx[["weights"]]) == nrow(X)
                )
                stopifnot("weights contain NA values" = !any(is.na(ctx[["weights"]])))
                stopifnot("weights must be non-negative" = all(ctx[["weights"]] >= 0))
                stopifnot("at least one weight must be positive" = any(ctx[["weights"]] > 0))
            }
            .aa_check_no_zero_rows(X)
        },
        preprocess = function(ctx) {
            X <- as.matrix(ctx[["x"]])
            row_norms <- sqrt(rowSums(X * X))
            X_gen <- .aa_unit_rows(X, row_norms)
            X_loss <- if (precision == "row_norm") X else X_gen
            if (!is.null(ctx[["weights"]])) {
                weights <- ctx[["weights"]] / mean(ctx[["weights"]])
                X_loss <- X_loss * sqrt(weights)
            }
            flip <- .aa_directional_hemisphere(X_gen, method = hemisphere)
            list(
                X = X,
                X_gen = X_gen,
                X_loss = X_loss,
                X_flip = flip[["X"]],
                hemisphere_direction = flip[["direction"]],
                row_norms = row_norms,
                precision = precision
            )
        },
        edge_case = function(ctx, prep) NULL,
        init = function(ctx, prep) {
            .aa_directional_init_vars(
                X_flip = prep[["X_flip"]],
                X_gen = prep[["X_gen"]],
                K = ctx[["K"]],
                init = ctx[["init"]],
                init_args = ctx[["init_args"]],
                eps = ctx[["eps"]]
            )
        },
        fit = function(ctx, prep, init_vars) {
            .aa_fit_directional(
                X_loss = prep[["X_loss"]],
                X_flip = prep[["X_flip"]],
                A0 = init_vars[["A0"]],
                B = init_vars[["B"]],
                S = init_vars[["S"]],
                max_iter = ctx[["max_iter"]],
                tol = ctx[["tol"]],
                tol_r2 = ctx[["tol_r2"]],
                eps = ctx[["eps"]],
                verbose = ctx[["verbose"]],
                step_size = step_size,
                max_iter_optimizer = max_iter_optimizer,
                step_shrinkage = step_shrinkage,
                max_no_update = as.integer(max_no_update)
            )
        },
        final_loss = .aa_final_loss,
        prepare_output = function(ctx, prep, fit) {
            directional_archetypes(
                coordinates = fit[["A"]],
                coefficients = fit[["B"]],
                compositions = fit[["S"]],
                loss = fit[["loss"]],
                converged = fit[["converged"]],
                call = ctx[["call"]],
                data = prep[["X"]],
                init = fit[["A0"]],
                generator_data = prep[["X_flip"]],
                hemisphere_direction = prep[["hemisphere_direction"]],
                row_norms = prep[["row_norms"]],
                precision = prep[["precision"]],
                fit_info = list(method = "directional", precision = prep[["precision"]]) %|p|%
                    list(
                        family = "watson",
                        robust = FALSE,
                        missing = FALSE,
                        delta = 0,
                        init = ctx[["init"]],
                        scaling = "none",
                        sample_weights = !is.null(ctx[["weights"]])
                    )
            )
        }
    )
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
    row_norms <- row_norms %||% sqrt(rowSums(X * X))
    safe <- ifelse(row_norms > eps, row_norms, 1)
    X / safe
}

.aa_align_rows <- function(Y, X) {
    sign <- ifelse(rowSums(X * Y) < 0, -1, 1)
    Y * sign
}

.aa_directional_hemisphere <- function(X_gen, method = c("pca", "none")) {
    method <- match.arg(method)
    if (method == "none") {
        return(list(X = X_gen, direction = NULL))
    }

    pc <- stats::prcomp(X_gen, center = FALSE, scale. = FALSE, rank. = 1L)
    direction <- pc[["rotation"]] # (p x 1) unit vector normal to the PCA hyperplane
    direction <- direction / norm(direction, "F")
    projection <- as.numeric(X_gen %*% direction)
    sign <- ifelse(projection < 0, -1, 1)
    list(X = X_gen * sign, direction = as.vector(direction))
}

.aa_directional_init_vars <- function(X_flip, X_gen, K, init, init_args, eps) {
    if (is.matrix(init) || inherits(init, "data.frame")) {
        init <- as.matrix(init)
        if (!identical(dim(init), c(K, ncol(X_flip)))) {
            stop("`init` matrix must have dimension `K x ncol(data)`.", call. = FALSE)
        }
        .aa_check_no_zero_rows(init)
        B <- fit_simplex(X_flip, .aa_unit_rows(init), eps = eps)
        A <- B %*% X_flip
        S <- fit_simplex(A, X_gen, eps = eps)
        rownames(A) <- rownames(init)
        rownames(B) <- rownames(A)
        colnames(S) <- rownames(A)
        return(list(A0 = A, B = B, S = S))
    }

    init_call <- c(list(X = X_flip, K = K, eps = eps), init_args)
    if (is.character(init)) {
        init_call[["method"]] <- init
        init_vars <- do.call(aa_init, init_call)
    } else if (is.function(init)) {
        init_vars <- do.call(init, init_call)
    } else {
        stop("`init` must be an aa_init method string, a function, or a matrix.",
            call. = FALSE
        )
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
        grad_S <- tcrossprod(grad_Y, A) # Equation 6
        grad_S <- grad_S - rowSums(grad_S * S) # pseudo-grad
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
        if (accepted_update) { # use latest version of Y, terms
            grad_Y <- .aa_directional_grad_Y(X_loss, Y, terms[["z"]], terms[["q"]])
        }
        grad_B <- crossprod(S, tcrossprod(grad_Y, X_flip)) # Eq 7: t(S) %*% grad_Y %*% t(X_flip)
        grad_B <- grad_B - rowSums(grad_B * B) # pseudo-grad
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
            for (nm in names(loss)) {
                loss[[nm]][j] <- loss[[nm]][i]
            }
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
    z <- rowSums(X_loss * Y) # (N x 1) vector of x_i^T xhat_i
    q <- pmax(rowSums(Y * Y), .Machine$double.eps) # (N x 1) vector of xhat_i^T xhat_i
    v <- z / sqrt(q) # (N x 1) v = z / ||xhat||
    Lw <- norm(v, "2")^2 # eq. 4
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
