# Infix Utilities -------------------------------------------------------------------------------

# Null-coalescing operator: x %||% y returns x when non-NULL, otherwise y.
# Imported from rlang; defined here only as documentation anchor.

# List-patch operator: defaults %|p|% overrides fills in missing keys of
# `defaults` with values from `overrides`, mirroring utils::modifyList but
# with the opposite arg order (first arg wins where present).
#
# defaults[names(overrides)] <- overrides; defaults
#
# @keywords internal
`%|p|%` <- function(defaults, overrides) {
    defaults[names(overrides)] <- overrides
    defaults
}

# Scalar / Vector Predicates -------------------------------------------------------------------

# A matrix object from base R or the Matrix package.
is_matrix <- function(x) is.matrix(x) || inherits(x, "Matrix")

# A tabular data object accepted at package API boundaries.
is_tabular <- function(x) is_matrix(x) || inherits(x, "data.frame")

# A single TRUE or FALSE value.
is_logical <- function(x) isTRUE(x) || isFALSE(x)

# A single finite numeric value.
is_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

# A single finite numeric value strictly greater than zero.
is_positive <- function(x) is_number(x) && x > 0

# A single finite numeric value greater than or equal to zero.
is_non_negative <- function(x) is_number(x) && x >= 0

# A single positive whole number (length-1, finite, integer-valued, >= 1).
is_count <- function(x, start_from = 1L) is_number(x) && x == as.integer(x) && x >= start_from

# A numeric vector or matrix with no non-finite values (NA, NaN, Inf, -Inf).
is_all_finite <- function(x) is.numeric(x) && all(is.finite(x))

# Numeric vector or matrix values that are finite, integer-valued, and >= start_from.
is_all_count <- function(x, start_from = 1L) {
    is_all_finite(x) && all(x == as.integer(x)) && all(x >= start_from)
}

is_all_positive <- function(x) is_all_finite(x) && all(x > 0)

is_all_non_negative <- function(x) is_all_finite(x) && all(x >= 0)

is_row_stochastic <- function(x, tol = sqrt(.Machine$double.eps)) {
    Sx <- tryCatch(rowSums(x), error = function(e) NULL)
    if (is.null(Sx)) {
        return(FALSE)
    }
    vals <- if (inherits(x, "Matrix")) methods::slot(x, "x") else as.numeric(x)
    non_negative <- all(is.finite(vals)) && all(vals >= -tol)
    non_negative && isTRUE(all.equal(Sx, rep(1, length(Sx)), check.attributes = FALSE, tolerance = tol))
}

# A single character string.
is_single_string <- function(x) is.character(x) && length(x) == 1L

# A single non-empty character string.
is_non_empty_string <- function(x) is_single_string(x) && nzchar(x)

# Weighting Function for Robust Archetypal Analysis ---------------------------------------------

.aa_robust_label <- function(robust) {
    if (identical(robust, FALSE)) {
        return(NULL)
    }
    if (identical(robust, TRUE)) {
        return("psi.bisquare")
    }
    if (is_non_empty_string(robust)) {
        return(robust)
    }
    if (is.function(robust)) {
        return("function")
    }
    NULL
}

.aa_resolve_robust <- function(robust, robust_args = list()) {
    if (identical(robust, FALSE)) {
        return(NULL)
    }

    if (!is.list(robust_args)) {
        stop("`robust_args` must be a list.", call. = FALSE)
    }
    if (any(names(robust_args) %in% c("u", "deriv"))) {
        stop("`robust_args` cannot include `u` or `deriv`.", call. = FALSE)
    }

    if (identical(robust, TRUE)) {
        robust <- "psi.bisquare"
    }

    error_msg <- "`robust` must be FALSE, TRUE, a MASS psi function name, or a function."
    if (is_non_empty_string(robust)) {
        .aa_require_namespace_for("MASS", "character `robust` psi functions")
        psi_fun <- tryCatch(
            get(robust, envir = asNamespace("MASS"), mode = "function", inherits = FALSE),
            error = function(e) NULL
        )
        if (!is.function(psi_fun)) {
            stop(error_msg, call. = FALSE)
        }
        return(list(fun = psi_fun, label = robust))
    }

    if (is.function(robust)) {
        return(list(fun = robust, label = "function"))
    }
    stop(error_msg, call. = FALSE)
}

.aa_weight_fun <- function(robust, robust_args = list()) {
    spec <- .aa_resolve_robust(robust, robust_args)
    if (is.null(spec)) {
        return(NULL)
    }

    function(row_rss) {
        row_resid <- sqrt(pmax(row_rss, 0))
        scale <- stats::mad(row_resid, center = 0)
        if (!is.finite(scale) || scale <= .Machine$double.eps) {
            return(rep(1, length(row_rss)))
        }

        weights <- do.call(
            spec[["fun"]],
            c(list(u = row_resid / scale, deriv = 0), robust_args)
        )
        .aa_check_row_weights(weights, length(row_rss))
    }
}

.aa_weight_rows <- function(X, row_weights = NULL) {
    if (is.null(row_weights)) {
        return(X)
    }
    X * row_weights
}

.aa_check_row_weights <- function(row_weights, n) {
    if (is.null(row_weights)) {
        return(NULL)
    }
    if (length(row_weights) != n) {
        stop("`row_weights` must have one value per row of `X`.", call. = FALSE)
    }
    if (!is_all_non_negative(row_weights)) {
        stop("`row_weights` must be finite and non-negative.", call. = FALSE)
    }
    row_weights
}


.aa_check_sample_weights <- function(weights, n, allow_null = TRUE) {
    if (is.null(weights)) {
        if (allow_null) {
            return(NULL)
        }
        stop("weights must not be NULL", call. = FALSE)
    }
    if (length(weights) != n) {
        fmt <- "Number of weights (%d) must equal number of rows in data (%d)"
        stop(sprintf(fmt, length(weights), n), call. = FALSE)
    }
    if (!is_all_non_negative(weights)) {
        stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
    if (!any(weights > 0)) {
        stop("at least one weight must be positive.", call. = FALSE)
    }
    if (any(weights == 0)) {
        warning("Some sample weights are zero.", call. = FALSE)
    }
    weights
}

# Mathematical Subroutines -----------------------------------------------------

# Compute squared Euclidean distance of each sample from center
.aa_dist2 <- function(X, center = FALSE) {
    x2 <- Matrix::rowSums(X * X)
    if (isFALSE(center)) {
        return(x2)
    }

    x_mean <- if (isTRUE(center)) colMeans(X) else as.numeric(center)
    stopifnot("center must match columns in X" = length(x_mean) == ncol(X))
    centered_x2 <- x2 - 2 * as.numeric(X %*% x_mean) + as.numeric(x_mean %*% x_mean)
    pmax(centered_x2, 0)
}


# Compute pairwise squared distances between rows of X and Y
.aa_pdist2 <- function(X, Y = NULL) {
    if (is.null(Y)) {
        return(as.matrix(stats::dist(X))^2)
    }

    # both x and y must be in the same dimensional space
    stopifnot(ncol(X) == ncol(Y))

    # squared distances from cosine law ||x||^2 + ||y||^2 − 2<x,y>
    D2 <- outer(rowSums(X * X), rowSums(Y * Y), "+") - 2 * tcrossprod(X, Y)
    pmax(D2, 0) # ensure non-negative distances
}

# Compute pairwise Manhattan distances between rows of X and Y
# TODO: consider using proxy::dist() for more efficient implementations
.aa_pdist_manhattan <- function(X, Y = NULL) {
    if (is.null(Y)) {
        return(as.matrix(stats::dist(X, method = "manhattan")))
    }

    stopifnot(ncol(X) == ncol(Y))
    D <- matrix(0, nrow = nrow(X), ncol = nrow(Y))
    for (j in seq_len(ncol(X))) {
        D <- D + abs(outer(as.vector(X[, j]), as.vector(Y[, j]), "-"))
    }
    D
}

# Centered log-ratio transform after row closure, with zero replacement.
.aa_clr <- function(X, zero_replace = sqrt(.Machine$double.eps)) {
    totals <- rowSums(X)
    if (any(totals <= 0)) {
        stop("Cannot apply clr transform to rows with non-positive total", call. = FALSE)
    }
    closed <- X / totals
    closed <- pmax(closed, zero_replace)
    log_closed <- log(closed)
    log_closed - rowMeans(log_closed)
}

# Otsu's method for 1-D data: finds the threshold maximising between-class
# variance.  Returns NA when the input is empty or has a single unique value.
.aa_otsu_threshold <- function(x, n_bins = 256L) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    rng <- range(x)
    if (rng[1L] == rng[2L]) {
        return(NA_real_)
    }
    breaks <- seq(rng[1L], rng[2L], length.out = n_bins + 1L)
    h <- graphics::hist(x, breaks = breaks, plot = FALSE)
    p <- h$counts / sum(h$counts)
    mids <- h$mids
    omega <- cumsum(p)
    mu <- cumsum(p * mids)
    mu_T <- mu[length(mu)]
    valid <- omega > 0 & omega < 1
    sigma_b2 <- rep(-Inf, length(mids))
    sigma_b2[valid] <- (mu_T * omega[valid] - mu[valid])^2 /
        (omega[valid] * (1 - omega[valid]))
    mids[which.max(sigma_b2)]
}

# Clustering efficiency between original matrix X and reconstructed matrix Y
# (based on the cluster means)
.aa_effic <- function(X, Y) {
    # X, Y: data matrices (N×M), rows = samples, cols = features
    Sx <- stats::cov(X) # (M×M)
    Sy <- stats::cov(Y)

    # Try using Cholesky to invert PD matrix
    cx <- tryCatch(chol(Sx), error = function(e) NULL)
    if (is.null(cx)) {
        warning("cov(X) is singular; using Moore-Penrose pseudo-inverse for efficiency.",
            call. = FALSE
        )
        .aa_require_namespace_for("MASS", "singular-covariance efficiency fallback")
        return(sum(diag(MASS::ginv(Sx) %*% Sy)))
    }

    # since trace(AB) = sum(A * t(B)) = sum(A * B) when symmetric
    sum(Sy * chol2inv(cx))
}

# Truncated symmetric eigendecomposition: top-k eigenvalues and eigenvectors of
# S. Uses RSpectra or irlba when available, otherwise falls back to base eigen().
# Returns list(d = eigenvalues, V = eigenvectors), with near-zero values pruned.
.aa_sym_eigen <- function(S, k) {
    if (requireNamespace("RSpectra", quietly = TRUE)) {
        ev <- RSpectra::eigs_sym(S, k = k, which = "LM")
        d <- ev$values
        V <- ev$vectors
    } else if (requireNamespace("irlba", quietly = TRUE)) {
        ev <- irlba::partial_eigen(S, n = k)
        d <- ev$values
        V <- ev$vectors
    } else {
        ev <- eigen(S, symmetric = TRUE)
        d <- ev$values[seq_len(k)]
        V <- ev$vectors[, seq_len(k), drop = FALSE]
    }
    pos <- d > .Machine$double.eps * max(abs(d))
    list(d = d[pos], V = V[, pos, drop = FALSE])
}

# Kernel PCA from a Gram matrix G (N x N) using generator weights H (K x N).
# Returns a list with $data (N x k training-sample projections) and
# $archetypes (K x k archetype projections).
.aa_kernel_kpca <- function(G, H, k = 2L) {
    cm <- colMeans(G)
    gm <- mean(G)
    Gc <- sweep(G - cm, 2L, cm, "-") + gm
    ev <- .aa_sym_eigen(Gc, k)
    dsq <- sqrt(ev$d)
    V <- ev$V
    k <- length(dsq)
    cnames <- paste0("KPCA", seq_len(k))
    X_proj <- sweep(V, 2L, dsq, "*")
    rownames(X_proj) <- rownames(G)
    colnames(X_proj) <- cnames
    KH <- H %*% G
    KH <- sweep(KH - rowMeans(KH), 2L, cm, "-") + gm
    A_proj <- sweep(KH %*% V, 2L, dsq, "/")
    rownames(A_proj) <- rownames(H)
    colnames(A_proj) <- cnames
    list(data = X_proj, archetypes = A_proj)
}

# Archetypes Fitting Subroutines -----------------------------------------------------

.aa_check_fit_controls <- function(ctx, n = nrow(ctx[["x"]])) {
    if (!is_count(ctx[["max_iter"]], start_from = 0L)) {
        stop("`max_iter` must be a non-negative integer.", call. = FALSE)
    }
    if (!is_positive(ctx[["tol"]])) {
        stop("`tol` must be positive.", call. = FALSE)
    }
    if (!(is_number(ctx[["tol_r2"]]) && ctx[["tol_r2"]] >= 0 && ctx[["tol_r2"]] <= 1)) {
        stop("`tol_r2` must be between 0 and 1.", call. = FALSE)
    }
    if (!(is_count(ctx[["K"]]) && ctx[["K"]] <= n)) {
        stop(
            paste(
                "`K` must be a positive integer less than or equal to",
                "the number of samples."
            ),
            call. = FALSE
        )
    }
    if (!is_non_negative(ctx[["eps"]])) {
        stop("`eps` must be non-negative.", call. = FALSE)
    }
    if (!is.list(ctx[["robust_args"]])) {
        stop("`robust_args` must be a list.", call. = FALSE)
    }
    .aa_resolve_robust(ctx[["robust"]], ctx[["robust_args"]])
    invisible(TRUE)
}

.aa_require_namespace_for <- function(pkg, feature) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        err <- simpleError(
            sprintf("Package `%s` is required for %s. Install it to use this path.", pkg, feature)
        )
        err[["pkg"]] <- pkg
        err[["feature"]] <- feature
        class(err) <- c("aa_missing_namespace_error", class(err))
        stop(err)
    }
    invisible(TRUE)
}

.aa_check_scale <- function(scale, p) {
    if (is_logical(scale)) {
        return(invisible(TRUE))
    }

    if (is.numeric(scale) && is.null(dim(scale))) {
        if (length(scale) != p) {
            stop("Vector `scale` must have one value per feature.", call. = FALSE)
        }
        if (!is_all_positive(scale)) {
            stop("Vector `scale` values must be finite and positive.", call. = FALSE)
        }
        return(invisible(TRUE))
    }

    is_scale_matrix <- is_matrix(scale)
    if (!is_scale_matrix) {
        stop(
            paste(
                "`scale` must be TRUE, FALSE, a numeric vector,",
                "or a dense or sparse matrix."
            ),
            call. = FALSE
        )
    }
    if (!identical(dim(scale), c(p, p))) {
        stop("Matrix `scale` must have one row and column per feature.", call. = FALSE)
    }

    values <- if (isS4(scale) && "x" %in% methods::slotNames(scale)) {
        scale@x
    } else {
        as.vector(scale)
    }
    if (!is_all_finite(values)) {
        stop("Matrix `scale` contains missing or non-finite values.", call. = FALSE)
    }
    is_symmetric <- if (inherits(scale, "Matrix")) {
        Matrix::isSymmetric(scale)
    } else {
        isSymmetric(scale)
    }
    if (!is_symmetric) {
        stop("Matrix `scale` must be symmetric.", call. = FALSE)
    }
    is_pd <- !inherits(try(chol(as.matrix(scale)), silent = TRUE), "try-error")
    if (!is_pd) {
        stop("Matrix `scale` must be positive definite.", call. = FALSE)
    }
    invisible(TRUE)
}


# Check if number of archetypes K corresponds to edge cases (1 or number of samples)
.aa_checks_edge_cases <- function(data, K, verbose = FALSE, M = NULL) {
    out <- NULL
    if (K == nrow(data)) { # X = A
        if (verbose) {
            message("K equals number of samples, returning identity archetypes")
        }
        out <- .aa_identity_fit(data, M = M)
    } else if (K == 1L) { # Archetype = mean of X
        if (verbose) message("K equals 1, returning mean archetype")
        out <- .aa_mean_fit(data, M = M)
    }
    out
}

# Edge case K == N: each sample is its own archetype
.aa_identity_fit <- function(X, M = NULL, call = NULL) {
    A <- X
    rownames(A) <- paste0("A", seq_len(nrow(X))) # remove row names for consistency
    identity <- if (inherits(X, "sparseMatrix")) Matrix::Diagonal(nrow(X)) else diag(nrow(X))
    B <- S <- identity
    dimnames(B) <- list(rownames(A), rownames(X))
    dimnames(S) <- list(rownames(X), rownames(A))
    loss <- data.frame(loss = 0, r2 = 1)

    archetypes(
        call         = NULL,
        data         = X,
        weights      = attr(X, "weights"),
        init         = A,
        A            = A,
        coefficients = B,
        compositions = S,
        loss         = loss,
        converged    = TRUE,
        feature_map  = .aa_identity_feature_map(A)
    )
}

# Edge case K == 1: single archetype at mean of X
.aa_mean_fit <- function(X, M = NULL) {
    n_obs <- if (is.null(M)) rep(nrow(X), ncol(X)) else colSums(M)
    x_sum <- colSums(X)
    x2 <- colSums(X * X)
    x_mean <- x_sum / n_obs
    x_mean[!is.finite(x_mean)] <- 0
    A <- matrix(
        x_mean,
        nrow = 1L,
        ncol = ncol(X),
        dimnames = list("A1", colnames(X))
    )
    B <- matrix(
        1 / nrow(X),
        nrow = 1L,
        ncol = nrow(X),
        dimnames = list("A1", rownames(X))
    )
    S <- matrix(
        1,
        nrow = nrow(X),
        ncol = 1L,
        dimnames = list(rownames(X), "A1")
    )

    xss <- sum(x2)
    rss <- sum(pmax(x2 - ifelse(n_obs > 0, x_sum * x_sum / n_obs, 0), 0))
    loss <- data.frame(loss = rss, r2 = if (xss > 0) 1 - rss / xss else 1)

    archetypes(
        call         = NULL,
        data         = X,
        weights      = attr(X, "weights"),
        init         = A,
        A            = A,
        coefficients = B,
        compositions = S,
        loss         = loss,
        converged    = TRUE,
        feature_map  = .aa_identity_feature_map(A)
    )
}

.aa_col_sds <- function(X) {
    if (!inherits(X, "sparseMatrix")) {
        return(matrixStats::colSds(X))
    }

    n <- nrow(X)
    x_mean <- colMeans(X)
    x2 <- colSums(X * X)
    x_var <- pmax((x2 - n * x_mean * x_mean) / max(n - 1L, 1L), 0)
    out <- sqrt(x_var)
    names(out) <- colnames(X)
    out
}

# Remove features (columns of X) with low variance
.aa_filter_low_variance <- function(X, sd_threshold) {
    sd_vals <- attr(X, "scaled:scale") %||% .aa_col_sds(X)
    mask <- sd_vals >= sd_threshold
    M <- sum(mask)
    attr(X, "scaled:scale") <- sd_vals
    if (M < ncol(X)) {
        # Throw Warning
        fmt <- "The following %d features are filtered out due to low variance: %s"
        dropped_features <- colnames(X)[!mask] %||% which(!mask)
        dropped_features <- paste(dropped_features, collapse = ", ")
        warning(sprintf(fmt, ncol(X) - M, dropped_features), call. = FALSE)

        # Filter X
        x_attrs <- attributes(X)
        x_attrs[["mask"]] <- mask
        x_attrs[["scaled:scale"]] <- sd_vals[mask]
        X <- X[, mask, drop = FALSE]
        if (inherits(X, "sparseMatrix")) {
            attr(X, "mask") <- mask
            attr(X, "scaled:scale") <- sd_vals[mask]
        } else {
            x_attrs[["dim"]] <- dim(X)
            x_attrs[["dimnames"]] <- dimnames(X)
            attributes(X) <- x_attrs
        }
    }
    X
}

.aa_pack_lower_tri <- function(x) {
    x <- as.matrix(x)
    if (nrow(x) != ncol(x)) {
        stop("`x` must be a square matrix.", call. = FALSE)
    }
    c(nrow(x), x[lower.tri(x, diag = TRUE)])
}

.aa_unpack_lower_tri <- function(x) {
    p <- as.integer(x[[1L]])
    expected <- p * (p + 1L) / 2L
    if (!is.finite(p) || p < 0L || length(x) != expected + 1L) {
        stop("Packed lower-triangular matrix has invalid dimensions.", call. = FALSE)
    }

    out <- matrix(0, nrow = p, ncol = p)
    out[lower.tri(out, diag = TRUE)] <- x[-1L]
    out
}

.aa_preprocess_missing <- function(data, sd_threshold, verbose, scale = TRUE) {
    if (is_matrix(scale)) {
        stop("matrix `scale` is not supported with `missing = TRUE`.", call. = FALSE)
    }

    requested_scale_true <- isTRUE(scale)
    scale_mode <- if (identical(scale, FALSE)) "none" else "vector"
    sparse <- inherits(data, "sparseMatrix")

    if (sparse) {
        data <- Matrix::drop0(data)
        M <- Matrix::drop0(data != 0)
        stats <- .aa_observed_col_stats(data, M)
        original_center <- stats[["mean"]]
        X <- data
    } else {
        data <- as.matrix(data)
        M <- !is.na(data)
        X <- data
        X[!M] <- 0
        stats <- .aa_observed_col_stats(X, M)
        original_center <- stats[["mean"]]
    }
    names(original_center) <- colnames(data)

    attr(X, "scaled:scale") <- stats[["sd"]]
    if (isTRUE(scale)) {
        scale <- stats[["sd"]]
        safe_scale <- ifelse(scale > 0, scale, 1)
        if (sparse) {
            X <- .aa_scale_sparse_observed(
                X,
                center = stats[["mean"]],
                scale = safe_scale
            )
        } else {
            X <- sweep(data, 2L, stats[["mean"]], "-")
            X[!M] <- 0
            X <- sweep(X, 2L, safe_scale, "/")
        }
        attr(X, "scaled:center") <- stats[["mean"]]
        attr(X, "scaled:scale") <- stats[["sd"]]
    }

    if (!sparse && mean(M) < 0.10) {
        M <- Matrix::Matrix(M, sparse = TRUE)
        X <- Matrix::drop0(Matrix::Matrix(X, sparse = TRUE))
        attr(X, "scaled:scale") <- stats[["sd"]]
        if (isTRUE(scale)) {
            attr(X, "scaled:center") <- stats[["mean"]]
        }
    }

    X <- .aa_filter_low_variance(X, sd_threshold)
    mask <- attr(X, "mask")
    if (!is.null(mask)) {
        M <- M[, mask, drop = FALSE]
    }

    if (!requested_scale_true && is.numeric(scale) && is.null(dim(scale))) {
        if (!is.null(mask)) {
            scale <- scale[mask]
        }
        scale_factor <- ifelse(scale > 0, scale, 1)
        x_attrs <- attributes(X)
        X <- if (inherits(X, "sparseMatrix")) {
            Matrix::colScale(X, 1 / scale_factor)
        } else {
            sweep(X, 2L, scale_factor, "/")
        }
        attributes(X) <- utils::modifyList(attributes(X), x_attrs)
        attr(X, "scale:factor") <- scale_factor
    } else if (requested_scale_true) {
        if (!is.null(mask)) {
            scale <- scale[mask]
        }
        attr(X, "scale:factor") <- ifelse(scale > 0, scale, 1)
    }

    attr(X, "scale:mode") <- scale_mode
    attr(X, "restore:center") <- original_center
    attr(X, "missing") <- TRUE
    if (inherits(M, "sparseMatrix")) {
        M <- Matrix::drop0(M)
    }

    list(X = X, M = M)
}

.aa_observed_col_stats <- function(X, M) {
    n <- colSums(M)
    x_sum <- colSums(X)
    x_mean <- x_sum / n
    x_mean[!is.finite(x_mean)] <- 0
    x2 <- colSums(X * X)
    x_var <- (x2 - n * x_mean * x_mean) / pmax(n - 1L, 1L)
    x_var[n <= 1L] <- 0
    x_sd <- sqrt(pmax(x_var, 0))
    names(x_mean) <- names(x_sd) <- colnames(X)
    list(mean = x_mean, sd = x_sd, n = n)
}

.aa_scale_sparse_observed <- function(X, center, scale) {
    entries <- Matrix::summary(X)
    if (length(entries[["i"]]) == 0L) {
        return(X)
    }
    vals <- (entries[["x"]] - center[entries[["j"]]]) / scale[entries[["j"]]]
    Matrix::drop0(Matrix::sparseMatrix(
        i = entries[["i"]],
        j = entries[["j"]],
        x = vals,
        dims = dim(X),
        dimnames = dimnames(X)
    ))
}

.aa_fit_space_col_sds <- function(X) {
    scale_mode <- attr(X, "scale:mode")
    data_scale <- attr(X, "scaled:scale")
    if (identical(scale_mode, "none") && !is.null(data_scale)) {
        return(data_scale)
    }

    if (identical(scale_mode, "vector")) {
        scale_factor <- attr(X, "scale:factor")
        if (!is.null(data_scale) && !is.null(scale_factor)) {
            return(data_scale / scale_factor)
        }
    }

    .aa_col_sds(X)
}

.aa_auto_bigM <- function(X, multiplier = 200) {
    feature_sd <- .aa_fit_space_col_sds(X)
    feature_scale <- sqrt(mean(feature_sd^2))
    if (!is.finite(feature_scale) || feature_scale <= .Machine$double.eps) {
        feature_scale <- 1
    }
    multiplier * feature_scale * sqrt(2 / max(ncol(X), 1L))
}

# Common subroutine to preprocess `init` matrix of archetype coordinates (A);
# maps `init` from original data space to preprocessed space of X.
.aa_preprocess_init <- function(init, X) {
    # If `init` is not a matrix, return it as is (it will be processed by the init function)
    if (!is_tabular(init)) {
        return(init)
    }

    iM <- attr(X, "bigM")
    n_features <- ncol(X) - !is.null(iM)
    mask <- attr(X, "mask")

    if (ncol(init) != n_features) {
        can_subset <- !is.null(mask) && ncol(init) == length(mask)
        if (!can_subset) {
            fmt <- "ncol(init) = %d does not match number of data features (%d)"
            stop(sprintf(fmt, ncol(init), n_features))
        }
    }

    # Scale `init` to match X
    init <- as.matrix(init)
    x_center <- attr(X, "scaled:center")
    if (!is.null(x_center)) {
        init <- sweep(init, 2L, x_center, "-")
    }

    # Filter `init` to match filtered features in X
    if (!is.null(mask) && ncol(init) == length(mask)) {
        init <- init[, mask, drop = FALSE]
    }

    if (identical(attr(X, "scale:mode"), "matrix")) {
        init <- init %*% .aa_unpack_lower_tri(attr(X, "scale:factor"))
    }
    if (identical(attr(X, "scale:mode"), "vector")) {
        init <- sweep(init, 2L, attr(X, "scale:factor"), "/")
    }

    # Add bigM column to init
    if (!is.null(iM)) {
        init <- cbind(
            matrix(attr(X, "bigM.value"), nrow = nrow(init), ncol = 1L),
            init
        )
        colnames(init)[iM] <- "bigM"
    }

    init
}

# Check that archetype names are valid and return them;
# if missing, generate default names A1, A2, ..., AK
.aa_init_names <- function(A) {
    nm <- rownames(A)
    if (is.null(nm)) {
        return(paste0("A", seq_len(nrow(A))))
    }

    if (any(is.na(nm))) {
        stop("Archetype names must not be missing.", call. = FALSE)
    }
    if (!all(nzchar(nm))) {
        stop("Archetype names must not be empty.", call. = FALSE)
    }
    if (anyDuplicated(nm)) {
        stop("Archetype names must be unique.", call. = FALSE)
    }
    nm
}

# Helper function to .aa_init_vars when `init` is a matrix of archetype coordinates (A)
.aa_matrix_init <- function(X, K, init, eps, L, delta = 0, tol = 1e-6) {
    # Check dimensionality of `init` matrix
    if (nrow(init) != K) {
        fmt <- "nrow(init) = %d does not match K (%d)"
        stop(sprintf(fmt, nrow(init), K))
    }
    if (ncol(init) != ncol(X)) {
        fmt <- "ncol(init) = %d does not match preprocessed data features (%d)"
        stop(sprintf(fmt, ncol(init), ncol(X)))
    }

    # Check that `init` coordinates are within the data hull (or within `delta` slack)
    nm <- .aa_init_names(init)
    a_lo <- max(1 - delta, ifelse(eps > 0, eps, 1e-8))
    a_hi <- 1 + delta
    B <- tryCatch(
        .aa_fit_qp(
            A = X,
            X = init,
            eps = eps,
            project = if (delta == 0) proj_l1 else NULL,
            row_sum_bounds = c(a_lo, a_hi),
            feature = "matrix-valued `init`"
        ),
        aa_missing_namespace_error = function(e) {
            if (identical(e[["pkg"]], "quadprog")) {
                msg <- paste(
                    "Matrix-valued `init` requires the `quadprog` package.",
                    "Install `quadprog` to use custom archetype coordinates."
                )
                stop(msg, call. = FALSE)
            }
            stop(e)
        }
    )
    A <- B %*% X
    err <- norm(A - init, type = "F")
    if (any(err > tol)) {
        ix <- which(err > tol)
        fmt <- paste(
            "Initial archetype coordinates outside the allowed data hull were",
            "projected; affected rows: %s"
        )
        warning(sprintf(fmt, paste(utils::head(ix, 10L), collapse = ", ")), call. = FALSE)
    }

    rownames(A) <- rownames(B) <- nm
    colnames(B) <- rownames(X)
    list(
        A = A,
        B = B,
        S = .aa_init_S(X, A, eps = eps),
        loss = list(loss = rep(NA_real_, L), r2 = rep(NA_real_, L))
    )
}

.aa_trace_row_rss <- function(row_xss, S, XAt, AAt) {
    pmax(row_xss - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt)), 0)
}

# Check convergence based on relative loss improvement and R2 threshold
.aa_check_convergence <- function(loss, i, tol, tol_r2, verbose) {
    j <- i + 1L
    converged <- with(
        loss,
        abs(loss[j] - loss[i]) < tol * loss[i] || r2[j] > tol_r2
    )
    if (verbose && i %% 10 == 0) {
        fmt <- "Iteration %d: loss = %.4f, R2 = %.3f"
        message(sprintf(fmt, i, loss[["loss"]][j], loss[["r2"]][j]))
    }
    converged
}
