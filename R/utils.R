# Weighting Function for Robust Archetypal Analysis ---------------------------------------------

# Tukey's bisquare row weights from squared row residual norms.
.aa_bisquare_weights <- function(row_rss, c = 4.685) {
    stopifnot("`c` must be positive" = length(c) == 1L && is.finite(c) && c > 0)
    row_resid <- sqrt(pmax(row_rss, 0))
    scale <- stats::mad(row_resid, center = 0)
    if (!is.finite(scale) || scale <= .Machine$double.eps)
        return(rep(1, length(row_rss)))

    u <- row_resid / scale
    ifelse(u <= c, (1 - (u / c)^2)^2, 0)
}

.aa_trivial_row_weights <- function(row_weights) {
    is.null(row_weights) || all(row_weights == 1)
}

.aa_weight_rows <- function(X, row_weights = NULL) {
    if (.aa_trivial_row_weights(row_weights))
        return(X)
    X * row_weights
}

.aa_check_row_weights <- function(row_weights, n) {
    stopifnot("row_weights must match rows in X" = length(row_weights) == n)
    stopifnot("row_weights contain NA values" = !any(is.na(row_weights)))
    stopifnot("row_weights must be non-negative" = all(row_weights >= 0))
    invisible(TRUE)
}

# Scale data without forcing sparse inputs through dense centering.
scale_safe <- function(X, center = !inherits(X, "sparseMatrix"), scale = TRUE) {
    is_sparse <- inherits(X, "sparseMatrix")
    X_attrs <- attributes(X)

    # Normalize center
    if (is.numeric(center)) {
        if (length(center) == 1L)
            center <- rep(center, ncol(X))
        stopifnot(length(center) == ncol(X))
        if (all(center == 0))
            center <- FALSE
    }

    if (isTRUE(center))
        center <- colMeans(X)

    if (!isFALSE(center)) {
        X <- sweep(as.matrix(X), 2L, center, "-")
        if (is_sparse)
            warning("Centering matrices breaks sparsity; consider using `center = FALSE`", call. = FALSE)
        is_sparse <- FALSE
    }

    if (!isTRUE(scale)) {
        if (!isFALSE(center))
            attr(X, "scaled:center") <- center
        return(X)
    }

    if (is.numeric(scale)) {
        if (length(scale) == 1L)
            scale <- rep(scale, ncol(X))
        stopifnot(length(scale) == ncol(X))
        x_scale <- scale
    } else {
        n <- nrow(X)
        x_mean <- colMeans(X)
        x2 <- colSums(X * X)
        x_var <- pmax((x2 - n * x_mean * x_mean) / max(n - 1L, 1L), 0)
        x_scale <- sqrt(x_var)
    }
    names(x_scale) <- colnames(X)
    safe_scale <- ifelse(x_scale > 0, x_scale, 1)

    X <- if (is_sparse) Matrix::colScale(X, 1 / safe_scale) else sweep(X, 2L, safe_scale, "/")
    if (!is.null(X_attrs[["dimnames"]]))
        dimnames(X) <- X_attrs[["dimnames"]]
    if (!isFALSE(center))
        attr(X, "scaled:center") <- center
    attr(X, "scaled:scale") <- x_scale
    X
}

# Mathematical Subroutines -----------------------------------------------------

# Compute squared Euclidean distance of each sample from center
.dist2 <- function(X, center = FALSE) {
    x2 <- rowSums(X * X)
    if (isFALSE(center))
        return(x2)

    x_mean <- if (isTRUE(center)) colMeans(X) else as.numeric(center)
    stopifnot("center must match columns in X" = length(x_mean) == ncol(X))
    centered_x2 <- x2 - 2*as.numeric(X %*% x_mean) + as.numeric(x_mean %*% x_mean)
    pmax(centered_x2, 0)
}


# Compute pairwise squared distances between columns of X and Y
.pdist2 <- function(X, Y) {
    # both x and y must be in the same dimensional space
    stopifnot(ncol(X) == ncol(Y))

    # squared distances from cosine law ||x||^2 + ||y||^2 − 2<x,y>
    D2 <- outer(rowSums(X * X), rowSums(Y * Y), "+") - 2*tcrossprod(X, Y)
    pmax(D2, 0)     # ensure non-negative distances
}

# Otsu's method for 1-D data: finds the threshold maximising between-class
# variance.  Returns NA when the input is empty or has a single unique value.
.otsu_threshold <- function(x, n_bins = 256L) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    rng <- range(x)
    if (rng[1L] == rng[2L]) return(NA_real_)
    breaks <- seq(rng[1L], rng[2L], length.out = n_bins + 1L)
    h <- graphics::hist(x, breaks = breaks, plot = FALSE)
    p <- h$counts / sum(h$counts)
    mids <- h$mids
    omega <- cumsum(p)
    mu    <- cumsum(p * mids)
    mu_T  <- mu[length(mu)]
    valid <- omega > 0 & omega < 1
    sigma_b2 <- rep(-Inf, length(mids))
    sigma_b2[valid] <- (mu_T * omega[valid] - mu[valid])^2 /
        (omega[valid] * (1 - omega[valid]))
    mids[which.max(sigma_b2)]
}

# Clustering efficiency between original matrix X and reconstructed matrix Y
# (based on the cluster means)
effic <- function(X, Y) {
    # X, Y: data matrices (N×M), rows = samples, cols = features
    N <- nrow(X)
    M <- ncol(X)
    Sx <- stats::cov(X)  # (M×M)

    # Try using Cholesky to invert PD matrix
    cx <- tryCatch(chol(Sx), error = function(e) NULL)
    if (is.null(cx))  # FALLBACK to general solution
        return(sum(diag(solve(Sx, stats::cov(Y)))))


    # Edge Case: # of samples smaller than the # of features
    if (N < M) {
        # avoid computing Sy and use the Mahalanobis approach instead.
        Yc <- sweep(Y, 2L, colMeans(Y), FUN = "-", check.margin = FALSE)
        Z  <- backsolve(cx, t(Yc))
        # remember cov(Y) uses denominator (N−1)
        return(norm(Z, type = "F")^2 / (N - 1))
    }
    # since trace(AB) = sum(A * t(B)) = sum(A * B) when symmetric
    sum(stats::cov(Y) * chol2inv(cx))
}

# Archetypes Fitting Subroutines -----------------------------------------------------

.aa_check_fit_controls <- function(ctx, n = nrow(ctx[["x"]])) {
    stopifnot("max_iter must be a non-negative integer" =
                  ctx[["max_iter"]] == as.integer(ctx[["max_iter"]]) &&
                      ctx[["max_iter"]] >= 0L)
    stopifnot("tol must be positive" = ctx[["tol"]] > 0)
    stopifnot("tol_r2 must be between (0, 1)" =
                  ctx[["tol_r2"]] >= 0 && ctx[["tol_r2"]] <= 1)
    stopifnot("K must be an integer" = ctx[["K"]] == as.integer(ctx[["K"]]))
    stopifnot("K must be an integer greater or equal to 1" = ctx[["K"]] >= 1L)
    stopifnot("K cannot be greater than number of samples" = ctx[["K"]] <= n)
    stopifnot("max_kappa must be >=1" = ctx[["max_kappa"]] >= 1)
    stopifnot("eps must be non-negative" = ctx[["eps"]] >= 0)
    stopifnot("robust must be TRUE or FALSE" =
                  is.logical(ctx[["robust"]]) && length(ctx[["robust"]]) == 1L &&
                      !is.na(ctx[["robust"]]))
    stopifnot("tukey_c must be positive" =
                  length(ctx[["tukey_c"]]) == 1L && is.finite(ctx[["tukey_c"]]) &&
                      ctx[["tukey_c"]] > 0)
    invisible(TRUE)
}

.aa_new_loss <- function(L) {
    list(loss = rep(NA_real_, L),
         r2   = rep(NA_real_, L),
         k_S  = rep(NA_real_, L),
         k_A  = rep(NA_real_, L))
}

.aa_check_scale <- function(scale, p) {
    if (isTRUE(scale) || identical(scale, FALSE))
        return(invisible(TRUE))

    if (is.numeric(scale) && is.null(dim(scale))) {
        stopifnot("vector scale must have one value per feature" = length(scale) == p)
        stopifnot("vector scale contains missing or non-finite values" =
                      !any(is.na(scale)) && all(is.finite(scale)))
        stopifnot("vector scale must be positive" = all(scale > 0))
        return(invisible(TRUE))
    }

    is_scale_matrix <- is.matrix(scale) || inherits(scale, "Matrix")
    stopifnot("scale must be TRUE, FALSE, a numeric vector, or a dense or sparse matrix" =
                  is_scale_matrix)
    stopifnot("matrix scale must have one row and column per feature" =
                  identical(dim(scale), c(p, p)))

    values <- if (isS4(scale) && "x" %in% methods::slotNames(scale)) {
        scale@x
    } else {
        as.vector(scale)
    }
    stopifnot("matrix scale must be numeric" = is.numeric(values))
    stopifnot("matrix scale contains missing or non-finite values" =
                  !any(is.na(values)) && all(is.finite(values)))
    is_symmetric <- if (inherits(scale, "Matrix")) {
        Matrix::isSymmetric(scale)
    } else {
        isSymmetric(scale)
    }
    stopifnot("matrix scale must be symmetric" = is_symmetric)
    is_pd <- !inherits(try(chol(as.matrix(scale)), silent = TRUE), "try-error")
    stopifnot("matrix scale must be positive definite" = is_pd)
    invisible(TRUE)
}


# Check if number of archetypes K corresponds to edge cases (1 or number of samples)
.aa_checks_edge_cases <- function(data, K, verbose = FALSE) {
    out <- NULL
    if (K == nrow(data)) { # X = A
        if (verbose)
            message("K equals number of samples, returning identity archetypes")
        out <- .identity_archetypes(data)
    } else if (K == 1L) { # Archetype = mean of X
        if (verbose) message("K equals 1, returning mean archetype")
        out <- .mean_archetype(data)
    }
    out
}

# Edge case K == N: each sample is its own archetype
.identity_archetypes <- function(X, call = NULL) {
    A <- X
    rownames(A) <- paste0("A", seq_len(nrow(X)))  # remove row names for consistency
    S <- B <- diag(nrow(X))
    colnames(B) <- rownames(S) <- rownames(X)
    rownames(B) <- colnames(S) <- rownames(A)
    loss <- data.frame(loss = 0, r2 = 1, k_S = 1, k_A = kappa(A))

    archetypes(
        call         = NULL,
        data         = X,
        init         = A,
        coordinates  = A,
        coefficients = B,
        compositions = S,
        loss         = loss,
        converged    = TRUE
    )
}

# Edge case K == 1: single archetype at mean of X
.mean_archetype <- function(X) {
    x_mean <- colMeans(X)
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

    xss <- norm(X, type = "F")^2
    rss <- xss - nrow(X) * as.numeric(x_mean %*% x_mean)
    loss <- data.frame(loss = rss, r2 = 1 - rss / xss, k_S = 1, k_A = 1)

    archetypes(
        call         = NULL,
        data         = X,
        init         = A,
        coordinates  = A,
        coefficients = B,
        compositions = S,
        loss         = loss,
        converged    = TRUE
    )
}

# Remove features (columns of X) with low variance
.filter_low_variance <- function(X, sd_threshold) {
    sd_vals <- attr(X, "scaled:scale")
    if (is.null(sd_vals))
        sd_vals <- matrixStats::colSds(X)
    mask <- sd_vals >= sd_threshold
    M <- sum(mask)
    attr(X, "scaled:scale") <- sd_vals
    if (M < ncol(X)) {
        # Throw Warning
        fmt <- "The following %d features are filtered out due to low variance: %s"
        dropped_features <- if (is.null(colnames(X))) which(!mask) else colnames(X)[!mask]
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
            attributes(X) <- x_attrs
        }
    }
    X
}


.aa_preprocess_missing <- function(data, sd_threshold, verbose, scale = TRUE) {
    if (is.matrix(scale) || inherits(scale, "Matrix"))
        stop("matrix `scale` is not supported with `missing = TRUE`.", call. = FALSE)

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
        if (isTRUE(scale))
            attr(X, "scaled:center") <- stats[["mean"]]
    }

    X <- .filter_low_variance(X, sd_threshold)
    mask <- attr(X, "mask")
    if (!is.null(mask))
        M <- M[, mask, drop = FALSE]

    if (!requested_scale_true && is.numeric(scale) && is.null(dim(scale))) {
        if (!is.null(mask))
            scale <- scale[mask]
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
        if (!is.null(mask))
            scale <- scale[mask]
        attr(X, "scale:factor") <- ifelse(scale > 0, scale, 1)
    }

    attr(X, "scale:mode") <- scale_mode
    attr(X, "restore:center") <- original_center
    attr(X, "missing") <- TRUE
    if (inherits(M, "sparseMatrix"))
        M <- Matrix::drop0(M)

    list(X = X, M = M, undo_scale = .aa_undo_scale)
}

.aa_undo_scale <- function(mat, X) {
    stopifnot(ncol(mat) == ncol(X))

    # Remove bigM if present
    iM <- attr(X, "bigM")
    x_names <- colnames(X)
    if (!is.null(iM))
        mat <- mat[, -iM, drop = FALSE]
    if (!is.null(iM))
        x_names <- x_names[-iM]
    mat <- as.matrix(mat)

    scale_mode <- attr(X, "scale:mode")
    if (identical(scale_mode, "matrix")) {
        scale_factor <- attr(X, "scale:factor")
        mat <- t(solve(t(scale_factor), t(mat)))
    } else if (identical(scale_mode, "vector")) {
        mat <- sweep(mat, 2L, attr(X, "scale:factor"), "*")
    }

    mask <- attr(X, "mask")
    if  (is.null(mask))
        mask <- rep(TRUE, ncol(mat))  # no filtering

    x_mean <- attr(X, "restore:center")
    if (is.null(x_mean)) {
        x_mean <- rep(0, length(mask))
        names(x_mean) <- x_names
    }
    out <- matrix(x_mean, nrow = nrow(mat), ncol = length(x_mean), byrow = TRUE,
                  dimnames = list(NULL, names(x_mean)))

    if (identical(scale_mode, "vector")) {
        x_center <- attr(X, "scaled:center")
        if (is.null(x_center)) {
            x_center <- rep(0, length(mask))
            names(x_center) <- names(x_mean)
        }
        if (!is.null(attr(X, "mask"))) {
            x_center <- x_center[mask]
        }
        out[, mask] <- matrix(x_center, nrow = nrow(mat), ncol = ncol(mat),
                              byrow = TRUE) + mat
    } else {
        out[, mask] <- mat
    }

    out
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
    if (length(entries[["i"]]) == 0L)
        return(X)
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
    if (identical(scale_mode, "none") && !is.null(data_scale))
        return(data_scale)

    if (identical(scale_mode, "vector")) {
        scale_factor <- attr(X, "scale:factor")
        if (!is.null(data_scale) && !is.null(scale_factor))
            return(data_scale / scale_factor)
    }

    matrixStats::colSds(X)
}

.aa_auto_bigM <- function(X, multiplier = 200) {
    feature_sd <- .aa_fit_space_col_sds(X)
    feature_scale <- sqrt(mean(feature_sd^2))
    if (!is.finite(feature_scale) || feature_scale <= .Machine$double.eps)
        feature_scale <- 1
    multiplier * feature_scale * sqrt(2 / max(ncol(X), 1L))
}

# Common subroutine to preprocess `init` matrix of archetype coordinates (A);
# maps `init` from original data space to preprocessed space of X.
.aa_preprocess_init <- function(init, X) {
    # If `init` is not a matrix, return it as is (it will be processed by the init function)
    if (!(is.matrix(init) || inherits(init, "data.frame")))
        return(init)

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
    if (!is.null(x_center))
        init <- sweep(init, 2L, x_center, "-")

    # Filter `init` to match filtered features in X
    if (!is.null(mask) && ncol(init) == length(mask))
        init <- init[, mask, drop = FALSE]

    if (identical(attr(X, "scale:mode"), "matrix"))
        init <- init %*% attr(X, "scale:factor")
    if (identical(attr(X, "scale:mode"), "vector"))
        init <- sweep(init, 2L, attr(X, "scale:factor"), "/")

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
    if (is.null(nm))
        return(paste0("A", seq_len(nrow(A))))

    stopifnot("Archetype names must not be missing" = !any(is.na(nm)))
    stopifnot("Archetype names must not be empty" = all(nzchar(nm)))
    stopifnot("Archetype names must be unique" = !anyDuplicated(nm))
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
    B <- fit_qp(
        A = X,
        X = init,
        eps = eps,
        project = if (delta == 0) proj_l1 else NULL,
        row_sum_bounds = c(a_lo, a_hi)
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
        S = .init_S(X, A, eps = eps),
        loss = .aa_new_loss(L)
    )
}

.aa_trace_row_rss <- function(row_xss, S, XAt, AAt) {
    pmax(row_xss - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt)), 0)
}

.aa_update_loss <- function(loss, i, loss_terms, verbose, max_kappa = 1,
                            k_A = c("exact", "gram")) {
    k_A <- match.arg(k_A)
    # Update loss metrics: add row "i" (current iteration) to the loss dataframe
    loss[["loss"]][i] <- loss_terms[["rss"]]
    loss[["r2"]][i]  <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

    # Compute condition numbers
    if ((i - 1L) %% 10 != 0)
        return(loss)  # only update kappa every 10 iterations for efficiency

    if (max_kappa > 1) {
        # S tends to be long and skinny (N >> K), so use rcond for better stability
        # A tends to be roughly square (K ~ M), so use kappa directly
        loss[["k_S"]][i] <- sqrt(1 / rcond(loss_terms[["StS"]]))
        if (identical(k_A, "gram")) {
            loss[["k_A"]][i] <- sqrt(1 / rcond(loss_terms[["AAt"]]))
        } else if (is.na(loss[["k_A"]][i])) {
            loss[["k_A"]][i] <- kappa(loss_terms[["A"]], exact = TRUE)
        }
    } else {
        loss[["k_S"]][i] <- max_kappa
        if (is.na(loss[["k_A"]][i])) loss[["k_A"]][i] <- max_kappa
    }
    loss
}

# Check convergence based on relative loss improvement and R2 threshold
.aa_check_convergence <- function(loss, i, tol, tol_r2, max_kappa, verbose) {
    j <- i + 1L  # save some typing...
    # Main
    converged <- with(
        loss,
        abs(loss[j] - loss[i]) < tol * loss[i] || r2[j] > tol_r2
    )
    # Warnings/Messages
    if (verbose && i %% 10 == 0) {
        fmt <- "Iteration %d: loss = %.4f, R2 = %.3f"
        message(sprintf(fmt, i, loss[["loss"]][j], loss[["r2"]][j]))
    }
    k_S <- loss[["k_S"]][j]
    k_A <- loss[["k_A"]][j]
    if ((!is.na(k_S) && k_S > max_kappa) || (!is.na(k_A) && k_A > max_kappa)) {
        fmt <- "Warning: Condition number exceeded max_kappa (k_S=%.1f, k_A=%.1f)"
        warning(sprintf(fmt, k_S, k_A), call. = FALSE)
    }
    converged
}
