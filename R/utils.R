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
            warning("Centering matrices breaks sparsity; consider using `center = FALSE`")
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

.aa_check_inputs <- function(data, tol, tol_r2, K, max_kappa, eps, robust, tukey_c, scale = TRUE) {
    stopifnot("data contains missing values" = !any(is.na(data)))
    stopifnot("tol must be positive" = tol > 0)
    stopifnot("tol_r2 must be between (0, 1)" = tol_r2 >= 0 && tol_r2 <= 1)
    stopifnot("K must be an integer" = K == as.integer(K))
    stopifnot("K must be an integer greater or equal to 1" = K >= 1L)
    stopifnot("K cannot be greater than number of samples" = K <= nrow(data))
    stopifnot("max_kappa must be >=1" = max_kappa >= 1)
    stopifnot("eps must be non-negative" = eps >= 0)
    stopifnot("robust must be TRUE or FALSE" =
                  is.logical(robust) && length(robust) == 1L && !is.na(robust))
    stopifnot("tukey_c must be positive" =
                  length(tukey_c) == 1L && is.finite(tukey_c) && tukey_c > 0)
    .aa_check_scale(scale, ncol(data))
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
    loss <- data.frame(rss = 0, r2 = 1, k_S = 1, k_A = kappa(A))

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
    loss <- data.frame(rss = rss, r2 = 1 - rss / xss, k_S = 1, k_A = 1)

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
        warning(sprintf(fmt, ncol(X) - M, dropped_features))

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


# Common subroutine to preprocess input data matrix:
# - optionally z-score or apply a matrix metric embedding
# - filter out low-variance features
# - apply user-provided sample weights (if any)
# - add bigM intercept term (if bigM > 0) to "force" the simplex constraint during nnls fit
# Returns a list with:
# - X: preprocessed data matrix with attributes to undo scaling and filtering
# - undo_scale: function to undo scaling and filtering of archetype coordinates
.aa_preprocess <- function(data, sd_threshold, weights, verbose, bigM = 0, scale = TRUE) {
    if (verbose) message("Preprocessing data...")

    original_center <- colMeans(data)
    names(original_center) <- colnames(data)
    scale_mode <- if (identical(scale, FALSE)) {
        "none"
    } else if (isTRUE(scale) || (is.numeric(scale) && is.null(dim(scale)))) {
        "vector"
    } else {
        "matrix"
    }

    X <- if (inherits(data, "sparseMatrix")) {
        data
    } else {
        as.matrix(data)
    }
    if (isTRUE(scale)) {
        n <- nrow(data)
        x_mean <- colMeans(data)
        x2 <- colSums(data * data)
        x_var <- pmax((x2 - n * x_mean * x_mean) / max(n - 1L, 1L), 0)
        scale <- sqrt(x_var)
        names(scale) <- colnames(data)

        if (!inherits(data, "sparseMatrix")) {
            X <- sweep(X, 2L, x_mean, "-")
            attr(X, "scaled:center") <- x_mean
        }
        attr(X, "scaled:scale") <- scale
    }

    # Filter out low-variance features
    X <- .filter_low_variance(X, sd_threshold)
    mask <- attr(X, "mask")

    if (identical(scale_mode, "matrix")) {
        if (!is.null(mask))
            scale <- scale[mask, mask, drop = FALSE]
        scale_factor <- t(chol(as.matrix(scale)))
        X <- as.matrix(X) %*% scale_factor
        retained_names <- names(original_center)
        if (!is.null(mask))
            retained_names <- retained_names[mask]
        colnames(X) <- retained_names
        attr(X, "mask") <- mask
        attr(X, "scale:factor") <- scale_factor
    } else if (identical(scale_mode, "vector")) {
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
    }
    attr(X, "scale:mode") <- scale_mode
    attr(X, "restore:center") <- original_center
    N <- nrow(X) # number of samples
    if (is.null(bigM))
        bigM <- .aa_auto_bigM(X)

    # Weight samples by user-provided "importance" weights
    if (!is.null(weights)) {
        # Check weights sanity
        if (length(weights) != N) {
            fmt <- "Number of weights (%d) must equal number of rows in data (%d)"
            stop(sprintf(fmt, length(weights), N))
        }
        stopifnot("Weights contain NA values" = !any(is.na(weights)))
        stopifnot("Weights must be non-negative" = all(weights >= 0))
        # Rescale and apply weights
        weights <- weights / mean(weights) # normalize weights to mean 1
        x_attrs <- attributes(X)
        X <- X * weights
        attributes(X) <- utils::modifyList(attributes(X), x_attrs)
        attr(X, "weights") <- weights # store weights in X attributes
    }

    if (bigM > 0) {
        # add bigM intercept term to "force" the simplex constraint during nnls fit
        x_attrs <- attributes(X)
        bigM_col <- matrix(bigM, nrow = N, ncol = 1L, dimnames = list(rownames(X), "bigM"))
        if (inherits(X, "sparseMatrix"))
            bigM_col <- as(bigM_col, "sparseMatrix")
        X <- cbind(bigM_col, X)
        # Restore attributes
        attr(X, "scaled:center") <- x_attrs[["scaled:center"]]
        attr(X, "scaled:scale")  <- x_attrs[["scaled:scale"]]
        attr(X, "scale:mode") <- x_attrs[["scale:mode"]]
        attr(X, "scale:factor") <- x_attrs[["scale:factor"]]
        attr(X, "restore:center") <- x_attrs[["restore:center"]]
        attr(X, "mask") <- x_attrs[["mask"]]
        attr(X, "bigM")  <- 1L
        attr(X, "bigM.value") <- bigM
    }

    # To undo scaling when returning archetypes
    undo_scale <- function(mat, X) {
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

    list(X = X, undo_scale = undo_scale)
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
        warning(sprintf(fmt, paste(utils::head(ix, 10L), collapse = ", ")))
    }

    rownames(A) <- rownames(B) <- nm
    colnames(B) <- rownames(X)
    list(
        A = A,
        B = B,
        S = .init_S(X, A, eps = eps),
        loss = list(rss = rep(NA_real_, L),
                    r2  = rep(NA_real_, L),
                    k_S = rep(NA_real_, L),
                    k_A = rep(NA_real_, L))
    )
}

# Common subroutine to initialize variables:
# archetype coordinates (A), coefficients (B), compositions (S) and
# loss metrics: rss, r2, condition numbers of S and A (k_S, k_A)
.aa_init_vars <- function(X, K, init, init_args, eps, max_iter, verbose, delta = 0) {
    if (verbose) message("Initializing archetypes...")
    L <- max_iter + 1L

    # `init` is fixed coordinates of archetypes: call .aa_matrix_init
    if (is.matrix(init) || inherits(init, "data.frame")) {
        # use provided coordinate matrix as initialization; ignore `init_args`
        init <- .aa_preprocess_init(init, X)
        if (length(init_args) > 0L) {
            warning("`init_args` are ignored when `init` is a matrix")
            init_args <- list()
        }
        init_vars <- .aa_matrix_init(X, K, init, eps, L, delta)
        return(init_vars)
    }

    # init_vars must be generated from data
    if (is.character(init)) {
        # use `aa_init` function with method specified by `init` string
        stopifnot("`init` must be a single string" = length(init) == 1L)
        init_args <- c(list(method = init), init_args)
        init <- aa_init
    } else if (!is.function(init)) {
        stop("`init` must be a function, a single string, or archetypes coordinate matrix")
    }
    init_vars <- do.call(init, args = c(list(X = X, K = K), init_args))
    rownames(init_vars[["A"]]) <- .aa_init_names(init_vars[["A"]])
    rownames(init_vars[["B"]]) <- rownames(init_vars[["A"]])
    init_vars[["S"]] <- .init_S(X, init_vars[["A"]], eps = eps)
    init_vars[["loss"]] <- list(rss = rep(NA_real_, L),
                                r2  = rep(NA_real_, L),
                                k_S = rep(NA_real_, L),
                                k_A = rep(NA_real_, L))
    init_vars
}

# Workhorse function to run archetypes fitting with method-specific arguments;
# Computes the loss metrics per row and for the whole dataset without forming
# the full residual matrix. Interemediate parts can be cached in the loop for efficiency.
# return_S_terms can skip computing the expensive StS and StX when not needed for the fit (nnls case)
# TODO: consider refactoring into: cache, weight, objective, s_terms
.aa_loss_terms <- function(X, A, S, weight_fun = NULL,
                           return_S_terms = TRUE,
                           xss = NULL,
                           rss = NULL,
                           row_xss = NULL,
                           row_rss = NULL,
                           row_weights = NULL,
                           StS = NULL,
                           StX = NULL,
                           AAt = NULL,
                           XAt = NULL) {
    # TODO: instead of creating copies consider correcting the affected terms (used to compute rss, xss)
    iM <- attr(X, "bigM")
    if (!is.null(iM)) {
        X <- X[, -iM, drop = FALSE]
        A <- A[, -iM, drop = FALSE]
    }
    update_row_weights <- !is.null(weight_fun)
    if (update_row_weights) {
        xss <- NULL
        rss <- NULL
        row_weights <- NULL
        StS <- NULL
        StX <- NULL
    } else if (!is.null(row_weights)) {
        .aa_check_row_weights(row_weights, nrow(X))
        if (all(row_weights == 1))
            row_weights <- NULL
    }

    # row_xss does not change between iterations, so compute it once and cache for efficiency
    if (is.null(row_xss)) row_xss <- rowSums(X * X)

    # A-terms
    if (is.null(AAt)) AAt <- tcrossprod(A)
    if (is.null(XAt)) XAt <- tcrossprod(X, A)

    # Residual terms
    if (is.null(row_rss))
        row_rss <- pmax(row_xss - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt)), 0)
    if (update_row_weights) {
        new_row_weights <- weight_fun(row_rss)
        .aa_check_row_weights(new_row_weights, nrow(X))
        if (!.aa_trivial_row_weights(new_row_weights)) {
            row_weights <- new_row_weights
            xss <- sum(row_weights * row_xss)  # weighted total sum of squares
        }
    }
    if (is.null(xss)) xss <- sum(.aa_weight_rows(row_xss, row_weights))
    if (is.null(rss)) rss <- sum(.aa_weight_rows(row_rss, row_weights))

    S_weighted <- .aa_weight_rows(S, row_weights)
    if (return_S_terms) {
        if (is.null(StX)) StX <- crossprod(S_weighted, X)
        if (is.null(StS)) StS <- crossprod(S_weighted, S)
    }
    if (!return_S_terms) {
        StX <- NULL
        StS <- NULL
    }

    list(
        rss = rss,
        xss = xss,
        row_xss = row_xss,
        row_rss = row_rss,
        row_weights = row_weights,
        StS = StS,
        StX = StX,
        S_weighted = S_weighted,
        AAt = AAt,
        XAt = XAt,
        A = A
    )
}


.aa_update_loss <- function(loss, i, loss_terms, verbose, max_kappa = 1) {
    # Update loss metrics: add row "i" (current iteration) to the loss dataframe
    loss[["rss"]][i] <- loss_terms[["rss"]]
    loss[["r2"]][i]  <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]

    # Compute condition numbers
    if ((i - 1L) %% 10 != 0)
        return(loss)  # only update kappa every 10 iterations for efficiency

    if (max_kappa > 1) {
        # S tends to be long and skinny (N >> K), so use rcond for better stability
        # A tends to be roughly square (K ~ M), so use kappa directly
        loss[["k_S"]][i] <- sqrt(1 / rcond(loss_terms[["StS"]]))
        if (is.na(loss[["k_A"]][i]))
            loss[["k_A"]][i] <- kappa(loss_terms[["A"]], exact = TRUE)
    } else {
        loss[["k_S"]][i] <- max_kappa
        if (is.na(loss[["k_A"]][i])) loss[["k_A"]][i] <- max_kappa
    }
    loss
}

# Check convergence based on relative RSS improvement and R2 threshold
.aa_check_convergence <- function(loss, i, tol, tol_r2, max_kappa, verbose) {
    j <- i + 1L  # save some typing...
    # Main
    converged <- with(
        loss,
        abs(rss[j] - rss[i]) < tol * rss[i] || r2[j] > tol_r2
    )
    # Warnings/Messages
    if (verbose && i %% 10 == 0) {
        fmt <- "Iteration %d: RSS = %.4f, R2 = %.3f"
        message(sprintf(fmt, i, loss[["rss"]][j], loss[["r2"]][j]))
    }
    k_S <- loss[["k_S"]][j]
    k_A <- loss[["k_A"]][j]
    if ((!is.na(k_S) && k_S > max_kappa) || (!is.na(k_A) && k_A > max_kappa)) {
        fmt <- "Warning: Condition number exceeded max_kappa (k_S=%.1f, k_A=%.1f)"
        warning(sprintf(fmt, k_S, k_A))
    }
    converged
}

# Format output of archetypes fitting into "archetypes" S3 object
.aa_prepare_output <- function(X, A, B, S,
                               i, loss, converged,
                               undo_scale, max_iter, verbose, delta = 0,
                               data = NULL, call = NULL, A0 = NULL) {

    # Format loss dataframe: keep only rows up to iteration "i" and reset row names
    j <- i + 1L
    loss <- as.data.frame(loss)[1:j, , drop = FALSE]
    rownames(loss) <- NULL

    # Undo scaling of archetype coordinates to return them in the original data space
    archetype_names <- if (!is.null(A0)) rownames(A0) else rownames(A)
    A <- undo_scale(A, X)
    if (!is.null(A0))
        A0 <- undo_scale(A0, X)

    # set row and column names for output matrices
    if (is.null(archetype_names))
        archetype_names <- paste0("A", seq_len(nrow(A)))
    rownames(A) <- rownames(B) <- colnames(S) <- archetype_names
    if (!is.null(A0))
        rownames(A0) <- archetype_names
    colnames(B) <- rownames(S) <- rownames(X)

    # Warn or message about convergence
    if (!converged) {
        fmt <- "Algorithm did not converge after %d iterations"
        warning(sprintf(fmt, max_iter))
    }

    if (verbose) {
        fmt <- ifelse(converged,
                      "Converged after %d iterations:",
                      "Final iteration %d:")
        fmt <- paste(fmt, "loss = %.4g, R2 = %.3f")
        message(sprintf(fmt, i, loss[j, "rss"], loss[j, "r2"]))
    }

    archetypes(
        call         = call,
        data         = data,
        init         = A0,
        coordinates  = A,
        coefficients = B,
        compositions = S,
        slack        = delta,
        loss         = loss,
        converged    = converged
    )
}
