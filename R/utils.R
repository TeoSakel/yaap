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

.aa_unit_weights <- function(row_rss) rep(1, length(row_rss))

.aa_check_row_weights <- function(row_weights, n) {
    stopifnot("row_weights must match rows in X" = length(row_weights) == n)
    stopifnot("row_weights contain NA values" = !any(is.na(row_weights)))
    stopifnot("row_weights must be non-negative" = all(row_weights >= 0))
    invisible(TRUE)
}

.aa_loss_terms <- function(X, A, S, weight_fun,
                           return_S_terms = TRUE,
                           xss = NULL,
                           rss = NULL,
                           x2 = NULL,
                           row_rss = NULL,
                           row_weights = NULL,
                           StS = NULL,
                           StX = NULL,
                           AAt = NULL,
                           XAt = NULL) {
    iM <- attr(X, "bigM")
    if (!is.null(iM)) {
        X <- X[, -iM, drop = FALSE]
        A <- A[, -iM, drop = FALSE]
    }
    if (is.null(AAt))
        AAt <- tcrossprod(A)
    if (is.null(XAt))
        XAt <- tcrossprod(X, A)
    if (is.null(x2))
        x2 <- matrixStats::rowSums2(X * X)
    if (is.null(row_rss)) {
        row_rss <- x2 - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt))
        row_rss <- pmax(row_rss, 0)
    }

    if (is.null(row_weights))
        row_weights <- weight_fun(row_rss)
    .aa_check_row_weights(row_weights, nrow(X))

    if (is.null(xss))
        xss <- sum(row_weights * x2)
    if (is.null(rss))
        rss <- sum(row_weights * row_rss)
    if (return_S_terms && (is.null(StX) || is.null(StS))) {
        S_weighted <- S * row_weights
        if (is.null(StX))
            StX <- crossprod(S_weighted, X)
        if (is.null(StS))
            StS <- crossprod(S_weighted, S)
    }
    if (!return_S_terms) {
        StX <- NULL
        StS <- NULL
    }

    list(
        rss = rss,
        xss = xss,
        x2 = x2,
        row_rss = row_rss,
        row_weights = row_weights,
        StS = StS,
        StX = StX,
        AAt = AAt,
        XAt = XAt,
        A = A
    )
}

# Compute squared Euclidean distance of each sample from center
.dist2 <- function(X, center = FALSE) {
    d <- scale(X, center = center, scale = FALSE)  # shift columns by a
    matrixStats::rowSums2(d*d)
}


# Compute pairwise squared distances between columns of X and Y
.pdist2 <- function(X, Y) {
    # both x and y must be in the same dimensional space
    stopifnot(ncol(X) == ncol(Y))

    # squared distances from cosine law ||x||^2 + ||y||^2 − 2<x,y>
    sx <- matrixStats::rowSums2(X * X)   # ||x||^2
    sy <- matrixStats::rowSums2(Y * Y)   # ||y||^2
    cp <- tcrossprod(X, Y)  # 2*<x,y>
    D2 <- outer(sx, sy, "+") - 2*cp
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
        return(sum(diag( solve(Sx, stats::cov(Y)) )))


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

# Blocks of code for archetypes fitting ---------------------------------------


.aa_check_inputs <- function(data, tol, tol_r2, K, max_kappa, eps) {
    stopifnot("data contains missing values" = !any(is.na(data)))
    stopifnot("tol must be positive" = tol > 0)
    stopifnot("tol_r2 must be between (0, 1)" = tol_r2 >= 0 && tol_r2 <= 1)
    stopifnot("K must be an integer" = K == as.integer(K))
    stopifnot("K must be an integer greater or equal to 1" = K >= 1L)
    stopifnot("K cannot be greater than number of samples" = K <= nrow(data))
    stopifnot("max_kappa must be >=1" = max_kappa >= 1)
    stopifnot("eps must be non-negative" = eps >= 0)
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
    return(out)
}

# Edge case K == N: each sample is its own archetype
.identity_archetypes <- function(X, call = NULL) {
    A <- X
    rownames(A) <- paste0("A", seq_len(nrow(X)))  # remove row names for consistency
    S <- B <- diag(nrow(X))
    colnames(B) <- rownames(S) <- rownames(X)
    rownames(B) <- colnames(S) <- rownames(A)
    loss <- data.frame(rss = 0, r2 = 1, k_S = 1, k_A = kappa(A))

    out <- archetypes(
        call         = NULL,
        data         = X,
        init         = A,
        coordinates  = A,
        coefficients = B,
        compositions = S,
        loss         = loss,
        converged    = TRUE
    )
    return(out)
}

# Edge case K == 1: single archetype at mean of X
.mean_archetype <- function(X) {
    x_mean <- matrixStats::colMeans2(X)
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
    loss <- data.frame(rss = rss, rsq = 1 - rss / xss, k_S = 1, k_A = 1)

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

# Filter out features (columns of X) with low variance
.filter_low_variance <- function(X, sd_threshold) {
    sd_vals <- attr(X, "scaled:scale")
    if (is.null(sd_vals))
        sd_vals <- apply(X, 2, stats::sd)
    mask <- sd_vals >= sd_threshold
    M <- sum(mask)
    if (M < ncol(X)) {
        # Throw Warning
        fmt <- "The following %d features are filtered out due to low variance: %s"
        dropped_features <- if (is.null(colnames(X))) which(!mask) else colnames(X)[!mask]
        dropped_features <- paste(dropped_features, collapse = ", ")
        warning(sprintf(fmt, ncol(X) - M, dropped_features))

        # Filter X
        x_attrs <- attributes(X)
        x_attrs[["mask"]] <- mask
        X <- X[, mask, drop = FALSE]
        attributes(X) <- x_attrs
    }
    X
}


.aa_preprocess <- function(data, sd_threshold, weights, verbose, bigM = 0) {

    if (verbose) message("Preprocessing data...")

    # Scale input matrix to 0 mean and unit variance
    X <- scale(as.matrix(data))

    # Filter out low-variance features
    X <- .filter_low_variance(X, sd_threshold)
    N <- nrow(X) # number of samples

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
        X <- X * weights
        attr(X, "weights") <- weights # store weights in X attributes
    }

    xss <- norm(X, type = "F")^2

    if (bigM > 0) {
        # add bigM intercept term to "force" the simplex constraint during nnls fit
        x_attrs <- attributes(X)
        X <- cbind(
            matrix(bigM, nrow = N, ncol = 1L, dimnames = list(rownames(X), "bigM")),
            X
        )
        # Restore attributes
        attr(X, "scaled:center") <- x_attrs[["scaled:center"]]
        attr(X, "scaled:scale")  <- x_attrs[["scaled:scale"]]
        attr(X, "bigM")  <- 1L
        attr(X, "bigM.value") <- bigM
    }

    # To undo scaling when returning archetypes
    undo_scale <- function(mat, X) {
        stopifnot(ncol(mat) == ncol(X))

        # Remove bigM if present
        iM <- attr(X, "bigM")
        if (!is.null(iM))
            mat <- mat[, -iM, drop = FALSE]

        # Undo scaling
        x_mean <- attr(X, "scaled:center")
        if (is.null(x_mean)) {
            x_mean <- rep(0, ncol(mat)) # no centering
            names(x_mean) <- colnames(X)
        }
        x_std <- attr(X, "scaled:scale")
        if (is.null(x_std))
            x_std <- rep(1, ncol(mat)) # no scaling
        mask <- attr(X, "mask")
        if  (is.null(mask))
            mask <- rep(TRUE, ncol(mat))  # no filtering

        M <- length(x_mean)
        out <- matrix(x_mean, nrow = nrow(mat), ncol = M, byrow = TRUE,
                      dimnames = list(NULL, names(x_mean)))
        out[, mask] <- out[, mask, drop = FALSE] + mat * x_std

        out
    }

    list(X = X, undo_scale = undo_scale, xss = xss)
}

.aa_preprocess_init <- function(init, data, X) {
    if (!is.matrix(init))
        return(init)

    if (ncol(init) != ncol(data)) {
        fmt <- "ncol(init) = %d does not match number of data features (%d)"
        stop(sprintf(fmt, ncol(init), ncol(data)))
    }

    init <- as.matrix(init)
    x_center <- attr(X, "scaled:center")
    x_scale <- attr(X, "scaled:scale")
    if (!is.null(x_center))
        init <- sweep(init, 2L, x_center, "-")
    if (!is.null(x_scale))
        init <- sweep(init, 2L, x_scale, "/")

    mask <- attr(X, "mask")
    if (!is.null(mask))
        init <- init[, mask, drop = FALSE]

    iM <- attr(X, "bigM")
    if (!is.null(iM)) {
        init <- cbind(
            matrix(attr(X, "bigM.value"), nrow = nrow(init), ncol = 1L),
            init
        )
        colnames(init)[iM] <- "bigM"
    }

    init
}

.aa_init_names <- function(A) {
    nm <- rownames(A)
    if (is.null(nm))
        return(paste0("A", seq_len(nrow(A))))

    stopifnot("Archetype names must not be missing" = !any(is.na(nm)))
    stopifnot("Archetype names must not be empty" = all(nzchar(nm)))
    stopifnot("Archetype names must be unique" = !anyDuplicated(nm))
    nm
}

.aa_matrix_init <- function(X, K, init, eps, delta = 0, tol = 1e-6) {
    if (nrow(init) != K) {
        fmt <- "nrow(init) = %d does not match K (%d)"
        stop(sprintf(fmt, nrow(init), K))
    }
    if (ncol(init) != ncol(X)) {
        fmt <- "ncol(init) = %d does not match preprocessed data features (%d)"
        stop(sprintf(fmt, ncol(init), ncol(X)))
    }

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
    list(A = A, B = B)
}

.aa_init_vars <- function(X, K, init, init_args, eps, max_iter, verbose, delta = 0) {

    if (verbose) message("Initializing archetypes...")
    L <- max_iter + 1L
    if (is.character(init)) {
        stopifnot("`init` must be a single string" = length(init) == 1L)
        init_args <- c(list(method = init), init_args)
        init <- aa_init
    } else if (is.matrix(init)) {
        if (length(init_args) > 0L) {
            warning("`init_args` are ignored when `init` is a matrix")
            init_args <- list()
        }
        init_vars <- .aa_matrix_init(X, K, init, eps, delta)
        init_vars[["S"]] <- .init_S(X, init_vars[["A"]], eps = eps)
        init_vars[["loss"]] <- list(rss = rep(NA_real_, L),
                                    r2  = rep(NA_real_, L),
                                    k_S = rep(NA_real_, L),
                                    k_A = rep(NA_real_, L))
        return(init_vars)
    } else {
        stopifnot("`init` must be a function, a single string, or a coordinate matrix" = is.function(init))
    }
    init_vars <- do.call(init, args = c(list(X = X, K = K), init_args))
    rownames(init_vars[["A"]]) <- .aa_init_names(init_vars[["A"]])
    rownames(init_vars[["B"]]) <- rownames(init_vars[["A"]])
    init_vars[["S"]] <- .init_S(X, init_vars[["A"]], eps = eps)
    init_vars[["loss"]] <- list(rss = rep(NA_real_, L),
                                r2  = rep(NA_real_, L),
                                k_S = rep(NA_real_, L),
                                k_A = rep(NA_real_, L))
    return(init_vars)
}

.aa_update_loss <- function(loss, i, loss_terms, verbose, max_kappa = 1) {
    loss[["rss"]][i] <- loss_terms[["rss"]]
    loss[["r2"]][i]  <- 1 - loss_terms[["rss"]] / loss_terms[["xss"]]
    if (i %% 10 != 0) {
        loss[["k_S"]][i] <- NA_real_
        loss[["k_A"]][i] <- NA_real_
    } else if (max_kappa > 1) {
        loss[["k_S"]][i] <- sqrt(1/rcond(loss_terms[["StS"]]))
        A <- loss_terms[["A"]]
        Gram <- if (nrow(A) < ncol(A)) loss_terms[["AAt"]] else crossprod(A)
        loss[["k_A"]][i] <- sqrt(1/rcond(Gram))
    } else {
        loss[["k_S"]][i] <- max_kappa
        loss[["k_A"]][i] <- max_kappa
    }
    return(loss)
}

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

.aa_prepare_output <- function(X, A, B, S,
                               i, loss, converged,
                               undo_scale, max_iter, verbose, delta = 0,
                               data = NULL, call = NULL, A0 = NULL) {

    j <- i + 1L
    loss <- as.data.frame(loss)[1:j, , drop = FALSE]
    rownames(loss) <- NULL

    archetype_names <- if (!is.null(A0)) rownames(A0) else rownames(A)
    A <- undo_scale(A, X)
    if (!is.null(A0))
        A0 <- undo_scale(A0, X)

    if (is.null(archetype_names))
        archetype_names <- paste0("A", seq_len(nrow(A)))
    rownames(A) <- rownames(B) <- colnames(S) <- archetype_names
    if (!is.null(A0))
        rownames(A0) <- archetype_names
    colnames(B) <- rownames(S) <- rownames(X)

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

    out <- archetypes(
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
    out
}
