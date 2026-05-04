# Tukey's Bisquare Weight Function
#
# Computes Tukey's bisquare (biweight) robust weights for outlier downweighting.
#
# @param resid matrix of residuals (same dimension as X)
# @param c tuning constant (default: 4.685 for 95% efficiency)
#
# @returns Matrix of weights between 0 and 1, same shape as resid
#
# @details The bisquare weight function is:
# w(u) = (1 - (u/c)^2)^2 if |u| <= c, 0 otherwise
# where u = row-wise standardized residuals.
#
bisquare0 <- function(resid, c = 4.685) {
    # Row-wise MAD
    mad_row <- apply(resid, 1, mad, na.rm = TRUE)
    mad_row[mad_row == 0] <- 1  # avoid division by zero

    # Standardize each row
    u <- resid / mad_row

    weights <- ifelse(abs(u) <= c, (1 - (u/c)^2)^2, 0)
    return(weights)
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
    Sx <- cov(X)  # (M×M)

    # Try using Cholesky to invert PD matrix
    cx <- tryCatch(chol(Sx), error = function(e) NULL)
    if (is.null(cx))  # FALLBACK to general solution
        return(sum(diag( solve(Sx, cov(Y)) )))


    # Edge Case: # of samples smaller than the # of features
    if (N < M) {
        # avoid computing Sy and use the Mahalanobis approach instead.
        Yc <- sweep(Y, 2L, colMeans(Y), FUN = "-", check.margin = FALSE)
        Z  <- backsolve(cx, t(Yc))
        # remember cov(Y) uses denominator (N−1)
        return(norm(Z, type = "F")^2 / (N - 1))
    }
    # since trace(AB) = sum(A * t(B)) = sum(A * B) when symmetric
     sum(cov(Y) * chol2inv(cx))
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
        sd_vals <- apply(X, 2, sd)
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

.aa_init_vars <- function(X, K, init, init_args, eps, max_iter, verbose) {

    if (verbose) message("Initializing archetypes...")
    L <- max_iter + 1L
    if (is.character(init)) {
        stopifnot("`init` must be a single string" = length(init) == 1L)
        init_args <- c(list(method = init), init_args)
        init <- aa_init
    } else {
        stopifnot("`init` must be a function or a single string" = is.function(init))
    }
    init_vars <- do.call(init, args = c(list(X = X, K = K), init_args))
    init_vars[["S"]] <- .init_S(X, init_vars[["A"]], eps = eps)
    init_vars[["loss"]] <- list(rss = rep(NA_real_, L),
                                r2  = rep(NA_real_, L),
                                k_S = rep(NA_real_, L),
                                k_A = rep(NA_real_, L))
    return(init_vars)
}

.aa_update_loss <- function(loss, i, verbose, max_kappa = 1, X = NULL, ...) {
    res <- if (is.null(X)) .aa_loss_pdg(...) else .aa_loss_nnls(X, ...)
    loss[["rss"]][i] <- res[["rss"]]
    loss[["r2"]][i]  <- 1 - res[["rss"]] / res[["xss"]]
    if (i %% 10 != 0) {
        loss[["k_S"]][i] <- NA_real_
        loss[["k_A"]][i] <- NA_real_
    } else if (max_kappa > 1) {
        loss[["k_S"]][i] <- sqrt(1/rcond(res[["StS"]]))
        A <- res[["A"]]
        Gram <- if (nrow(A) < ncol(A)) res[["AAt"]] else crossprod(A)
        loss[["k_A"]][i] <- sqrt(1/rcond(Gram))
    } else {
        loss[["k_S"]][i] <- max_kappa
        loss[["k_A"]][i] <- max_kappa
    }
    return(loss)
}

.aa_loss_nnls <- function(X, xss, A, S,...) {
    # Compute RSS without forming the full residual matrix,
    # using the cosine law: ||X||^2 + ||AS||^2 - 2*trace(S %*% t(A) %*% X)
    iM <- attr(X, "bigM")
    if (!is.null(iM)) {
        # if bigM column is present, we need to remove it from the loss computation
        A <- A[, -iM, drop = FALSE]
        X <- X[, -iM, drop = FALSE]
    }
    AAt <- tcrossprod(A)
    StX <- crossprod(S, X)
    StS <- crossprod(S)
    rss <- xss - 2 * sum(A * StX) + sum(StS * AAt)
    list(rss = rss, xss = xss, StS = StS, AAt = AAt, A = A)
}

.aa_loss_pdg <- function(rss, xss, StS, AAt, A = NULL, ...) {
    # Just reformat the arguments into a consistent list format for the loss update function
    list(rss = rss, xss = xss, StS = StS, AAt = AAt, A = A)
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
    loss <- head(as.data.frame(loss), j)
    rownames(loss) <- NULL

    A <- undo_scale(A, X)
    if (!is.null(A0))
        A0 <- undo_scale(A0, X)

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
