#' Project rows of matrix onto the probability simplex
#'
#' @description These functions project each row of a numeric matrix onto the
#' probability simplex or the L1-norm ball respectively by implementing Condat's
#' algorithm.
#'
#' @param mat A numeric matrix where each row will be projected.
#' @param eps A small positive number to ensure numerical stability (default: 1e-8).
#'
#' @return A row-stochastic matrix.
#'
#' @examples
#' mat <- matrix(runif(12), nrow = 4)
#'
#' proj_simplex(mat)
#'
#' proj_l1(mat)
#'
#' @details
#' The values of the input matrix are first clipped at `eps` to be non-negative
#' and then projected onto the simplex or L1-norm ball.
#'
#' @references
#' Condat, L. (2016). Fast projection onto the simplex and the L_1 ball.
#' Mathematical Programming, 158(1), 575-585.
#' \url{http://dx.doi.org/10.1007/s10107-015-0946-6}
#'
#' @name simplex_projection
NULL

#' @rdname simplex_projection
#' @export
proj_simplex <- function(mat, eps = 0) {
    stopifnot("eps must be non-negative" = eps >= 0)
    proj_row <- function(v) {
        sv <- sum(v)
        if (abs(sv - 1) <= length(v) * eps) return(v / sv)
        u     <- sort(v, decreasing = TRUE)
        cssv  <- cumsum(u)
        rho   <- max(which(u > (cssv - 1) / seq_along(u)))
        theta <- (cssv[rho] - 1) / rho
        pmax(v - theta, 0)
    }
    out <- pmax(mat, eps)
    for (i in seq_len(nrow(out)))  # TODO: implement in parallel
        out[i, ] <- proj_row(out[i, ])
    out
}


#' @rdname simplex_projection
#' @export
proj_l1 <- function(mat, eps = 0) {
    stopifnot("eps must be non-negative" = eps >= 0)
    out <- pmax(mat, eps)      # project to 1st orthant
    out <- out / rowSums(out)  # l1-normalize each row
    out
}


#' Fit Data to Convex Hull defined by Archetypes
#'
#' @description
#' This function fits new data points to a convex hull of the given archetype
#' coordinates by minimizing the square distance of the projection.
#'
#' @param A Numeric matrix (K x M) of archetype coordinates.
#' @param X Numeric matrix (N x M) of new data points to fit.
#' @param eps Numeric scalar used for numerical stability; ensures non-negativity of fit.
#' @param method Character string specifying the fitting method:
#'   - `"nnls"`: Non-negative least squares with subsequent projection onto simplex.
#'   - `"QP"`: Quadratic programming approach to directly fit onto simplex.
#' @param project Function used to project the `nnls` fit onto the simplex  (default: `proj_l1`).
#' @param lambda Numeric scalar used for numerical stability in `QP` solver (default: 1e-8).
#'
#' @details
#' The function solves the constrained least-squares problem that minimizes
#' `||X - S*A||_2` such that `S >= 0` and `rowSums(S) = 1`.
#'
#' The method `"nnls"` fits `S` via non-negative least squares using the `nnls`
#' package. The results are then projected onto the simplex to ensure
#' row-stochasticity via the `project` method (by default \code{\link{proj_l1}}).
#'
#' The method `"QP"` solves the full constrained least squares problem via
#' quadratic programming using the `quadprog` package, enforcing both non-negativity
#' and row-stochasticity constraints directly. This method may be more accurate
#' but also slower for large data sets.
#'
#' `lambda` acts as a $L_2$ regularization parameter to ensure numerical stability
#' of the quadratic programming solver.
#'
#' @return Numeric row-stochastic matrix (N x K) of fitted compositions.
#'
#' @export
fit_simplex <- function(A, X, method = c("nnls", "QP"), eps = 0, project = proj_l1, lambda = 1e-8) {

    # Prepare Inputs
    X <- if (is.vector(X)) matrix(X, nrow = 1L) else as.matrix(X)
    A <- as.matrix(A)

    # Check inputs
    stopifnot("eps must be non-negative" = eps >= 0)
    stopifnot("lambda must be non-negative" = lambda >= 0)
    stopifnot("Incompatible dimensions between archetypes and new data" = ncol(A) == ncol(X))

    # Main
    method <- match.arg(method)
    S <- switch(
        method,
        nnls = project(fit_nnls(X, t(A)), eps = eps),  # TODO: add bigM?
        QP   = fit_qp(A, X, eps, proj_l1, lambda)
    )

    # Prepare output
    colnames(S) <- rownames(A)
    rownames(S) <- rownames(X)
    S
}

fit_qp <- function(A,
                   X,
                   eps,
                   project = NULL,
                   lambda = 1e-8,
                   row_sum_bounds = c(1, 1),
                   ...) {

    N <- nrow(X)
    K <- nrow(A)
    stopifnot("row_sum_bounds must be a length-two numeric vector" = length(row_sum_bounds) == 2L)
    stopifnot("row_sum_bounds must be non-negative and ordered" =
                  row_sum_bounds[1L] >= 0 && row_sum_bounds[1L] <= row_sum_bounds[2L])

    # The QP program formulation is:
    # min 0.5 s' Q s - d' s  with constraints t(Amat) %*% s >= bvec
    # The original objective is min_S ||X - S A||^2_F
    # This expands to min 2*S' A' A S- 2 A' X S + const
    # Q = 2 A'A,
    Dmat <- tcrossprod(A)
    diag(Dmat) <- diag(Dmat) + lambda # for numerical stability
    if (diff(row_sum_bounds) == 0) {
        # Equality 1' s = bound via meq = 1, plus inequalities s >= 0.
        meq  <- 1L
        Amat <- cbind(matrix(1, nrow = K, ncol = 1L), diag(K))
        bvec <- c(row_sum_bounds[1L], rep(0, K))
    } else {
        # Inequalities lower <= 1' s <= upper, plus s >= 0.
        meq  <- 0L
        Amat <- cbind(rep(1, K), rep(-1, K), diag(K))
        bvec <- c(row_sum_bounds[1L], -row_sum_bounds[2L], rep(0, K))
    }

    # Main Loop
    S <- matrix(eps, nrow = N, ncol = K)
    for (i in seq_len(N)) {
        dvec   <- A %*% X[i, ]
        S[i, ] <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
    }
    # Make sure to project onto simplex if we only enforced equality constraints (i.e., no upper bound)
    if (meq == 1L && !is.null(project))
        S <- project(S, eps) * row_sum_bounds[1L]
    S
}


#' One-hot encode a vector
#'
#' `onehot()` converts a vector into a **one-hot encoded** matrix: each distinct
#' value corresponds to a column, and each row has a `1` in the matching column
#' and `0` elsewhere. This is useful for turning categorical values into numbers
#' for modeling and data preprocessing.
#'
#' - **Numeric/integer input:** values are treated as 1-based column indices.
#'   The number of columns (`nc`) defaults to the largest value in the input,
#'   but you can set `nc` manually if you need more columns (e.g., reserving
#'   space for unseen categories).
#' - **Factor input:** each factor level becomes a column; the number of columns
#'   is fixed to the number of levels and the columns are named by those levels.
#' - **Character input:** automatically converted to a factor before encoding.
#'
#' Missing values produce rows of all zeros.
#'
#' @param ind A vector to encode. Can be numeric/integer, factor, or character.
#' @param sparse Logical; if `TRUE`, return a memory-efficient sparse matrix.
#' @param ... Additional arguments passed to the specific input handler.
#'
#' @return A one-hot encoded matrix:
#'   - dense `matrix` when `sparse = FALSE`
#'   - sparse `Matrix::dgCMatrix` when `sparse = TRUE`
#'
#' @examples
#' # Numeric vector: columns correspond to indices 1..max(ind)
#' onehot(c(1, 2, 1, 3))
#'
#' # Factor: columns follow factor levels (and are named)
#' f <- factor(c("red", "blue", "red", NA, "green"), levels = c("red","blue","green"))
#' onehot(f)
#'
#' # Character: automatically treated as a factor
#' onehot(c("cat","dog","cat"), sparse = TRUE)
#'
#' # Reserving extra columns for future/unseen categories (numeric input only)
#' onehot(c(1, 2, 1), sparse = FALSE, nc = 5)
#'
#' @export
onehot <- function(ind, sparse = FALSE, ...) {
    UseMethod("onehot")
}

#' @describeIn onehot Numeric indices
#' @param nc Number of columns (for numeric input only). If `NULL`, uses `max(ind)`.
#' @exportS3Method onehot default
onehot.default <- function(ind, sparse = FALSE, nc = NULL, ...) {
    nr <- length(ind)
    if (is.null(nc)) nc <- max(ind, na.rm = TRUE)
    dnames <- list(names(ind), NULL)

    keep <- which(!is.na(ind))
    if (length(keep) == 0) {
        # Edge Case: All input is NA return zero matrix
        if (sparse) {
            B <- Matrix::sparseMatrix(
                i = integer(0),
                j = integer(0),
                x = numeric(0),
                dims = c(nr, nc),
                dimnames = dnames
            )
        } else {
            B <- matrix(0, nrow = nr, ncol = nc, dimnames = dnames)
        }
        return(B)
    }

    # Main case
    if (sparse) {
        B <- Matrix::sparseMatrix(
            i = keep,
            j = ind[keep],
            x = 1,
            dims = c(nr, nc),
            dimnames = dnames
        )
    } else {
        B <- matrix(0, nrow = nr, ncol = nc, dimnames = dnames)
        B[cbind(keep, ind[keep])] <- 1
    }
    B
}

#' @describeIn onehot Factor method: `nc` is fixed to the number of levels
#' @exportS3Method onehot factor
onehot.factor <- function(ind, sparse = FALSE, ...) {
    lvl <- levels(ind)
    ind <- as.integer(ind)
    nc  <- length(lvl)
    out <- onehot.default(ind, sparse = sparse, nc = nc)
    dimnames(out)[[2]] <- lvl
    out
}

#' @describeIn onehot Character method: delegates to `factor`
#' @exportS3Method onehot character
onehot.character <- function(ind, sparse = FALSE, ...) {
    onehot.factor(as.factor(ind), sparse = sparse)
}
