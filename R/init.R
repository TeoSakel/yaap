#' Archetypal Analysis Initialization Functions
#'
#' @description Various initialization methods for Archetypal Analysis (see Details).
#'
#' @param X a numeric matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes to be initialized
#' @param method initialization method. One of `"uniform_archetypes"`,
#'   `"furthest_first"`, `"kmeans_pp"`, `"furthest_sum"`, `"coreset_initfn"`,
#'   `"aa_pp"`, or `"aa_pp_mc"` (default: `"furthest_sum"`).
#' @param sparse whether `B` should be a sparse matrix (default: same as X)
#' @param m optional batch size for coreset initialization
#' @param batch_size optional batch size for MCMC approximation of the AA++ initialization
#' @param ... additional arguments (used for compatibility with user-defined functions)
#'
#' @return a named list containing two matrices: the archetype coordinates `A`
#' and a row-stochastic matrix `B` such that `A = B %*% X`.
#'
#' @details
#'
#' `aa_init()` runs the selected initialization method and formats the result
#' as a named list containing the archetype coordinates `A` and row-stochastic
#' matrix `B`. User-defined initialization functions passed directly to
#' `archetypes_nnls()` or `archetypes_pgd()` must follow that same interface.
#'
#' The following initialization methods are available:
#'
#' - `"uniform_archetypes"`: selects uniformly at random `K` archetypes from `X`.
#' - `"furthest_first"`: selects the first archetype randomly and then greedily
#'   selects the point furthest from the current set of archetypes.
#' - `"kmeans_pp"`: a soft version of `furthest_first` where points are sampled
#'   in proportion to their distance from the current set of archetypes instead
#'   of greedily picking the furthest point every time.
#' - `"furthest_sum"`: selects the first archetype randomly and then greedily
#'   selects the next archetypes that maximizes the sum of distances of all
#'   points from the current set of archetypes (see Mørup & Hansen 2012).
#' - `"coreset_initfn"`: initializes archetypes using a coreset of size `m` sampled
#'   from the data matrix `X` and then applies `furthest_sum` on the coreset
#'   (see Mair & Brefeld 2019).
#' - `"aa_pp"`: AA++ is a probabilistic initialization method similar to `kmeans++`
#'   but instead using the distances to the current set of archetypes it uses
#'   distances to their convex hull (see Mair & Sjölund 2023).
#' - `"aa_pp_mc"`: MCMC approximation of AA++ initialization. Each archetype is
#'   sampled by performing AA++ on a sub-sample of size `batch_size` from the
#'   data matrix `X` (see Mair & Sjölund 2023).
#'
#' @references
#' Mørup, M., & Hansen, L. K. (2012).
#' Archetypal analysis for machine learning and data mining.
#' Neurocomputing, 80, 54-63.
#' \url{https://doi.org/10.1016/j.neucom.2011.06.033}
#'
#' Mair, S., & Sjölund, J. (2023).
#' Archetypal analysis++: Rethinking the initialization strategy.
#' \url{https://arxiv.org/abs/2301.13748}
#'
#' Mair, S., & Brefeld, U. (2019).
#' Coresets for archetypal analysis.
#' Advances in Neural Information Processing Systems, 32.
#' \url{https://proceedings.neurips.cc/paper_files/paper/2019/file/7f278ad602c7f47aa76d1bfc90f20263-Paper.pdf}
#'
#' @export
aa_init <- function(X,
                    K,
                    method = "furthest_sum",
                    sparse = inherits(X, "sparseMatrix"),
                    m = NULL,
                    batch_size = NULL,
                    ...) {

    # Input checks -----------------------------------------------------------
    stopifnot("K must be greater than 1" = K >= 1)
    stopifnot("K must be an integer" = K == as.integer(K))
    stopifnot("Number of samples must be at least K" = nrow(X) >= K)


    stopifnot("`method` must be a single string" = is.character(method) && length(method) == 1L)
    method <- match.arg(
        method,
        c(
            "uniform_archetypes",
            "furthest_first",
            "kmeans_pp",
            "furthest_sum",
            "coreset_initfn",
            "aa_pp",
            "aa_pp_mc"
        )
    )

    if (method == "coreset_initfn") {
        stopifnot("`m` must be supplied for `coreset_initfn`" = !is.null(m))
        stopifnot("`m` must be an integer" = m == as.integer(m))
        stopifnot("`m` must be at least K" = m >= K)
        stopifnot("Number of samples must be at least `m`" = nrow(X) >= m)
    } else if (method == "aa_pp_mc") {
        batch_size <- ifelse(is.null(batch_size), m, batch_size)
        stopifnot("`batch_size` must be supplied for `aa_pp_mc`" = !is.null(batch_size))
        stopifnot("`batch_size` must be an integer" = batch_size == as.integer(batch_size))
        stopifnot("`batch_size` must be at least K" = batch_size >= K)
        stopifnot("Number of samples must be at least `batch_size`" = nrow(X) >= batch_size)
    }

    # Main code ----------------------------------------------------------------

    # Edge case for K=1: return the point closest to the mean (the "archemean")
    if (K == 1) {
        ind <- which.min(.dist2(X, center = TRUE))
        return(.ind_to_init(X, ind, sparse = sparse))
    }

    ind <- switch(
        method,
        uniform_archetypes = uniform_archetypes(X, K, ...),
        furthest_first     = furthest_first(X, K, ...),
        kmeans_pp          = kmeans_pp(X, K, ...),
        furthest_sum       = furthest_sum(X, K, ...),
        aa_pp              = aa_pp(X, K, ...),
        aa_pp_mc           = aa_pp_mc(X, K, batch_size = batch_size, ...),
        coreset_initfn     = coreset_initfn(X, K, m = m, ...)
    )
    .ind_to_init(X, ind, sparse = sparse)
}

uniform_archetypes <- function(X, K, ...) sample(nrow(X), K, replace = FALSE)

furthest_first <- function(X, K, ...) {

    b <- integer(K)  # indices of archetypes

    # 1) randomly select the first archetype
    dists <- .dist2(X, center = TRUE)  # distances from the mean
    b[1L] <- .sample_distal_points(dists, 1L)

    # 2) compute next K-1 archetypes by selecting the furthest from current set
    for (k in seq_len(K - 1L)) {
        dists <- .dist_to_nearest_archetype(X, b[1:k])
        b[k + 1L] <- which.max(dists)
    }

    b
}


kmeans_pp <- function(X, K, sparse = inherits(X, "sparseMatrix"), ...) {

    b <- integer(K)  # indices of archetypes

    # 1) randomly select the first archetype
    dists <- .dist2(X, center = TRUE)  # distances from the mean
    b[1L] <- .sample_distal_points(dists, 1L)

    # 2) compute next K-1 archetypes by sampling from the points furthest from the current set
    for (k in seq_len(K - 1)) {
        dists <- .dist_to_nearest_archetype(X, b[1:k])
        b[k + 1L] <- .sample_distal_points(dists, 1L)
    }

    b
}

furthest_sum <- function(X, K, ...) {

    b <- integer(K)  # indices of archetypes

    # 1) randomly select the first archetype
    b[1L] <- sample(nrow(X), 1L)

    # 2) compute initial distances from that first point
    dists <- .dist2(X, X[b[1L], ])
    initial_dists <- dists

    # 3) select k−1 points by adding up distances
    select_max <- function(dists, archetypes) {
        dists[archetypes] <- 0  # current archetypes cannot be selected again
        which.max(dists)
    }

    for (k in seq_len(K - 1)) {
        b[k + 1L] <- select_max(dists, b[1:k])
        dists <- dists + .dist2(X, X[b[k + 1L], , drop = FALSE])
    }

    # 4) “forget” the very first random pick and select new first archetype
    dists <- dists - initial_dists
    dists[b[-1]] <- 0  # do not select any archetype from from 2:K
    b[1] <- which.max(dists)
    b
}

coreset_initfn <- function(X, K, m, ...) {
    # Coresets for Archetypal Analysis - Mair and Brefeld, 2019
    # https://github.com/smair/archetypalanalysis-coreset/blob/master/code/coresets.py
    # https://github.com/smair/archetypalanalysis-coreset/blob/master/code/experiments.py
    # m: cardinality of coreset

    q <- .dist2(X, center = TRUE)  # distances from the mean
    coreset <- .sample_distal_points(q, m)
    b <- furthest_sum(X[coreset, , drop = FALSE], K, ...)
    coreset[b]
}

aa_pp <- function(X, K, sparse = inherits(X, "sparseMatrix"), ...) {
    # A++ initialization for Archetypal Analysis - Mair and Brefeld, 2019

    # if K is 2 AA++ reduces to kmeans++
    if (K == 2) return(kmeans_pp(X, K, ...))

    b <- integer(K)  # indices of archetypes

    # 1) randomly select the first archetype
    b[1L] <- sample(nrow(X), 1L)

    # 2) sample the second archetype from points distal to the first
    b[2L] <- .sample_distal_points(.dist2(X, X[b[1], ]), 1L)

    # 3) compute the first K-2 archetypes by iteratively running AA and sampling
    # from the points distal to the current archetype convex-hull.
    for (k in 3:K) {
        A <- X[b[1:(k - 1)], , drop = FALSE]  # current archetypes
        S <- fit_nnls(X, t(A))
        res <- X - S %*% A
        dists <- rowSums(res * res)  # squared residuals
        b[k] <- .sample_distal_points(dists, 1L)
    }
    b
}

aa_pp_mc <- function(X, K, batch_size = m, m = NULL, ...) {
    # AA++ initialization for Archetypal Analysis - Mair and Brefeld, 2019
    # approximate AA++ initialization by sampling m points each time

    # if K is 2 AA++ reduces to kmeans++
    if (K == 2) return(kmeans_pp(X, K, ...))

    # if n equal m it reduces to AA++
    if (nrow(X) == batch_size) return(aa_pp(X, K, ...))

    # 1) randomly select the first archetype
    b <- integer(K)  # indices of archetypes
    b[1] <- sample(nrow(X), 1L)

    # 2) sample the second archetype from points distal to the first
    b[2] <- .sample_distal_points(.dist2(X, X[b[1], ]), 1L)

    # 3) compute the first K-2 archetypes by iteratively running AA.
    for (k in 3:K) {
        A <- X[b[1:(k - 1)], , drop = FALSE]  # current archetypes

        batch <- sample(nrow(X), batch_size, replace = TRUE)
        S <- fit_nnls(X[batch, , drop = FALSE], t(A))
        res <- X[batch, , drop = FALSE] - S %*% A
        dists <- rowSums(res * res)  # squared residuals

        ib <- 1L
        for (j in seq_along(dists))
            if (dists[j] > stats::runif(1) * dists[ib])
                ib <- j
        b[k] <- batch[ib]
    }
    b
}

# Compute the distance to the nearest archetype for each sample
# X is a matrix of samples
# ind is a vector of indices selecting the archetypes from X
.dist_to_nearest_archetype <- function(X, ind) {
    A <- X[ind, , drop = FALSE]  # archetypes
    dists <- .pdist2(A, X)
    matrixStats::colMins(dists)
}

# Sample points proportionally to their distance from a reference point
.sample_distal_points <- function(dists, size = 1) {
    N <- length(dists)
    sample(N, size = size, replace = TRUE, prob = dists)
}

# Initialize variables for Archetypal Analysis
.ind_to_init <- function(X, ind, sparse) {
    # make sure ind is positive indices selecting rows
    stopifnot(mode(ind) %in% c("numeric", "logical", "character"))
    if (mode(ind) == "logical") {
        stopifnot(
            "Logical indices length must be equal to number of rows in X" = length(ind) == nrow(X)
        )
        ind <- which(ind)  # convert logical to indices
    } else if (mode(ind) == "character") {
        ind <- match(ind, rownames(X), nomatch = 0L)
        stopifnot("Some archetype names do not match rows in X" = all(ind > 0L))
    } else if (all(ind <= 0)) {
        ind <- setdiff(seq_len(nrow(X)), -ind)
    }

    ind <- ind[ind > 0]
    nm <- if (!is.null(names(ind))) {
        names(ind)
    } else if (!is.null(rownames(X))) {
        rownames(X)[ind]
    } else {
        paste0("A", seq_along(ind))
    }

    A <- X[ind, , drop = FALSE]  # Archetypes
    B <- onehot(ind, sparse = sparse, nrow(X))  # Row-stochastic matrix
    rownames(A) <- rownames(B) <- nm
    list(A = A, B = B)
}

.init_S <- function(X, A, eps = 0) {
    S <- proj_l1(1 / .pdist2(X, A), eps = eps)  # init S by similarity score
    S[is.nan(S)] <- 1 # NaNs = Inf/Inf for the archetypes
    S
}
