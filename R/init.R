#' Archetypal Analysis Initialization Functions
#'
#' @description Various initialization methods for Archetypal Analysis (see Details).
#'
#' @param X a numeric matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes to be initialized
#' @param method initialization method. One of `"uniform_archetypes"`,
#'   `"furthest_first"`, `"kmeans_pp"`, `"furthest_sum"`, `"coreset_initfn"`,
#'   `"aa_pp"`, `"aa_pp_mc"`, or `"hull_outmost"` (default: `"furthest_sum"`).
#' @param sparse whether `B` should be a sparse matrix (default: same as X)
#' @param m optional batch size for coreset initialization
#' @param batch_size optional batch size for MCMC approximation of the AA++ initialization
#' @param hull_method strategy used by `"hull_outmost"`. One of `"full"`,
#'   `"projected"`, or `"partitioned"` (default: `"full"`).
#' @param projected_dim projection dimension used by `"hull_outmost"` when
#'   `hull_method = "projected"` (default: `2L`).
#' @param n_partitions number of row partitions used by `"hull_outmost"` when
#'   `hull_method = "partitioned"` (default: `10L`).
#' @param n_projection_max optional maximum number of random projections used by
#'   `"hull_outmost"` for `hull_method = "projected"`. If `NULL`, all
#'   feature projections are evaluated.
#' @param use_unique_candidates whether `"hull_outmost"` should de-duplicate
#'   hull candidates before vote tallying (default: `FALSE`).
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
#' - `"hull_outmost"`: computes hull candidates using one of the
#'   `hull_method` strategies (`"full"`, `"projected"`, or `"partitioned"`)
#'   and then selects `K` archetypes via an outmost-vote ranking. This family of
#'   hull-based initializations is adapted from the \pkg{archetypal} package
#'   (Mouselimis et al., 2025).
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
#' Mouselimis, L., et al. (2025).
#' \pkg{archetypal}: Archetypal Analysis with Principal Convex Hull Analysis.
#' R package version 1.3.1.
#' \url{https://cran.r-project.org/package=archetypal}
#'
#' @export
aa_init <- function(X,
                    K,
                    method = "furthest_sum",
                    sparse = inherits(X, "sparseMatrix"),
                    m = NULL,
                    batch_size = NULL,
                    hull_method = c("full", "projected", "partitioned"),
                    projected_dim = 2L,
                    n_partitions = 10L,
                    n_projection_max = NULL,
                    use_unique_candidates = FALSE,
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
            "aa_pp_mc",
            "hull_outmost"
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
    } else if (method == "hull_outmost") {
        hull_method <- match.arg(hull_method, c("full", "projected", "partitioned"))

        stopifnot("`projected_dim` must be an integer" = projected_dim == as.integer(projected_dim))
        stopifnot("`projected_dim` must be in [1, ncol(X)]" =
                      projected_dim >= 1L && projected_dim <= ncol(X))

        stopifnot("`n_partitions` must be an integer" = n_partitions == as.integer(n_partitions))
        stopifnot("`n_partitions` must be greater than zero" = n_partitions >= 1L)

        stopifnot("`use_unique_candidates` must be a single logical" =
                      is.logical(use_unique_candidates) && length(use_unique_candidates) == 1L)
        stopifnot("`use_unique_candidates` must not be NA" = !is.na(use_unique_candidates))

        if (!is.null(n_projection_max)) {
            stopifnot("`n_projection_max` must be an integer" =
                          n_projection_max == as.integer(n_projection_max))
            stopifnot("`n_projection_max` must be greater than zero" = n_projection_max >= 1L)
        }
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
        hull_outmost       = hull_outmost(
            X,
            K,
            hull_method = hull_method,
            projected_dim = projected_dim,
            n_partitions = n_partitions,
            n_projection_max = n_projection_max,
            use_unique_candidates = use_unique_candidates,
            ...
        ),
        coreset_initfn     = coreset_initfn(X, K, m = m, ...)
    )
    .ind_to_init(X, ind, sparse = sparse)
}

uniform_archetypes <- function(X, K, ...) sample(nrow(X), K, replace = FALSE)

furthest_first <- function(X, K, distances = NULL, center_dists = NULL, ...) {

    b <- integer(K)  # indices of archetypes

    # 1) randomly select the first archetype
    dists <- center_dists %||% .dist2(X, center = TRUE)
    b[1L] <- .sample_distal_points(dists, 1L)

    # 2) compute next K-1 archetypes by selecting the furthest from current set
    for (k in seq_len(K - 1L)) {
        dists <- .dist_to_nearest_archetype(X, b[1:k], distances = distances)
        b[k + 1L] <- which.max(dists)
    }

    b
}


kmeans_pp <- function(X, K, sparse = inherits(X, "sparseMatrix"),
                      distances = NULL, center_dists = NULL, ...) {

    b <- integer(K)  # indices of archetypes

    # 1) randomly select the first archetype
    dists <- center_dists %||% .dist2(X, center = TRUE)
    b[1L] <- .sample_distal_points(dists, 1L)

    # 2) compute next K-1 archetypes by sampling from the points furthest from the current set
    for (k in seq_len(K - 1)) {
        dists <- .dist_to_nearest_archetype(X, b[1:k], distances = distances)
        b[k + 1L] <- .sample_distal_points(dists, 1L)
    }

    b
}

furthest_sum <- function(X, K, distances = NULL, refinement_steps = 10L, ...) {
    stopifnot("`refinement_steps` must be a single non-negative integer" =
                  length(refinement_steps) == 1L &&
                      !is.na(refinement_steps) &&
                      refinement_steps == as.integer(refinement_steps) &&
                      refinement_steps >= 0L)

    N <- nrow(X)
    refinement_steps <- as.integer(refinement_steps)
    effective_refinement_steps <- min(refinement_steps, max(N - K, 0L))
    n_iterations <- K + effective_refinement_steps

    # 1) randomly select the first archetype
    ind_t <- sample(N, 1L)
    b <- ind_t  # queue of active archetype indices
    sum_dist <- numeric(N)
    eligible <- rep(TRUE, N)
    eligible[ind_t] <- FALSE

    get_dists <- function(ind) {
        if (is.null(distances))
            return(.dist2(X, X[ind, , drop = FALSE]))
        distances[, ind]
    }

    # 2) greedily add candidates by cumulative distance, then refine by
    # removing the oldest active candidate before selecting the next one.
    for (k in seq_len(n_iterations)) {
        if (length(b) >= K) {
            drop_ind <- b[1L]
            sum_dist <- sum_dist - get_dists(drop_ind)
            eligible[drop_ind] <- TRUE
            b <- b[-1L]
        }

        sum_dist <- sum_dist + get_dists(ind_t)
        candidate_rows <- which(eligible)
        ind_t <- candidate_rows[which.max(sum_dist[candidate_rows])]

        b <- c(b, ind_t)
        eligible[ind_t] <- FALSE
    }

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
    eps <- ifelse(sparse, 1e-8, 0)
    for (k in 3:K) {
        A <- X[b[1:(k - 1)], , drop = FALSE]  # current archetypes
        S <- proj_l1(fit_nnls(X, t(A)), eps = eps)
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
    eps <- ifelse(inherits(X, "sparseMatrix"), 1e-8, 0)
    for (k in 3:K) {
        A <- X[b[1:(k - 1)], , drop = FALSE]  # current archetypes

        batch <- sample(nrow(X), batch_size, replace = TRUE)
        S <- proj_l1(fit_nnls(X[batch, , drop = FALSE], t(A)), eps = eps)
        res <- X[batch, , drop = FALSE] - S %*% A
        dists <- rowSums(res * res)  # squared residuals
        dists[!is.finite(dists)] <- 0

        ib <- 1L
        for (j in seq_along(dists))
            if (dists[j] > stats::runif(1) * dists[ib])
                ib <- j
        b[k] <- batch[ib]
    }
    b
}

hull_outmost <- function(X,
                         K,
                         hull_method = c("full", "projected", "partitioned"),
                         projected_dim = 2L,
                         n_partitions = 10L,
                         n_projection_max = NULL,
                         use_unique_candidates = FALSE,
                         ...) {
    hull_method <- match.arg(hull_method)

    # Choose how to construct the candidate hull set; selection logic is shared below.
    # Returns a vector of row indices of X that are candidates for selection as archetypes.
    candidate_rows <- switch(
        hull_method,
        full = .aa_build_hull_candidates_full(X),
        projected = .aa_build_hull_candidates_projected(
            X,
            projected_dim = projected_dim,
            n_projection_max = n_projection_max
        ),
        partitioned = .aa_build_hull_candidates_partitioned(X, n_partitions = n_partitions)
    )

    tallies <- .aa_tally_outmost_votes(
        X,
        candidate_rows,
        use_unique_candidates = use_unique_candidates
    )
    .aa_select_from_votes(
        vote_order = tallies[["ranking"]],
        K = K,
        fallback_pool = tallies[["fallback_pool"]],
        N = nrow(X)
    )
}

.aa_build_hull_candidates_full <- function(X) {
    N <- nrow(X)
    D <- ncol(X)
    X_dense <- as.matrix(X)

    # Full-hull strategy: use the geometric envelope in the original feature space.
    if (D == 1L) {
        return(as.integer(unique(c(which.min(X_dense[, 1L]), which.max(X_dense[, 1L])))))
    }

    if (D == 2L) {
        return(as.integer(unique(grDevices::chull(X_dense))))
    }

    if (!requireNamespace("geometry", quietly = TRUE)) {
        stop(
            "`hull_outmost` with dimensions > 2 requires package `geometry`. ",
            "Install `geometry` or use `hull_method = 'projected'` with `projected_dim <= 2`."
        )
    }

    # facet: a simplex of D vertices that defines a face of the convex hull.
    # The vertices are returned as row indices of X.
    facets <- tryCatch(
        geometry::convhulln(X_dense, options = "Fx"),
        error = function(e) NULL
    )
    if (is.null(facets))
        return(seq_len(N))

    if (is.null(dim(facets))) {
        idx <- as.integer(facets)
    } else {
        idx <- as.integer(unique(as.vector(facets)))
    }
    idx <- idx[idx >= 1L & idx <= N]
    as.integer(unique(idx))
}

.aa_build_hull_candidates_projected <- function(X,
                                                projected_dim = 2L,
                                                n_projection_max = NULL) {
    D <- ncol(X)
    projected_dim <- as.integer(projected_dim)
    stopifnot("`projected_dim` must be in [1, ncol(X)]" = projected_dim >= 1L && projected_dim <= D)

    if (projected_dim == D)
        return(.aa_build_hull_candidates_full(X))

    # Projected-hull strategy: aggregate envelope points across low-dimensional views.
    # Views are constructed by projecting (selecting) on random subsets of features.
    combs <- utils::combn(D, projected_dim)
    n_combs <- ncol(combs)
    selected <- seq_len(n_combs)
    if (!is.null(n_projection_max) && n_projection_max < n_combs)
        selected <- sample.int(n_combs, size = n_projection_max, replace = FALSE)

    rows <- integer(0)
    for (j in selected) {
        cols <- combs[, j]
        local <- .aa_build_hull_candidates_full(X[, cols, drop = FALSE])
        rows <- c(rows, local)
    }
    rows
}

.aa_build_hull_candidates_partitioned <- function(X, n_partitions = 10L) {
    N <- nrow(X)
    n_partitions <- min(as.integer(n_partitions), N)

    # Partitioned-hull strategy: approximate global extremes from local hull champions.
    part_id <- rep(seq_len(n_partitions), length.out = N) # 1, 2, ..., n_partitions, 1, 2, ...
    shuffled_rows <- sample.int(N, N, replace = FALSE)
    rows <- integer(0)

    for (p in seq_len(n_partitions)) {
        grp <- shuffled_rows[part_id == p]
        if (length(grp) == 0L)
            next
        if (length(grp) <= 2L) {
            rows <- c(rows, grp)
            next
        }

        local <- .aa_build_hull_candidates_full(X[grp, , drop = FALSE])
        rows <- c(rows, grp[local])
    }

    rows
}

.aa_tally_outmost_votes <- function(X,
                                    candidate_rows,
                                    use_unique_candidates = FALSE) {
    N <- nrow(X)
    candidate_rows <- as.integer(candidate_rows)
    candidate_rows <- candidate_rows[candidate_rows >= 1L & candidate_rows <= N]
    if (length(candidate_rows) == 0L)
        stop("No hull candidates available for `hull_outmost`.")

    if (use_unique_candidates)
        candidate_rows <- unique(candidate_rows)

    Xc <- as.matrix(X[candidate_rows, , drop = FALSE])
    if (nrow(Xc) == 1L) {
        ranking <- candidate_rows
        return(list(ranking = ranking, fallback_pool = unique(candidate_rows)))
    }

    dm <- as.matrix(stats::dist(Xc, method = "euclidean", diag = FALSE, upper = FALSE))
    voted_rows <- candidate_rows[max.col(dm, ties.method = "first")]

    fallback_pool <- sort(unique(candidate_rows))
    freq_table <- table(factor(voted_rows, levels = fallback_pool))
    freq <- as.integer(freq_table)

    ranking <- fallback_pool[order(-freq, fallback_pool)]
    list(ranking = ranking, fallback_pool = fallback_pool)
}

.aa_select_from_votes <- function(vote_order, K, fallback_pool, N) {
    # Select the K most frequently voted rows, breaking ties by row index order.
    out <- as.integer(utils::head(unique(vote_order), K))
    # pad with rows from the fallback pool (also ordered by row index)
    if (length(out) < K) {
        pad <- setdiff(as.integer(fallback_pool), out)
        out <- c(out, utils::head(pad, K - length(out)))
    }
    # pad with any remaining rows in order of row index.
    if (length(out) < K) {
        pad <- setdiff(seq_len(N), out)
        out <- c(out, utils::head(pad, K - length(out)))
    }
    as.integer(out)
}

# Compute the distance to the nearest archetype for each sample
# X is a matrix of samples
# ind is a vector of indices selecting the archetypes from X
.dist_to_nearest_archetype <- function(X, ind, distances = NULL) {
    if (!is.null(distances))
        return(matrixStats::rowMins(distances[, ind, drop = FALSE]))
    A <- X[ind, , drop = FALSE]  # archetypes
    dists <- .pdist2(A, X)
    matrixStats::colMins(dists)
}

# Sample points proportionally to their distance from a reference point
.sample_distal_points <- function(dists, size = 1) {
    N <- length(dists)
    p <- as.numeric(dists)
    p[!is.finite(p)] <- 0
    p[p < 0] <- 0
    if (sum(p) <= 0)
        return(sample(N, size = size, replace = TRUE))

    sample(N, size = size, replace = TRUE, prob = p)
}

.aa_normalize_row_indices <- function(ind, n, row_names = NULL) {
    stopifnot(mode(ind) %in% c("numeric", "logical", "character"))
    if (mode(ind) == "logical") {
        stopifnot("Logical indices length must be equal to number of rows" = length(ind) == n)
        ind <- which(ind)
    } else if (mode(ind) == "character") {
        ind <- match(ind, row_names, nomatch = 0L)
        stopifnot("Some row names do not match" = all(ind > 0L))
    } else if (all(ind <= 0)) {
        ind <- setdiff(seq_len(n), -ind)
    }
    ind <- ind[ind > 0]
    stopifnot("indices must be valid sample rows" = all(ind <= n))
    ind
}

# Initialize variables for Archetypal Analysis
.ind_to_init <- function(X, ind, sparse) {
    # make sure ind is positive indices selecting rows
    ind <- .aa_normalize_row_indices(ind, nrow(X), rownames(X))
    nm <- names(ind) %||% paste0("A", seq_along(ind))

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

# Internal random initialization for directional AA.
#
# This mirrors the reference directional AA MATLAB code, which initializes both
# the archetype generator C and the composition matrix S with exponential random
# entries and then L1-normalizes them. It is intentionally not exposed as an
# `aa_init()` method yet: Euclidean AA initializers usually select or construct
# concrete data-space archetypes, while this helper initializes the directional
# generator coefficients directly.
.aa_directional_random_init <- function(X_flip, K, eps = 0) {
    N <- nrow(X_flip)

    # B is paper C^T in YAAAP's row-oriented convention. Each row is a convex
    # combination over samples and therefore generates one directional archetype.
    B <- matrix(stats::rexp(K * N), nrow = K, ncol = N)
    B <- proj_l1(B, eps = eps)

    # S is paper S^T. Each row is the composition of one sample in the
    # directional archetype hull.
    S <- matrix(stats::rexp(N * K), nrow = N, ncol = K)
    S <- proj_l1(S, eps = eps)

    A <- B %*% X_flip
    list(A0 = A, B = B, S = S)
}
