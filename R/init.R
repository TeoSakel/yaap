#' Archetypal Analysis Initialization Functions
#'
#' Various initialization methods for Archetypal Analysis (see Details).
#'
#' @param X a numeric matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes to be initialized
#' @param method initialization method. One of `"random"`, `"dirichlet"`,
#'   `"furthest_first"`, `"kmeans_pp"`, `"furthest_sum"`, `"aa_pp"`, or
#'   `"hull_outmost"` (default: `"furthest_sum"`).
#' @param sparse whether `B` should be a sparse matrix (default: same as X)
#' @param batch_size optional number of candidate rows to sample.
#' @param batch_type candidate sampling strategy, `"distal"` or `"uniform"`.
#' @param batch_replace whether to sample candidate batches with replacement.
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
#' - `"random"`: selects uniformly at random `K` rows of `X` as archetypes.
#' - `"dirichlet"`: draws each archetype as a random convex combination of all
#'   data points by sampling `B` row-wise from a Dirichlet(alpha, ..., alpha)
#'   distribution (Gamma(alpha, 1) weights normalized to unit sum). `alpha = 1`
#'   (default) gives a uniform distribution over all data points; `alpha < 1`
#'   concentrates mass near `K` data points, yielding sparser rows; `alpha > 1`
#'   concentrates mass toward the centroid. Accepts optional `alpha` argument
#'   via `...`.
#' - `"furthest_first"`: selects the first archetype randomly and then greedily
#'   selects the point furthest from the current set of archetypes.
#' - `"kmeans_pp"`: a soft version of `furthest_first` where points are sampled
#'   in proportion to their distance from the current set of archetypes instead
#'   of greedily picking the furthest point every time.
#' - `"furthest_sum"`: selects the first archetype randomly and then greedily
#'   selects the next archetypes that maximizes the sum of distances of all
#'   points from the current set of archetypes (see Mørup & Hansen 2012).
#' - `"aa_pp"`: AA++ is a probabilistic initialization method similar to `kmeans_pp`
#'   but instead using the distances to the current set of archetypes it uses
#'   distances to their convex hull (see Mair & Sjölund 2023).
#' - `"hull_outmost"`: computes hull candidates using one of the
#'   `hull_method` strategies (`"full"`, `"projected"`, or `"partitioned"`)
#'   and then selects `K` archetypes via an outmost-vote ranking. This family of
#'   hull-based initializations is adapted from the \pkg{archetypal} package
#'   (Mouselimis et al., 2025).
#'
#' Supplying `batch_size` applies the selected method to a candidate batch
#' rather than to all rows. For all methods except `"aa_pp"`, one batch is
#' sampled up front. For `"aa_pp"`, a fresh batch is sampled at each
#' approximation step; this is a variant of the Monte Carlo AA++ approximation
#' rather than the exact `"aa_pp_mc"` scheme. `batch_replace = FALSE` by default.
#'
#' The default `batch_type = "distal"` is the coreset sampling strategy of
#' Mair & Brefeld (2019): rows farther from the data center are more likely to
#' be candidates, which saves memory and time while focusing on the boundary
#' where archetypes are expected to lie. Use `batch_type = "uniform"` for an
#' unbiased candidate batch.
#'
#' @references
#' Mørup, M., & Hansen, L. K. (2012).
#' Archetypal analysis for machine learning and data mining.
#' Neurocomputing, 80, 54-63.
#' \doi{10.1016/j.neucom.2011.06.033}
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
                    batch_size = NULL,
                    batch_type = c("distal", "uniform"),
                    batch_replace = FALSE,
                    hull_method = c("full", "projected", "partitioned"),
                    projected_dim = 2L,
                    n_partitions = 10L,
                    n_projection_max = NULL,
                    use_unique_candidates = FALSE,
                    ...) {
    # Input checks -----------------------------------------------------------
    if (!is_count(K)) {
        stop("`K` must be a positive integer.", call. = FALSE)
    }
    if (nrow(X) < K) {
        stop("Number of samples must be at least `K`.", call. = FALSE)
    }


    if (!is_single_string(method)) {
        stop("`method` must be a single string.", call. = FALSE)
    }
    method <- match.arg(
        method,
        c(
            "random",
            "dirichlet",
            "furthest_first",
            "kmeans_pp",
            "furthest_sum",
            "aa_pp",
            "hull_outmost"
        )
    )
    batch_type <- match.arg(batch_type)
    if (!is_logical(batch_replace)) {
        stop("`batch_replace` must be TRUE or FALSE.", call. = FALSE)
    }
    batch_size <- .aa_validate_batch_size(
        batch_size,
        n = nrow(X),
        K = K,
        replace = batch_replace
    )

    if (method == "hull_outmost") {
        hull_method <- match.arg(hull_method, c("full", "projected", "partitioned"))

        if (!(is_count(projected_dim) && projected_dim <= ncol(X))) {
            stop("`projected_dim` must be between 1 and `ncol(X)`.", call. = FALSE)
        }
        if (!is_count(n_partitions)) {
            stop("`n_partitions` must be a positive integer.", call. = FALSE)
        }

        if (!is_logical(use_unique_candidates)) {
            stop("`use_unique_candidates` must be TRUE or FALSE.", call. = FALSE)
        }

        if (!is.null(n_projection_max)) {
            if (!is_count(n_projection_max)) {
                stop("`n_projection_max` must be a positive integer.", call. = FALSE)
            }
        }
    }

    # Main code ----------------------------------------------------------------

    # Edge case for K=1: return the point closest to the mean (the "archemean")
    if (K == 1) {
        ind <- which.min(.aa_dist2(X, center = TRUE))
        return(.aa_ind_to_init(X, ind, sparse = sparse))
    }

    batch <- .aa_sample(
        X,
        size = batch_size,
        type = batch_type,
        replace = batch_replace
    )
    X_init <- X[batch, , drop = FALSE]

    result <- switch(method,
        random = uniform_archetypes(X_init, K, ...),
        dirichlet = dirichlet(X_init, K, ...),
        furthest_first = furthest_first(X_init, K, ...),
        kmeans_pp = kmeans_pp(X_init, K, ...),
        furthest_sum = furthest_sum(X_init, K, ...),
        aa_pp = aa_pp(
            X,
            K,
            batch_size = batch_size,
            batch_type = batch_type,
            batch_replace = batch_replace,
            ...
        ),
        hull_outmost = hull_outmost(
            X_init,
            K,
            hull_method = hull_method,
            projected_dim = projected_dim,
            n_partitions = n_partitions,
            n_projection_max = n_projection_max,
            use_unique_candidates = use_unique_candidates,
            ...
        )
    )
    if (is.list(result)) {
        if (is.null(batch_size)) {
            return(result)
        }
        return(.aa_expand_init_batch(result, X, batch, sparse = sparse))
    }
    if (method == "aa_pp") {
        return(.aa_ind_to_init(X, result, sparse = sparse))
    }
    .aa_ind_to_init(X, batch[result], sparse = sparse)
}

.aa_validate_batch_size <- function(batch_size, n, K, replace = FALSE) {
    if (is.null(batch_size)) {
        return(NULL)
    }
    if (!is_count(batch_size, start_from = K)) {
        stop("`batch_size` must be a positive integer at least `K`.", call. = FALSE)
    }
    if (!replace) {
        if (batch_size > n) {
            stop(
                paste(
                    "`batch_size` must be no larger than the number of samples",
                    "when `batch_replace = FALSE`."
                ),
                call. = FALSE
            )
        }
    }
    as.integer(batch_size)
}

.aa_sample <- function(x,
                       size = NULL,
                       type = "distal",
                       replace = FALSE) {
    nr <- nrow(x)
    has_rows <- !is.null(nr)
    n <- if (has_rows) nr else length(x)
    if (is.null(size)) {
        return(seq_len(n))
    }
    prob <- if (!has_rows && identical(type, "distal")) {
        x
    } else if (identical(type, "distal")) {
        .aa_dist2(x, center = TRUE)
    } else {
        NULL
    }
    sample(n, size = size, replace = replace, prob = prob)
}

.aa_expand_init_batch <- function(init, X, batch, sparse = FALSE) {
    B_batch <- init[["B"]]
    B <- matrix(0, nrow = nrow(B_batch), ncol = nrow(X))
    B[, batch] <- B_batch
    rownames(B) <- rownames(B_batch)
    colnames(B) <- rownames(X)
    if (sparse) {
        B <- as(B, "sparseMatrix")
    }
    A <- as.matrix(B %*% X)
    rownames(A) <- rownames(B)
    list(A = A, B = B)
}

uniform_archetypes <- function(X, K, ...) sample(nrow(X), K, replace = FALSE)

furthest_first <- function(X, K, distances = NULL, center_dists = NULL, ...) {
    b <- integer(K) # indices of archetypes

    # 1) randomly select the first archetype
    dists <- center_dists %||% .aa_dist2(X, center = TRUE)
    b[1L] <- .aa_sample(dists, size = 1L, replace = TRUE)

    # 2) compute next K-1 archetypes by selecting the furthest from current set
    # TODO: add a refinement loop similar to furthest_sum
    for (k in seq_len(K - 1L)) {
        dists <- .aa_dist_to_nearest_archetype(X, b[1:k], distances = distances)
        b[k + 1L] <- which.max(dists)
    }

    b
}


kmeans_pp <- function(X, K, sparse = inherits(X, "sparseMatrix"),
                      distances = NULL, center_dists = NULL, ...) {
    b <- integer(K) # indices of archetypes

    # 1) randomly select the first archetype
    dists <- center_dists %||% .aa_dist2(X, center = TRUE)
    b[1L] <- .aa_sample(dists, size = 1L, replace = TRUE)

    # 2) compute next K-1 archetypes by sampling from the points furthest from the current set
    for (k in seq_len(K - 1)) {
        dists <- .aa_dist_to_nearest_archetype(X, b[1:k], distances = distances)
        b[k + 1L] <- .aa_sample(dists, size = 1L, replace = TRUE)
    }

    b
}

furthest_sum <- function(X, K, distances = NULL, refinement_steps = 10L, ...) {
    if (!is_count(refinement_steps, start_from = 0L)) {
        stop("`refinement_steps` must be a single non-negative integer.", call. = FALSE)
    }

    N <- nrow(X)
    effective_refinement_steps <- min(refinement_steps, max(N - K, 0L))
    n_iterations <- K + effective_refinement_steps

    # 1) randomly select the first archetype
    ind_t <- sample(N, 1L)
    b <- ind_t # queue of active archetype indices
    sum_dist <- numeric(N)
    eligible <- rep(TRUE, N)
    eligible[ind_t] <- FALSE

    get_dists <- function(ind) {
        if (is.null(distances)) {
            return(.aa_dist2(X, X[ind, , drop = FALSE]))
        }
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

aa_pp <- function(X,
                  K,
                  sparse = inherits(X, "sparseMatrix"),
                  batch_size = NULL,
                  batch_type = c("distal", "uniform"),
                  batch_replace = FALSE,
                  ...) {
    # A++ initialization for Archetypal Analysis - Mair and Brefeld, 2019

    # if K is 2 AA++ reduces to kmeans++
    if (K == 2 && is.null(batch_size)) {
        return(kmeans_pp(X, K, ...))
    }

    b <- integer(K) # indices of archetypes

    # 1) randomly select the first archetype
    first_batch <- .aa_sample(
        X,
        size = batch_size,
        type = batch_type,
        replace = batch_replace
    )
    b[1L] <- sample(first_batch, 1L)

    # 2) sample the second archetype from points distal to the first
    available <- setdiff(seq_len(nrow(X)), b[1L])
    second_batch_size <- if (is.null(batch_size)) NULL else min(batch_size, length(available))
    second_batch <- .aa_sample(
        X[available, , drop = FALSE],
        size = second_batch_size,
        type = batch_type,
        replace = batch_replace
    )
    second_candidates <- available[second_batch]
    b[2L] <- second_candidates[
        .aa_sample(.aa_dist2(X[second_candidates, , drop = FALSE], X[b[1L], ]),
            size = 1L,
            replace = TRUE
        )
    ]
    if (K == 2) {
        return(b)
    }

    # 3) compute the first K-2 archetypes by iteratively running AA and sampling
    # from the points distal to the current archetype convex-hull.
    eps <- ifelse(sparse, 1e-8, 0)
    for (k in 3:K) {
        A <- X[b[1:(k - 1)], , drop = FALSE] # current archetypes

        available <- setdiff(seq_len(nrow(X)), b[1:(k - 1)])
        current_batch_size <- if (is.null(batch_size)) NULL else min(batch_size, length(available))
        batch <- .aa_sample(
            X[available, , drop = FALSE],
            size = current_batch_size,
            type = batch_type,
            replace = batch_replace
        )
        candidates <- available[batch]
        X_candidates <- X[candidates, , drop = FALSE]
        S <- proj_l1(.aa_solve_nnls(X_candidates, t(A)), eps = eps)
        res <- X_candidates - S %*% A
        dists <- rowSums(res * res) # squared residuals
        dists[!is.finite(dists)] <- 0
        b[k] <- candidates[
            .aa_sample(dists, size = 1L, replace = TRUE)
        ]
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
    candidate_rows <- switch(hull_method,
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
    if (is.null(facets)) {
        return(seq_len(N))
    }

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
    if (!(is_count(projected_dim) && projected_dim <= D)) {
        stop("`projected_dim` must be an integer between 1 and `ncol(X)`.", call. = FALSE)
    }

    if (projected_dim == D) {
        return(.aa_build_hull_candidates_full(X))
    }

    # Projected-hull strategy: aggregate envelope points across low-dimensional views.
    # Views are constructed by projecting (selecting) on random subsets of features.
    combs <- utils::combn(D, projected_dim)
    n_combs <- ncol(combs)
    selected <- seq_len(n_combs)
    if (!is.null(n_projection_max) && n_projection_max < n_combs) {
        selected <- sample.int(n_combs, size = n_projection_max, replace = FALSE)
    }

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
        if (length(grp) == 0L) {
            next
        }
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
    if (length(candidate_rows) == 0L) {
        stop("No hull candidates available for `hull_outmost`.")
    }

    if (use_unique_candidates) {
        candidate_rows <- unique(candidate_rows)
    }

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
.aa_dist_to_nearest_archetype <- function(X, ind, distances = NULL) {
    if (!is.null(distances)) {
        return(matrixStats::rowMins(distances[, ind, drop = FALSE]))
    }
    A <- X[ind, , drop = FALSE] # archetypes
    dists <- .aa_pdist2(A, X)
    matrixStats::colMins(dists)
}

# Sample points proportionally to their distance from a reference point
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
.aa_ind_to_init <- function(X, ind, sparse) {
    # make sure ind is positive indices selecting rows
    ind <- .aa_normalize_row_indices(ind, nrow(X), rownames(X))
    nm <- names(ind)

    A <- X[ind, , drop = FALSE] # Archetypes
    B <- onehot(ind, sparse = sparse, nrow(X)) # Row-stochastic matrix
    rownames(A) <- rownames(B) <- nm
    list(A = A, B = B)
}

.aa_init_S <- function(X, A, eps = 0) {
    S <- proj_l1(1 / .aa_pdist2(X, A), eps = eps) # init S by similarity score
    S[is.nan(S)] <- 1 # NaNs = Inf/Inf for the archetypes
    S
}

dirichlet <- function(X, K, alpha = 1, ...) {
    if (!is_positive(alpha)) {
        stop("`alpha` must be a single positive number.", call. = FALSE)
    }
    N <- nrow(X)
    B <- matrix(stats::rgamma(K * N, shape = alpha), nrow = K, ncol = N)
    B <- proj_l1(B, eps = 0)
    A <- B %*% X
    list(A = A, B = B)
}
