#' Consistency Between Archetypal Analysis Fits
#'
#' Measures how similar two archetypal analysis fits are by comparing their
#' memberships, archetype definitions, or archetype locations.
#'
#' @param x An object.
#' @param y A fitted archetypes object to compare with `x`.
#' @param ... Arguments passed to methods.
#'
#' @details
#' Consistency is useful when checking whether an AA solution is stable across
#' random starts, resampled data, or nearby model choices. A high score means the
#' two fits describe the data in a similar way, even if archetype labels have
#' changed. Compare `compositions` to ask whether samples are assigned to the
#' same archetypal mixtures, `coefficients` to ask whether archetypes are built
#' from the same samples, and `coordinates` to ask whether the archetype
#' locations are close in the input space.
#'
#' For `what = "compositions"` or `what = "coefficients"`, the two components
#' are treated as row-stochastic membership matrices and compared with
#' normalized mutual information (NMI). For matrices `X` and `Y`, the joint
#' distribution is `Pxy = crossprod(X, Y) / nrow(X)`. Mutual information is
#' computed from `Pxy` and its row and column marginals, then normalized as
#' `2 * MI(X, Y) / (MI(X, X) + MI(Y, Y))`. This makes the score insensitive to
#' permutations of the archetype labels.
#'
#' For `what = "coordinates"`, archetypes are greedily matched by repeatedly selecting the
#' closest remaining pair across the two fits, then a similarity-style score is
#' computed by scaling the mean matched squared distances by the average
#' squared pairwise distance between observations in `data`. This denominator
#' is computed through the equivalent identity `2 * sum(colVars(data))`,
#' avoiding explicit pairwise distances. The score is useful for checking
#' whether two fits place archetypes in similar regions of the input space.
#' Values can be negative when matched archetypes are farther apart than the
#' typical data-to-data squared distance. It returns `NA` when `x` has more
#' archetypes than `y`, because not every archetype in `x` can be matched.
#'
#' @references
#' J. L. Hinrich, S. E. Bardenfleth, R. E. Røge, N. W. Churchill, K. H.
#' Madsen, and M. Mørup, "Archetypal analysis for modeling multisubject
#' fMRI data," *IEEE Journal of Selected Topics in Signal Processing*,
#' vol. 10, no. 7, pp. 1160-1171, 2016.
#'
#' @export
consistency <- function(x, y, ...) {
    UseMethod("consistency")
}

#' @rdname consistency
#' @param what Component to compare: compositions, coefficients, or coordinates.
#' @param data Optional input data used as the variance reference for coordinate
#'   consistency. Defaults to `x$data`.
#'
#' @return A numeric scalar.
#'
#' @exportS3Method
consistency.archetypes <- function(x,
                                   y,
                                   what = c("compositions", "coefficients", "coordinates"),
                                   data = NULL,
                                   ...) {
    if (missing(y)) {
        stop("`y` must be supplied when comparing archetypes objects.", call. = FALSE)
    }
    what <- match.arg(what)
    switch(what,
        compositions = .aa_nmi(
            .aa_consistency_weights(x, "compositions"),
            .aa_consistency_weights(y, "compositions")
        ),
        coefficients = .aa_nmi(
            .aa_consistency_weights(x, "coefficients"),
            .aa_consistency_weights(y, "coefficients")
        ),
        coordinates = {
            data <- data %||% x[["data"]]
            .aa_coordinate_consistency(x, y, data)
        }
    )
}

.aa_consistency_weights <- function(fit, what) {
    x <- switch(what,
        compositions = compositions(fit),
        coefficients = t(fit[["coefficients"]])
    )
    x <- as.matrix(x)
    row_total <- rowSums(x)
    ok <- is.finite(row_total) & row_total > 0
    if (any(!ok)) {
        stop("Consistency weights contain rows with zero or non-finite total.", call. = FALSE)
    }
    x / row_total
}

.aa_nmi <- function(x, y) {
    if (nrow(x) != nrow(y)) {
        stop("Consistency matrices must have the same number of rows.", call. = FALSE)
    }
    if (!(is_row_stochastic(x) && is_row_stochastic(y))) {
        stop("Consistency matrices must be row-stochastic.", call. = FALSE)
    }
    n <- nrow(x)
    pxy <- crossprod(x, y) / n # joint distribution of archetype co-membership
    pxx <- crossprod(x) / n # joint distribution of archetype co-membership in x
    pyy <- crossprod(y) / n # joint distribution of archetype co-membership in y
    2 * .aa_mi_from_joint(pxy) / (.aa_mi_from_joint(pxx) + .aa_mi_from_joint(pyy))
}

.aa_mi_from_joint <- function(pxy) {
    total <- sum(pxy)
    if (!isTRUE(all.equal(total, 1, tolerance = 1e-8))) {
        pxy <- pxy / total
    }
    px <- rowSums(pxy)
    py <- colSums(pxy)
    expected <- outer(px, py)
    keep <- pxy > 0 & expected > 0
    sum(pxy[keep] * log(pxy[keep] / expected[keep]))
}

.aa_coordinate_consistency <- function(x, y, data) {
    ax <- .aa_consistency_coordinates(x)
    ay <- .aa_consistency_coordinates(y)
    kx <- nrow(ax)
    ky <- nrow(ay)
    if (kx > ky) {
        return(NA_real_)
    }
    data <- .aa_consistency_data_matrix(data)
    if (ncol(ax) != ncol(data) || ncol(ay) != ncol(data)) {
        stop("Coordinate consistency requires input-space coordinates.", call. = FALSE)
    }

    d2 <- .aa_greedy_coordinate_d2(ax, ay)
    denom <- 2 * sum(matrixStats::colVars(data))
    if (!is.finite(denom) || denom <= 0) {
        stop("Coordinate consistency requires data with positive column variance.", call. = FALSE)
    }
    1 - mean(d2) / denom
}

.aa_consistency_coordinates <- function(fit) {
    if (inherits(fit, "kernel_archetypes")) {
        coords <- coordinates(fit)
        if (is.null(coords)) {
            stop("Kernel archetypes require `coordinates` for coordinate consistency.",
                call. = FALSE
            )
        }
        return(as.matrix(coords))
    }
    coords <- coordinates(fit)
    if (is.null(coords)) {
        stop("Coordinate consistency requires fitted coordinates.", call. = FALSE)
    }
    if (inherits(coords, "fd")) {
        return(.aa_fd_to_matrix(coords))
    }
    as.matrix(coords)
}

.aa_consistency_data_matrix <- function(data) {
    if (inherits(data, "fd")) {
        data <- .aa_fd_to_matrix(data)
    }
    data <- as.matrix(data)
    if (!is.numeric(data)) {
        stop("Coordinate consistency requires numeric input data.", call. = FALSE)
    }
    data
}

.aa_greedy_coordinate_d2 <- function(ax, ay) {
    if (nrow(ax) > nrow(ay)) {
        stop("`ay` must have at least as many rows as `ax`.", call. = FALSE)
    }

    available_x <- seq_len(nrow(ax))
    available_y <- seq_len(nrow(ay))
    d2 <- numeric(nrow(ax))
    distances <- .aa_pdist2(ax, ay)
    for (i in seq_len(nrow(ax))) {
        remaining <- distances[available_x, available_y, drop = FALSE]
        best <- which(remaining == min(remaining), arr.ind = TRUE)[1L, ]
        d2[i] <- remaining[best[1L], best[2L]]
        available_x <- available_x[-best[1L]]
        available_y <- available_y[-best[2L]]
    }
    d2
}
