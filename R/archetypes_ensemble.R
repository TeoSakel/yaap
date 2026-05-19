#' Tune Archetypal Analysis Models
#'
#' Fits a grid of archetypal analysis models over numbers of archetypes,
#' replicates, and tuning arguments, evaluates each fit, and stores the result
#' as an archetypes ensemble.
#'
#' @param data Data passed to [run_aa()].
#' @param fit_method Fitting method passed as `method` to [run_aa()].
#' @param K Integer vector with the numbers of archetypes to fit.
#' @param eval Evaluation function. Defaults to [stats::AIC()].
#' @param eval_name Optional name for the evaluation metric column.
#' @param nrep Number of replicates per grid point.
#' @param direction Preferred optimization direction, `"minimize"` or
#'   `"maximize"`. Defaults to `"minimize"` for AIC.
#' @param ... Additional arguments passed to [run_aa()]. Atomic vectors of
#'   length greater than one are expanded as tuning dimensions. Plain list
#'   arguments are kept fixed; wrap a list in [I()] to expand it as a tuning
#'   dimension.
#'
#' @return An object of class `archetypes_ensemble`.
#'
#' @export
tune_archetypes <- function(data,
                            fit_method = c("pgd", "nnls", "kernel", "directional", "paa"),
                            K,
                            eval = stats::AIC,
                            eval_name = NULL,
                            nrep = 1L,
                            direction = NULL,
                            ...) {
    fit_method <- match.arg(fit_method)
    if (missing(K))
        stop("`K` must be supplied.", call. = FALSE)
    if (!is.numeric(K) || length(K) < 1L || any(!is.finite(K)) ||
        any(K != as.integer(K)) || any(K < 1L)) {
        stop("`K` must be a positive integer vector.", call. = FALSE)
    }
    K <- as.integer(K)
    stopifnot("`nrep` must be a single positive integer." = is_count(nrep))
    if (!is.function(eval))
        stop("`eval` must be a function.", call. = FALSE)

    eval_name <- eval_name %||% .aa_eval_name(substitute(eval))
    stopifnot("`eval_name` must be a single non-empty string" = is_non_empty_string(eval_name))
    reserved_metrics <- c("model_id", "K", "replicate", "converged", "loss", "r2")
    if (eval_name %in% reserved_metrics)
        stop("`eval_name` must not overwrite an internal metrics column.", call. = FALSE)

    direction <- match.arg(direction %||% "minimize", c("minimize", "maximize"))

    dots <- list(...)
    if (length(dots) > 0L && (is.null(names(dots)) || any(!nzchar(names(dots)))))
        stop("All tuning arguments in `...` must be named.", call. = FALSE)
    if ("method" %in% names(dots))
        stop("Use `fit_method`, not `method`, when tuning archetypes.", call. = FALSE)
    grid_spec <- .aa_tune_grid_spec(K = K, nrep = nrep, dots = dots)
    grid <- grid_spec[["grid"]]
    varying_names <- grid_spec[["varying_names"]]
    fixed_args <- grid_spec[["fixed_args"]]

    fits <- vector("list", nrow(grid))
    names(fits) <- grid[["model_id"]]
    metrics <- grid
    metrics[["converged"]] <- NA
    metrics[["loss"]] <- NA_real_
    metrics[["r2"]] <- NA_real_
    metrics[[eval_name]] <- NA_real_

    for (i in seq_len(nrow(grid))) {
        fit_args <- fixed_args
        for (nm in varying_names)
            fit_args[[nm]] <- grid[[nm]][[i]]

        fit <- do.call(
            run_aa,
            c(
                list(x = data, K = grid[["K"]][[i]], method = fit_method),
                fit_args
            )
        )
        fits[[i]] <- fit
        fit_metrics <- .aa_fit_metrics(fit)
        metrics[["converged"]][i] <- fit_metrics[["converged"]]
        metrics[["loss"]][i] <- fit_metrics[["loss"]]
        metrics[["r2"]][i] <- fit_metrics[["r2"]]
        metrics[[eval_name]][i] <- .aa_eval_fit(eval, fit)
    }

    .aa_new_archetypes_ensemble(
        data = data,
        fits = fits,
        grid = grid,
        metrics = metrics,
        prefer_metric = eval_name,
        prefer_direction = direction,
        eval_fun = eval,
        eval_name = eval_name,
        call = match.call()
    )
}

.aa_new_archetypes_ensemble <- function(data,
                                    fits,
                                    grid,
                                    metrics,
                                    prefer_metric,
                                    prefer_direction,
                                    eval_fun,
                                    eval_name,
                                    call = NULL) {
    stopifnot("`fits` must be a list" = is.list(fits))
    stopifnot("`grid` must be a data frame" = inherits(grid, "data.frame"))
    stopifnot("`metrics` must be a data frame" = inherits(metrics, "data.frame"))
    stopifnot("`fits` and `metrics` must have the same length" =
                  length(fits) == nrow(metrics))
    if (!prefer_metric %in% names(metrics))
        stop("`prefer_metric` must be a column in `metrics`.", call. = FALSE)
    prefer_direction <- match.arg(prefer_direction, c("minimize", "maximize"))

    structure(
        list(
            data = data,
            fits = fits,
            grid = grid,
            metrics = metrics,
            prefer_metric = prefer_metric,
            prefer_direction = prefer_direction,
            eval_fun = eval_fun,
            eval_name = eval_name,
            call = call
        ),
        class = "archetypes_ensemble"
    )
}

#' Select the Best Model
#'
#' @param x An object.
#' @param ... Arguments passed to methods.
#'
#' @export
best <- function(x, ...)
    UseMethod("best")

#' @rdname best
#' @param metric Metric column to optimize. Defaults to the ensemble preferred
#'   metric.
#' @param direction Optimization direction. Defaults to the ensemble preferred
#'   direction.
#'
#' @return A list with `metrics` and `fit`.
#'
#' @exportS3Method
best.archetypes_ensemble <- function(x,
                                     metric = NULL,
                                     direction = NULL,
                                     ...) {
    metric <- metric %||% x[["prefer_metric"]]
    stopifnot("`metric` must be a single non-empty string" = is_non_empty_string(metric))
    stopifnot("`metric` must be a column in `x$metrics`" = metric %in% names(x[["metrics"]]))

    direction <- direction %||% x[["prefer_direction"]]
    direction <- match.arg(direction, c("minimize", "maximize"))

    values <- x[["metrics"]][[metric]]
    if (!is.numeric(values))
        stop("`metric` must refer to a numeric column.", call. = FALSE)
    if (all(is.na(values)))
        stop("All values for `metric` are `NA`.", call. = FALSE)

    which_best <- if (identical(direction, "minimize")) which.min(values) else which.max(values)

    list(
        metrics = x[["metrics"]][which_best, , drop = FALSE],
        fit = x[["fits"]][[which_best]]
    )
}

#' Consistency Between Archetypal Analysis Fits
#'
#' @param x An object.
#' @param y Optional object to compare with `x`.
#' @param ... Arguments passed to methods.
#'
#' @details
#' For `what = "compositions"` or `what = "coefficients"`, consistency is
#' computed as normalized mutual information (NMI) between the two row-stochastic
#' weight matrices. For matrices `X` and `Y`, the joint distribution is
#' `Pxy = crossprod(X, Y) / nrow(X)`. Mutual information is computed from `Pxy`
#' and its row and column marginals, and normalized as
#' `2 * MI(X, Y) / (MI(X, X) + MI(Y, Y))`.
#'
#' For `what = "coordinates"`, consistency is defined when `K_x <= K_y`.
#' Archetypes in `x` are greedily matched to their nearest unmatched archetype
#' in `y` by squared Euclidean distance, and the score is
#' `1 - mean(d2) / mean(matrixStats::colVars(data))`. When `K_x > K_y`, the
#' score is `NA`.
#'
#' @references
#' J. L. Hinrich, S. E. Bardenfleth, R. E. Røge, N. W. Churchill, K. H.
#' Madsen, and M. Mørup, "Archetypal analysis for modeling multisubject
#' fMRI data," *IEEE Journal of Selected Topics in Signal Processing*,
#' vol. 10, no. 7, pp. 1160-1171, 2016.
#'
#' @export
consistency <- function(x, y, ...)
    UseMethod("consistency")

#' @rdname consistency
#' @param what Component to compare for archetypes objects: compositions,
#'   coefficients, or coordinates.
#' @param data Optional input data used as the variance reference for coordinate
#'   consistency. Defaults to `x$data`.
#'
#' @return A numeric scalar for two archetypes objects, and a numeric pairwise
#'   consistency matrix for an archetypes ensemble.
#'
#' @exportS3Method
consistency.archetypes <- function(x,
                                   y,
                                   what = c("compositions", "coefficients", "coordinates"),
                                   data = NULL,
                                   ...) {
    if (missing(y))
        stop("`y` must be supplied when comparing archetypes objects.", call. = FALSE)
    what <- match.arg(what)
    switch(
        what,
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

#' @rdname consistency
#' @exportS3Method
consistency.kernel_archetypes <- consistency.archetypes

#' @rdname consistency
#' @exportS3Method
consistency.archetypes_ensemble <- function(x,
                                            y,
                                            what = c("compositions", "coefficients", "coordinates"),
                                            ...) {
    if (!missing(y) && is_single_string(y))
        what <- y
    else if (!missing(y))
        stop("`y` is not used when computing pairwise ensemble consistency.",
             call. = FALSE)
    what <- match.arg(what)
    fits <- x[["fits"]]
    n <- length(fits)
    out <- matrix(NA_real_, nrow = n, ncol = n)
    dimnames(out) <- list(names(fits), names(fits))

    if (what %in% c("compositions", "coefficients")) {
        for (i in seq_len(n)) {
            for (j in i:n) {
                score <- if (i == j) 1 else consistency(fits[[i]], fits[[j]], what = what, ...)
                out[i, j] <- out[j, i] <- score
            }
        }
        return(out)
    }

    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            out[i, j] <- consistency(
                fits[[i]],
                fits[[j]],
                what = "coordinates",
                data = x[["data"]],
                ...
            )
        }
    }
    out
}

#' @exportS3Method
summary.archetypes_ensemble <- function(object, ...)
    object[["metrics"]]

.aa_eval_name <- function(expr) {
    label <- paste(deparse(expr), collapse = "")
    if (identical(label, "stats::AIC") || identical(label, "AIC"))
        return("AIC")
    label <- sub("^.*::", "", label)
    if (!nzchar(label) || identical(label, "eval"))
        label <- "metric"
    label
}

.aa_tune_grid_spec <- function(K, nrep, dots) {
    dims <- list(K = K, replicate = seq_len(nrep))
    fixed_args <- list()
    varying_names <- character()

    for (nm in names(dots)) {
        value <- dots[[nm]]
        if (inherits(value, "AsIs") && is.list(value)) {
            dims[[nm]] <- unclass(value)
            varying_names <- c(varying_names, nm)
        } else if (is.atomic(value) && length(value) > 1L) {
            dims[[nm]] <- value
            varying_names <- c(varying_names, nm)
        } else {
            fixed_args[[nm]] <- value
        }
    }

    if (any(vapply(dims, length, integer(1L)) < 1L))
        stop("Tuning dimensions must not be empty.", call. = FALSE)

    idx <- expand.grid(
        lapply(dims, seq_along),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    names(idx) <- names(dims)

    grid <- data.frame(
        model_id = sprintf("model_%03d", seq_len(nrow(idx))),
        stringsAsFactors = FALSE
    )
    for (nm in names(dims)) {
        values <- lapply(idx[[nm]], function(i) dims[[nm]][[i]])
        if (all(vapply(values, .aa_is_scalar_atomic, logical(1L)))) {
            grid[[nm]] <- unlist(values, recursive = FALSE, use.names = FALSE)
        } else {
            grid[[nm]] <- I(values)
        }
    }
    grid[["K"]] <- as.integer(grid[["K"]])
    grid[["replicate"]] <- as.integer(grid[["replicate"]])

    list(grid = grid, fixed_args = fixed_args, varying_names = varying_names)
}

.aa_is_scalar_atomic <- function(x) {
    is.atomic(x) && length(x) == 1L
}

.aa_fit_metrics <- function(fit) {
    loss <- fit[["loss"]]
    final_loss <- NA_real_
    final_r2 <- NA_real_
    if (inherits(loss, "data.frame") && nrow(loss) > 0L) {
        if ("loss" %in% names(loss))
            final_loss <- utils::tail(loss[["loss"]], 1L)
        if ("r2" %in% names(loss))
            final_r2 <- utils::tail(loss[["r2"]], 1L)
    }
    c(
        converged = isTRUE(fit[["converged"]]),
        loss = final_loss,
        r2 = final_r2
    )
}

.aa_eval_fit <- function(eval, fit) {
    value <- eval(fit)
    if (!is.numeric(value) || length(value) != 1L || is.na(value))
        stop("Evaluation function must return a single non-NA numeric value.", call. = FALSE)
    unname(value)
}

.aa_consistency_weights <- function(fit, what) {
    x <- switch(
        what,
        compositions = compositions(fit),
        coefficients = t(fit[["coefficients"]])
    )
    x <- as.matrix(x)
    row_total <- rowSums(x)
    ok <- is.finite(row_total) & row_total > 0
    if (any(!ok))
        stop("Consistency weights contain rows with zero or non-finite total.", call. = FALSE)
    x / row_total
}

.aa_nmi <- function(x, y) {
    stopifnot("Consistency matrices must have the same number of rows." = nrow(x) == nrow(y))
    stopifnot("Consistency matrices must be row-stochastic." = is_row_stochastic(x) && is_row_stochastic(y))
    n <- nrow(x)
    pxy <- crossprod(x, y) / n  # joint distribution of archetype co-membership
    pxx <- crossprod(x) / n  # joint distribution of archetype co-membership in x
    pyy <- crossprod(y) / n  # joint distribution of archetype co-membership in y
    2 * .aa_mi_from_joint(pxy) / (.aa_mi_from_joint(pxx) + .aa_mi_from_joint(pyy))
}

.aa_mi_from_joint <- function(pxy) {
    total <- sum(pxy)
    if (!isTRUE(all.equal(total, 1, tolerance = 1e-8)))
        pxy <- pxy / total
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
    if (kx > ky)
        return(NA_real_)
    data <- .aa_consistency_data_matrix(data)
    if (ncol(ax) != ncol(data) || ncol(ay) != ncol(data))
        stop("Coordinate consistency requires input-space coordinates.", call. = FALSE)

    d2 <- .aa_greedy_coordinate_d2(ax, ay)
    denom <- mean(matrixStats::colVars(data))
    if (!is.finite(denom) || denom <= 0)
        stop("Coordinate consistency requires data with positive column variance.", call. = FALSE)
    1 - mean(d2) / denom
}

.aa_consistency_coordinates <- function(fit) {
    if (inherits(fit, "kernel_archetypes")) {
        coords <- coordinates(fit)
        if (is.null(coords))
            stop("Kernel archetypes require `coordinates` for coordinate consistency.",
                 call. = FALSE)
        return(as.matrix(coords))
    }
    coords <- coordinates(fit)
    if (is.null(coords))
        stop("Coordinate consistency requires fitted coordinates.", call. = FALSE)
    if (inherits(coords, "fd"))
        return(.aa_fd_to_matrix(coords))
    as.matrix(coords)
}

.aa_consistency_data_matrix <- function(data) {
    if (inherits(data, "fd"))
        data <- .aa_fd_to_matrix(data)
    data <- as.matrix(data)
    if (!is.numeric(data))
        stop("Coordinate consistency requires numeric input data.", call. = FALSE)
    data
}

.aa_greedy_coordinate_d2 <- function(ax, ay) {
    available <- seq_len(nrow(ay))
    d2 <- numeric(nrow(ax))
    distances <- .aa_pdist2(ax, ay)
    for (i in seq_len(nrow(ax))) {
        best <- which.min(distances[i, available])
        d2[[i]] <- distances[i, available[[best]]]
        available <- available[-best]
    }
    d2
}
