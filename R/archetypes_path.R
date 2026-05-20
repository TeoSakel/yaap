#' Fit an Archetypes Path Across K
#'
#' Fits a sequence of archetypal analysis models over candidate numbers of
#' archetypes. The input data are checked and preprocessed once, then the
#' selected solver's fit step is run independently for each candidate `K`.
#'
#' @param x numeric matrix, data frame, formula, or supported specialized input.
#' @param formula formula selecting variables from `data`. The response, when
#'   present, is ignored.
#' @param data optional data frame supplying variables for formula input.
#' @param K candidate numbers of archetypes. A single value expands to `1:K`;
#'   a vector is sorted and deduplicated.
#' @param ... additional arguments passed to [run_aa()].
#'
#' @return An object of class `archetypes_path` containing one `archetypes` fit
#'   per candidate K. Use `[[` with a numeric index or model name such as
#'   `"K3"` to extract one fitted model with the shared data restored. Use `[`
#'   to subset the path while preserving the `archetypes_path` container. Use
#'   `$` for metadata such as `K`, `method`, `family`, and `data`.
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' path <- archetypes_path(as.matrix(toy), K = 3, max_iter = 20)
#' screeplot(path)
#' path[["K2"]]
#'
#' @export
archetypes_path <- function(x, ...) {
    UseMethod("archetypes_path")
}

#' @rdname archetypes_path
#' @param subset optional expression selecting rows before fitting formula input.
#' @param na.action function controlling missing-value handling for formula
#'   input. Defaults to [stats::na.omit()].
#' @method archetypes_path formula
#' @export
archetypes_path.formula <- function(formula,
                                    data = NULL,
                                    K,
                                    ...,
                                    subset,
                                    na.action) {

    path_args <- .aa_path_dots(...)
    call <- path_args[["call"]] %||% match.call()
    if (is.null(path_args[["call"]]))
        call[[1L]] <- quote(archetypes_path)

    terms <- stats::delete.response(stats::terms(formula, data = data))
    mf_call <- match.call(expand.dots = FALSE)
    keep <- match(c("formula", "data", "subset", "na.action"), names(mf_call), 0L)
    mf_call <- mf_call[c(1L, keep)]
    mf_call[[1L]] <- quote(stats::model.frame)
    mf_call[["formula"]] <- terms
    mf <- eval(mf_call, parent.frame())

    terms <- attr(mf, "terms")
    attr(terms, "intercept") <- 0L
    X <- stats::model.matrix(terms, mf)
    intercept <- colnames(X) == "(Intercept)"
    if (any(intercept))
        X <- X[, !intercept, drop = FALSE]
    if (ncol(X) == 0L)
        stop("Formula input must select at least one predictor column.", call. = FALSE)

    out <- do.call(
        .aa_archetypes_path_default,
        c(list(x = X, K = K, call = call), path_args[["dots"]]),
        quote = TRUE
    )
    out[["formula"]] <- formula
    out[["terms"]] <- terms
    out
}

#' @rdname archetypes_path
#' @method archetypes_path fd
#' @export
archetypes_path.fd <- function(x, K, ...) {
    data <- x
    .aa_require_namespace_for("fda", "functional-data archetypes paths")

    method_args <- list(...)
    if ("scale" %in% names(method_args)) {
        msg <- paste("`scale` is computed from the fd basis and",
                     "cannot be supplied to `archetypes_path.fd()`.")
        stop(msg, call. = FALSE)
    }

    coefs <- t(stats::coef(data))
    stopifnot("`fd` must have a 2-D coefficient matrix." = length(dim(coefs)) == 2L)

    G <- fda::eval.penalty(data[["basis"]], Lfdobj = 0)

    path_args <- .aa_path_dots(...)
    call <- path_args[["call"]] %||% match.call()
    if (is.null(path_args[["call"]]))
        call[[1L]] <- quote(archetypes_path)

    out <- do.call(
        .aa_archetypes_path_default,
        c(list(x = coefs, K = K, call = call, scale = G), path_args[["dots"]]),
        quote = TRUE
    )
    out[["data"]] <- data
    out
}

#' @rdname archetypes_path
#' @inheritParams run_aa
#' @method archetypes_path default
#' @export
archetypes_path.default <- function(x,
                                    K,
                                    method = c("pgd", "nnls", "kernel", "directional", "paa"),
                                    family = "gaussian",
                                    init = NULL,
                                    init_args = list(),
                                    weights = NULL,
                                    scale = FALSE,
                                    robust = FALSE,
                                    robust_args = list(),
                                    sd_threshold = 1e-6,
                                    max_iter = 100L,
                                    tol = 1e-6,
                                    tol_r2 = 0.9999,
                                    max_kappa = 1000,
                                    eps = NULL,
                                    verbose = FALSE,
                                    missing = NULL,
                                    nrep = 1L,
                                    ...) {
    path_args <- .aa_path_dots(...)
    call <- path_args[["call"]] %||% match.call()
    if (is.null(path_args[["call"]]))
        call[[1L]] <- quote(archetypes_path)
    do.call(
        .aa_archetypes_path_default,
        c(
            list(
                x = x,
                K = K,
                call = call,
                method = method,
                family = family,
                init = init,
                init_args = init_args,
                weights = weights,
                scale = scale,
                robust = robust,
                robust_args = robust_args,
                sd_threshold = sd_threshold,
                max_iter = max_iter,
                tol = tol,
                tol_r2 = tol_r2,
                max_kappa = max_kappa,
                eps = eps,
                verbose = verbose,
                missing = missing,
                nrep = nrep
            ),
            path_args[["dots"]]
        ),
        quote = TRUE
    )
}

.aa_archetypes_path_default <- function(x,
                                        K,
                                        call,
                                        method = c("pgd", "nnls", "kernel", "directional", "paa"),
                                        family = "gaussian",
                                        init = NULL,
                                        init_args = list(),
                                        weights = NULL,
                                        scale = FALSE,
                                        robust = FALSE,
                                        robust_args = list(),
                                        sd_threshold = 1e-6,
                                        max_iter = 100L,
                                        tol = 1e-6,
                                        tol_r2 = 0.9999,
                                        max_kappa = 1000,
                                        eps = NULL,
                                        verbose = FALSE,
                                        missing = NULL,
                                        nrep = 1L,
                                        ...) {
    K       <- .aa_path_normalize_K(K, expand_scalar = TRUE)
    data    <- if (inherits(x, "data.frame")) as.matrix(x) else x
    eps     <- eps %||% ifelse(inherits(data, "sparseMatrix"), 0, 1e-8)
    missing <- missing %||% any(is.na(data))

    .aa_fit_path_engine(
        call = call,
        x = data,
        K = K,
        method = method,
        family = family,
        init = init,
        init_args = init_args,
        weights = weights,
        scale = scale,
        robust = robust,
        robust_args = robust_args,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        missing = missing,
        nrep = nrep,
        ...
    )
}

.aa_new_archetypes_path <- function(models,
                                    data,
                                    K,
                                    call,
                                    method,
                                    family = "gaussian") {
    stopifnot("`models` must be a list" = is.list(models))
    stopifnot("`K` must match `models`" = length(K) == length(models))
    if (!all(vapply(models, inherits, logical(1L), what = "archetypes")))
        stop("`models` must contain only archetypes fits.", call. = FALSE)
    names(models) <- names(models) %||% paste0("K", K)
    structure(
        list(
            models = models,
            data = data,
            K = K,
            call = call,
            method = method,
            family = family
        ),
        class = "archetypes_path"
    )
}

.aa_fit_path_engine <- function(call,
                                x,
                                K,
                                method = c("pgd", "nnls", "kernel", "directional", "paa"),
                                family = "gaussian",
                                init = NULL,
                                init_args = list(),
                                weights = NULL,
                                scale = FALSE,
                                robust = FALSE,
                                robust_args = list(),
                                sd_threshold = 1e-6,
                                max_iter = 100L,
                                tol = 1e-6,
                                tol_r2 = 0.9999,
                                max_kappa = 1000,
                                eps = ifelse(inherits(x, "sparseMatrix"), 0, 1e-8),
                                verbose = FALSE,
                                missing = any(is.na(x)),
                                nrep = 1L,
                                data = NULL,
                                ...) {
    K <- .aa_path_normalize_K(K, expand_scalar = FALSE)
    if (is_tabular(init))
        stop("`init` cannot be a matrix when fitting an archetypes path.", call. = FALSE)

    setup_call <- call
    setup_call[["K"]] <- max(K)
    setup <- .aa_fit_engine_setup(
        call = setup_call,
        x = x,
        K = max(K),
        method = method,
        family = family,
        init = init,
        init_args = init_args,
        weights = weights,
        scale = scale,
        robust = robust,
        robust_args = robust_args,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        missing = missing,
        nrep = nrep,
        data = data,
        ...
    )

    models <- vector("list", length(K))
    names(models) <- paste0("K", K)
    for (i in seq_along(K)) {
        ctx <- setup[["ctx"]]
        ctx[["K"]] <- K[[i]]
        ctx[["call"]][["K"]] <- K[[i]]
        fit <- .aa_fit_prepared(ctx, setup[["block"]], setup[["prep"]])
        fit[["data"]] <- NULL
        models[[i]] <- fit
    }

    .aa_new_archetypes_path(
        models = models,
        data = x,
        K = K,
        call = call,
        method = setup[["ctx"]][["method"]],
        family = family
    )
}


# Methods for archetypes_path objects ------------------------------------------------

#' @method [[ archetypes_path
#' @export
`[[.archetypes_path` <- function(x, i, ...) {
    raw <- unclass(x)

    # Extract fit
    models <- raw[["models"]]
    if (is_count(i)) {
        stopifnot("Model index is out of bounds." = i <= length(models))
        idx <- as.integer(i)
    } else if (is_single_string(i)) {
        idx <- match(i, names(models))
        if (is.na(idx))
            stop(sprintf("Unknown archetypes path model name: %s", i), call. = FALSE)
    } else {
        stop("`i` must be a single model index or model name.", call. = FALSE)
    }
    fit <- models[[idx]]

    # Inflate the fit with shared data and metadata from the path
    fit[["data"]] <- raw[["data"]]
    fit[["call"]][["K"]] <- raw[["K"]][[idx]]
    if (!is.null(raw[["formula"]]))
        fit[["formula"]] <- raw[["formula"]]
    if (!is.null(raw[["terms"]]))
        fit[["terms"]] <- raw[["terms"]]
    fit
}

#' @method [ archetypes_path
#' @export
`[.archetypes_path` <- function(x, i, ...) {
    raw <- unclass(x)
    models <- raw[["models"]]
    if (missing(i))
        i <- seq_along(models)

    idx <- if (is_all_count(i)) {
        stopifnot("Model index is out of bounds." = all(i <= length(models)))
        as.integer(i)
    } else if (is.character(i) && length(i) > 0L && !anyNA(i)) {
        match(i, names(models))
    } else if (is.logical(i) && !anyNA(i)) {
        stopifnot("Logical mask length must be equal to the number of models." =
                      length(i) == length(models))
        which(i)
    } else {
        stop("`i` must be model indices, model names, or a logical model mask.", call. = FALSE)
    }
    if (anyNA(idx))
        stop("Unknown archetypes path model.", call. = FALSE)

    out <- .aa_new_archetypes_path(
        models = models[idx],
        data = raw[["data"]],
        K = raw[["K"]][idx],
        call = raw[["call"]],
        method = raw[["method"]],
        family = raw[["family"]]
    )
    extra <- setdiff(names(raw), c("models", "data", "K", "call", "method", "family"))
    out[extra] <- raw[extra]
    out
}

#' @exportS3Method
length.archetypes_path <- function(x) {
    length(unclass(x)[["models"]])
}

#' @exportS3Method
names.archetypes_path <- function(x) {
    names(unclass(x)[["models"]])
}

#' @exportS3Method
print.archetypes_path <- function(x, ...) {
    raw <- unclass(x)
    cat("Archetypes path:\n")
    cat("Candidate K:", paste(raw[["K"]], collapse = ", "), "\n")
    cat("Method:", raw[["method"]], "\n")
    cat("Number of models:", length(x), "\n")
    invisible(x)
}

#' Scree Plot for an Archetypes Path
#'
#' Draws or prepares a model-selection curve over candidate `K` values.
#'
#' @param x An `archetypes_path` object.
#' @param y Metric to plot. `NULL` defaults to `"AIC"` for Euclidean Gaussian
#'   `pgd` and `nnls` paths, and to `"r2"` otherwise. A character value may be
#'   `"AIC"` or a column from each fit's `loss` table. A function is called on
#'   each extracted fit and must return a single numeric value.
#' @param plot Logical. Should the plot be drawn?
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns a data frame with one row per candidate K.
#'
#' @exportS3Method
screeplot.archetypes_path <- function(x, y = NULL, plot = TRUE, ...) {
    stopifnot("`plot` must be TRUE or FALSE." = is.logical(plot))
    raw <- unclass(x)

    # Extract metrics from fits
    fits      <- lapply(seq_along(raw[["models"]]), function(i) x[[i]])
    metric    <- .aa_path_metric_name(raw, y)
    values    <- vapply(fits, .aa_path_metric_value, numeric(1L), y = y,      metric = metric)
    loss      <- vapply(fits, .aa_path_metric_value, numeric(1L), y = "loss", metric = "loss")
    r2        <- vapply(fits, .aa_path_metric_value, numeric(1L), y = "r2",   metric = "r2")
    n_iter    <- vapply(fits, function(fit) nrow(fit[["loss"]]) - 1L, integer(1L))
    converged <- vapply(fits, function(fit) isTRUE(fit[["converged"]]), logical(1L))

    out <- data.frame(
        K = raw[["K"]],
        metric = rep(metric, length(values)),
        value = values,
        loss = loss,
        r2 = r2,
        n_iter = n_iter,
        converged = converged
    )

    if (plot) {
        args <- list(
            x = out[["K"]],
            y = out[["value"]],
            type = "b",
            pch = 19,
            xlab = "Number of archetypes (K)",
            ylab = metric
        ) %|p|% list(...)
        do.call(graphics::plot, args)
    }
    invisible(out)
}


# Internal helpers for archetypes_path methods ----------------------------------

.aa_path_dots <- function(...) {
    dots <- list(...)
    call <- dots[[".aa_call"]]
    dots[[".aa_call"]] <- NULL
    list(call = call, dots = dots)
}

.aa_path_normalize_K <- function(K, expand_scalar = FALSE) {
    stopifnot("`expand_scalar` must be TRUE or FALSE" = is_logical(expand_scalar))
    if (length(K) < 1L || !is_all_count(K))
        stop("`K` must be a positive integer vector.", call. = FALSE)
    K <- as.integer(K)
    if (expand_scalar && length(K) == 1L)
        K <- seq_len(K)
    sort(unique(K))
}

.aa_path_metric_name <- function(x, y) {
    if (is.null(y)) {
        method <- x[["method"]]
        family <- x[["family"]] %||% "gaussian"
        if (method %in% c("pgd", "nnls") && identical(family, "gaussian"))
            return("AIC")
        return("r2")
    }
    if (is.function(y))
        return("function")
    if (!is_non_empty_string(y))
        stop("`y` must be NULL, a single loss-column name, 'AIC', or a function.", call. = FALSE)
    y
}

.aa_path_metric_value <- function(fit, y, metric) {
    value <- if (is.function(y)) {
        y(fit)
    } else if (identical(metric, "AIC")) {
        tryCatch(AIC(fit), error = function(e) NA_real_, warning = function(w) NA_real_)
    } else {
        loss <- fit[["loss"]]
        if (is.null(loss[[metric]])) {
            fmt <- "`y` must name a column in the loss table or be 'AIC'; unknown metric '%s'."
            stop(sprintf(fmt, metric), call. = FALSE)
        }
        utils::tail(loss[[metric]], 1L)
    }
    if (!(is_number(value) || (identical(length(value), 1L) && is.numeric(value) && is.na(value))))
        stop("`y` function must return a single numeric value.", call. = FALSE)
    as.numeric(value)
}
