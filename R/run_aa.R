#' Run Archetypal Analysis
#'
#' Common entry point for fitting archetypal analysis models. Dispatches to the
#' solver selected by `method`. For solver-specific arguments, theoretical
#' background, and class-specific return slots, see the individual solver pages.
#'
#' @param x numeric matrix (rows = samples, columns = features), or an object
#'   with a class-specific `run_aa()` method.
#' @param formula formula selecting variables from `data`. The response, when
#'   present, is ignored.
#' @param data optional data frame supplying variables for formula input.
#' @param K number of archetypes. When `K` has length greater than one,
#'   `run_aa()` returns an `archetypes_path()` object with one fit per value.
#' @param ... additional arguments passed to the selected solver.
#'
#' @returns An object of class `archetypes` with components:
#' \describe{
#'   \item{`coordinates`}{(K x M) archetype coordinates in the original feature space.}
#'   \item{`coefficients`}{(K x N) weights expressing each archetype as a convex combination of samples.}
#'   \item{`compositions`}{(N x K) row-stochastic weights expressing each sample as a convex combination of archetypes.}
#'   \item{`loss`}{data frame of per-iteration metrics. All fitters include `loss` and `r2`; additional diagnostic columns are method-specific.}
#'   \item{`converged`}{logical convergence flag.}
#'   \item{`data`}{original data passed to the fitter.}
#'   \item{`call`}{the matched call.}
#'   \item{`family`}{observation family string (e.g. `"gaussian"`).}
#'   \item{`init`}{initial archetype coordinates (when available).}
#'   \item{`slack`, `weights`}{optional relaxation and sample-weight parameters.}
#'   \item{`feature_map`}{internal metadata mapping new data into the fitted optimization geometry for prediction.}
#' }
#'
#' @seealso
#' Solvers: [archetypes_pgd()], [archetypes_nnls()], [archetypes_kernel_pgd()],
#'   [archetypes_directional()], [archetypes_paa()].
#' Post-fit: [plot.archetypes()], [predict.archetypes()],
#'   [fitted.archetypes()], [residuals.archetypes()], [anames()].
#' Model selection: [AIC.archetypes()].
#' Tidy output: [tidy.archetypes()], [glance.archetypes()], [augment.archetypes()].
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' run_aa(as.matrix(toy), K = 3)
#' run_aa(as.matrix(toy), K = 3, method = "nnls")
#' run_aa(Species ~ ., data = iris, K = 3)
#'
#' \dontrun{
#' # Functional data example based on fda::growth
#' library(fda)
#' data(growth, package = "fda")
#'
#' hgtm <- t(growth$hgtm)
#' basis_fd <- create.bspline.basis(c(1, ncol(hgtm)), 10)
#' temp_fd <- Data2fd(
#'     argvals = seq_len(ncol(hgtm)),
#'     y = growth$hgtm,
#'     basisobj = basis_fd
#' )
#'
#' fit_fd <- run_aa(temp_fd, K = 3, max_iter = 20)
#' arch_fd <- coordinates(fit_fd)
#' }
#'
#' @export
run_aa <- function(x, ...) {
    UseMethod("run_aa")
}

#' @rdname run_aa
#' @param subset optional expression selecting rows before fitting formula input.
#' @param na.action function controlling missing-value handling for formula
#'   input. Defaults to [stats::na.omit()].
#'
#' @exportS3Method
run_aa.formula <- function(formula,
                           data = NULL,
                           K,
                           ...,
                           subset,
                           na.action) {

    call <- match.call()
    call[[1L]] <- quote(run_aa)

    if (length(K) > 1L) {
        # Repoint the call so archetypes_path S3 dispatch sees the original
        # input class; .aa_call preserves run_aa as the recorded public call.
        path_call <- call
        path_call[[1L]] <- quote(archetypes_path)
        path_call[[".aa_call"]] <- base::call("quote", call)
        return(eval(path_call, parent.frame()))
    }

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

    fit <- .aa_fit_engine(
        call = call,
        x = X,
        K = K,
        ...
    )
    fit[["formula"]] <- formula
    fit[["terms"]] <- terms
    fit
}

#' @rdname run_aa
#' @param method fitting method. One of `"pgd"`, `"nnls"`, `"kernel"`,
#'   `"directional"`, or `"paa"` (default: `"pgd"`).
#' @param family observation family passed to `method = "paa"`. Defaults to
#'   `"gaussian"`.
#' @param init initialization method for archetype starting coordinates.
#'   Accepts a function, a method name string, or a numeric (K x M) coordinate
#'   matrix. `NULL` selects `"furthest_sum"` for all methods except
#'   `"directional"`, which defaults to `"dirichlet"`. When a matrix is
#'   supplied it must have dimension `K x ncol(x)`; row names are used as
#'   archetype names. Available method strings:
#' \describe{
#'   \item{`"furthest_sum"`}{greedily maximises the sum of distances from the
#'     current archetype set (Mørup & Hansen 2012). Default for most methods.}
#'   \item{`"furthest_first"`}{greedy farthest-point selection.}
#'   \item{`"kmeans_pp"`}{probabilistic farthest-point selection (soft
#'     furthest-first).}
#'   \item{`"random"`}{uniformly random sample of K rows.}
#'   \item{`"dirichlet"`}{random convex combinations sampled from a
#'     Dirichlet distribution. Default for `"directional"`.}
#'   \item{`"aa_pp"`}{AA++ initialization (Mair & Sjölund 2023). Pass
#'     `batch_size` in `init_args` to use a Monte Carlo-inspired variant.}
#'   \item{`"hull_outmost"`}{hull-candidate outmost-vote ranking.}
#' }
#'   Any initializer can receive `batch_size`, `batch_type`, and
#'   `batch_replace` through `init_args`; `batch_type = "distal"` implements
#'   coreset-style candidate sampling (Mair & Brefeld 2019).
#'   See `vignette("initialization", package = "yaap")` for a comparison.
#' @param init_args list of additional arguments for the initialization function.
#' @param weights optional numeric vector of sample weights (default: `NULL`).
#'   Internally scaled to mean 1 and square-rooted before use.
#' @param scale common `run_aa()` scaling argument, present for consistency
#'   across method dispatch. Only Euclidean Gaussian methods (`"pgd"` and
#'   `"nnls"`) use it: `FALSE` (default) leaves columns on their original
#'   scale, `TRUE` applies z-score standardization, a positive numeric vector
#'   divides by user-supplied scale factors, and a symmetric positive-definite
#'   matrix applies the corresponding feature metric. Specialized methods
#'   (`"kernel"`, `"directional"`, and `"paa"`) define their own geometry or
#'   likelihood; non-`FALSE` values are ignored with a warning.
#' @param robust robust row reweighting selector. Use `FALSE` for ordinary
#'   squared error, `TRUE` for `"psi.bisquare"`, a MASS psi function name,
#'   or a custom psi function. See [MASS::rlm()] for psi details; `method =
#'   "MM"` is not supported because it is not applicable to AA.
#' @param robust_args list of tuning arguments passed to the robust psi function.
#' @param sd_threshold threshold for feature standard deviation below which
#'   columns are dropped before fitting (default: 1e-6).
#' @param max_iter maximum number of outer iterations (default: 100).
#' @param tol convergence tolerance on the residual sum of squares (default: 1e-6).
#' @param tol_r2 convergence tolerance on R\eqn{^2} (default: 0.9999).
#' @param max_kappa maximum condition number warning threshold for `method = "nnls"` (default: 1000).
#' @param eps small positive number for numerical stability
#'   (default: 0 for sparse input, 1e-8 for dense).
#' @param verbose whether to print progress messages (default: `FALSE`).
#' @param nrep number of random restarts; the best fit (lowest final loss)
#'   is returned (default: 1).
#' @param ... additional arguments passed to the selected solver. See
#'   [archetypes_pgd()], [archetypes_nnls()], [archetypes_kernel_pgd()],
#'   [archetypes_directional()], and [archetypes_paa()] for method-specific
#'   parameters.
#'
#' @exportS3Method
run_aa.default <- function(x,
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
    data <- if (inherits(x, "data.frame")) as.matrix(x) else x

    call <- match.call()
    call[[1L]] <- quote(run_aa)

    if (length(K) > 1L) {
        # Repoint the call so archetypes_path S3 dispatch sees the original
        # input class; .aa_call preserves run_aa as the recorded public call.
        path_call <- call
        path_call[[1L]] <- quote(archetypes_path)
        path_call[[".aa_call"]] <- base::call("quote", call)
        return(eval(path_call, parent.frame()))
    }

    eps <- eps %||% ifelse(inherits(data, "sparseMatrix"), 0, 1e-8)
    missing <- missing %||% any(is.na(data))

    .aa_fit_engine(
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

#' @rdname run_aa
#' @details
#' `run_aa.fd()` fits archetypal analysis to the coefficient matrix of an
#' `fda::fd` object. The feature scaling is set to the basis inner-product
#' matrix from `fda::eval.penalty(data$basis, Lfdobj = 0)`, so
#' `scale` cannot be supplied to this method. The internal optimization-space
#' archetype coordinates are basis coefficients; `coordinates()` returns them in
#' the original fd representation.
#'
#' @exportS3Method
run_aa.fd <- function(x, K, ...) {
    data <- x
    .aa_require_namespace_for("fda", "functional-data archetypal analysis")

    method_args <- list(...)
    if ("scale" %in% names(method_args))
        stop("`scale` is computed from the fd basis and cannot be supplied to `run_aa.fd()`.",
             call. = FALSE)

    call <- match.call()
    call[[1L]] <- quote(run_aa)

    if (length(K) > 1L) {
        # Repoint the call so archetypes_path S3 dispatch sees the original
        # input class; .aa_call preserves run_aa as the recorded public call.
        path_call <- call
        path_call[[1L]] <- quote(archetypes_path)
        path_call[[".aa_call"]] <- base::call("quote", call)
        return(eval(path_call, parent.frame()))
    }

    coefs <- stats::coef(data)
    if (length(dim(coefs)) != 2L)
        stop("`run_aa.fd()` currently supports fd objects with a 2-D coefficient matrix.",
             call. = FALSE)

    B <- t(coefs)
    G <- fda::eval.penalty(data[["basis"]], Lfdobj = 0)

    fit <- .aa_fit_engine(
        call = call,
        x = B,
        K = K,
        scale = G,
        ...
    )

    fit[["data"]] <- data
    fit
}

.aa_fit_engine_setup <- function(call,
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
    method <- match.arg(method, c("pgd", "nnls", "kernel", "directional", "paa"))
    init <- init %||% ifelse(identical(method, "directional"), "dirichlet", "furthest_sum")
    stopifnot("`missing` must be TRUE or FALSE" = is_logical(missing))
    stopifnot("`nrep` must be a single positive integer" = is_count(nrep))
    nrep <- as.integer(nrep)
    method_args <- list(...)
    precomputed_kernel <- identical(method, "kernel") &&
        identical(method_args[["kernel"]], "precomputed")
    .aa_check_x(x, missing = missing, validate_values = !precomputed_kernel)
    ctx <- list(
        call = call,
        x = x,
        data = data,
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
        nrep = nrep
    )

    block <- switch(
        method,
        pgd         = .aa_pgd_block(ctx, ...),
        nnls        = .aa_nnls_block(ctx, ...),
        kernel      = .aa_kernel_block(ctx, ...),
        directional = .aa_directional_block(ctx, ...),
        paa         = .aa_paa_block(ctx, ...)
    )
    block[["check"]](ctx)
    prep <- block[["preprocess"]](ctx)
    list(ctx = ctx, block = block, prep = prep)
}

.aa_fit_prepared <- function(ctx, block, prep) {
    out <- block[["edge_case"]](ctx, prep)  # check for K = 1 or K == N
    if (!is.null(out)) return(out)

    # Run multiple restarts and return the best fit
    # TODO: add parallelization option for multiple restarts
    best_fit <- NULL
    best_loss <- Inf
    for (.nrep_i in seq_len(ctx[["nrep"]])) {
        init_vars <- block[["init"]](ctx, prep)
        fit <- block[["fit"]](ctx, prep, init_vars)
        current_loss <- block[["final_loss"]](fit)
        if (is.finite(current_loss) && current_loss < best_loss) {
            best_fit <- fit
            best_loss <- current_loss
        }
    }
    block[["prepare_output"]](ctx, prep, best_fit)
}

.aa_fit_engine <- function(call,
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
    setup <- .aa_fit_engine_setup(
        call = call,
        x = x,
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
        data = data,
        ...
    )
    .aa_fit_prepared(setup[["ctx"]], setup[["block"]], setup[["prep"]])
}

.aa_final_loss <- function(fit) {
    fit[["loss"]][["loss"]][fit[["i"]] + 1L]
}

.aa_check_x <- function(x, missing = FALSE, validate_values = TRUE) {
    stopifnot("data must be a matrix-like object" = is_tabular(x))
    if (inherits(x, "sparseMatrix")) {
        if (!methods::is(x, "dMatrix"))
            stop("data must be numeric", call. = FALSE)
        values <- x@x
    } else {
        values <- as.matrix(x)
        stopifnot("data must be numeric" = is.numeric(values))
    }
    if (!validate_values)
        return(invisible(TRUE))
    if (missing) {
        stopifnot("observed data contains non-finite values" =
                      all(is.na(values) | is.finite(values)))
    } else {
        stopifnot("data contains missing values" = !any(is.na(values)))
        stopifnot("data contains non-finite values" = all(is.finite(values)))
    }
    invisible(TRUE)
}

.aa_check_missing_route <- function(ctx) {
    if (!identical(ctx[["method"]], "pgd") && ctx[["missing"]])
        stop("`missing = TRUE` is only supported for `method = 'pgd'`.", call. = FALSE)
    if (ctx[["missing"]] && !identical(ctx[["robust"]], FALSE))
        stop("`robust` is not supported with `missing = TRUE`.", call. = FALSE)
    if (ctx[["missing"]] && !is.null(ctx[["weights"]]))
        stop("`weights` are not supported with `missing = TRUE`.", call. = FALSE)
    if (ctx[["missing"]] && is_matrix(ctx[["scale"]]))
        stop("matrix `scale` is not supported with `missing = TRUE`.", call. = FALSE)
    invisible(TRUE)
}

.aa_check_max_no_update <- function(max_no_update) {
    stopifnot("max_no_update must be a positive integer" = is_count(max_no_update))
    invisible(TRUE)
}

.aa_check_projected_gradient_controls <- function(step_size,
                                                  max_iter_optimizer,
                                                  step_shrinkage,
                                                  max_no_update) {
    stopifnot("step_size must be positive" = is_positive(step_size))
    stopifnot("max_iter_optimizer must be a positive integer" = is_count(max_iter_optimizer))
    stopifnot("step_shrinkage must be between (0, 1)" =
                  is_number(step_shrinkage) && step_shrinkage > 0 && step_shrinkage < 1)
    .aa_check_max_no_update(max_no_update)
    invisible(TRUE)
}

.aa_euclidean_check <- function(ctx) {
    .aa_check_missing_route(ctx)
    .aa_check_fit_controls(ctx)
    .aa_check_scale(ctx[["scale"]], ncol(ctx[["x"]]))
}

.aa_euclidean_preprocess <- function(ctx, bigM = 0) {
    data <- ctx[["x"]]
    sd_threshold <- ctx[["sd_threshold"]]
    weights <- ctx[["weights"]]
    scale <- ctx[["scale"]]

    if (ctx[["verbose"]]) message("Preprocessing data...")

    if (ctx[["missing"]])
        return(.aa_preprocess_missing(data, sd_threshold, ctx[["verbose"]], scale = scale))

    if (inherits(data, "sparseMatrix"))
        data <- Matrix::drop0(data)

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

    X <- .aa_filter_low_variance(X, sd_threshold)
    mask <- attr(X, "mask")

    if (identical(scale_mode, "matrix")) {
        retained_names <- names(original_center)
        if (!is.null(mask)) {
            retained_names <- retained_names[mask]
            scale <- scale[mask, mask, drop = FALSE]
        }
        scale_factor <- t(chol(as.matrix(scale)))  # lower triangular
        X <- as.matrix(X) %*% scale_factor
        colnames(X) <- retained_names
        attr(X, "mask") <- mask
        # Compressed storage
        attr(X, "scale:factor") <- .aa_pack_lower_tri(scale_factor)
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
    N <- nrow(X)
    if (is.null(bigM))
        bigM <- .aa_auto_bigM(X)

    if (!is.null(weights)) {
        if (length(weights) != N) {
            fmt <- "Number of weights (%d) must equal number of rows in data (%d)"
            stop(sprintf(fmt, length(weights), N))
        }
        stopifnot("Weights contain NA values" = !any(is.na(weights)))
        stopifnot("Weights must be non-negative" = all(weights >= 0))
        stopifnot("at least one weight must be positive" = any(weights > 0))
        if (any(weights == 0))
            warning("Some sample weights are zero.", call. = FALSE)
        weights <- weights / mean(weights)
        x_attrs <- attributes(X)
        X <- X * sqrt(weights)  # sqrt will be undone during square loss computation
        attributes(X) <- utils::modifyList(attributes(X), x_attrs)
        attr(X, "weights") <- weights
    }

    if (bigM > 0) {
        x_attrs <- attributes(X)
        bigM_col <- matrix(bigM, nrow = N, ncol = 1L, dimnames = list(rownames(X), "bigM"))
        if (inherits(X, "sparseMatrix"))
            bigM_col <- as(bigM_col, "sparseMatrix")
        X <- cbind(bigM_col, X)
        attr(X, "scaled:center") <- x_attrs[["scaled:center"]]
        attr(X, "scaled:scale") <- x_attrs[["scaled:scale"]]
        attr(X, "scale:mode") <- x_attrs[["scale:mode"]]
        attr(X, "scale:factor") <- x_attrs[["scale:factor"]]
        attr(X, "restore:center") <- x_attrs[["restore:center"]]
        attr(X, "mask") <- x_attrs[["mask"]]
        attr(X, "bigM") <- 1L
        attr(X, "bigM.value") <- bigM
    }

    list(X = X)
}

.aa_euclidean_edge_case <- function(ctx, prep) {
    out <- .aa_checks_edge_cases(
        prep[["X"]],
        ctx[["K"]],
        ctx[["verbose"]],
        M = prep[["M"]]
    )
    if (is.null(out)) return(NULL)
    .aa_euclidean_output(
        ctx,
        prep,
        list(
            A0 = out[["init"]],
            A = coordinates(out),
            B = out[["coefficients"]],
            S = out[["compositions"]],
            delta = out[["slack"]],
            i = nrow(out[["loss"]]) - 1L,
            loss = out[["loss"]],
            converged = out[["converged"]]
        ),
        fit_info = list(method = ctx[["method"]])
    )
}

.aa_euclidean_init <- function(ctx, prep, delta = 0) {
    X <- prep[["X"]]
    init <- ctx[["init"]]
    init_args <- ctx[["init_args"]]
    if (ctx[["verbose"]]) message("Initializing archetypes...")
    L <- ctx[["max_iter"]] + 1L

    if (is_tabular(init)) {
        init <- .aa_preprocess_init(init, X)
        if (length(init_args) > 0L) {
            warning("`init_args` are ignored when `init` is a matrix", call. = FALSE)
            init_args <- list()
        }
        return(.aa_matrix_init(X, ctx[["K"]], init, ctx[["eps"]], L, delta))
    }

    if (is_non_empty_string(init)) {
        init_args <- c(list(method = init), init_args)
        init <- aa_init
    } else if (!is.function(init)) {
        stop("`init` must be a function, a single non-empty string, or archetypes coordinate matrix")
    }

    init_vars <- do.call(init, args = c(list(X = X, K = ctx[["K"]]), init_args))
    rownames(init_vars[["A"]]) <- .aa_init_names(init_vars[["A"]])
    rownames(init_vars[["B"]]) <- rownames(init_vars[["A"]])
    init_vars[["S"]] <- .aa_init_S(X, init_vars[["A"]], eps = ctx[["eps"]])
    init_vars[["loss"]] <- list(loss = rep(NA_real_, L), r2 = rep(NA_real_, L))
    init_vars
}

.aa_euclidean_output <- function(ctx, prep, fit, fit_info = list()) {
    X <- if (!is.null(prep[["output_X"]])) prep[["output_X"]] else prep[["X"]]
    feature_map <- .aa_euclidean_feature_map(X)
    feature_map[["eps"]] <- ctx[["eps"]]

    j <- fit[["i"]] + 1L
    loss <- as.data.frame(fit[["loss"]])[seq_len(j), , drop = FALSE]
    rownames(loss) <- NULL

    A_feature <- .aa_drop_bigM_column(fit[["A"]], attr(X, "bigM"))
    A <- fit[["A"]]
    A0 <- fit[["A0"]]
    archetype_names <- if (!is.null(A0)) rownames(A0) else rownames(A)
    A <- .aa_feature_map_inverse(feature_map, A_feature)
    if (!is.null(A0))
        A0 <- .aa_feature_map_inverse(feature_map, .aa_drop_bigM_column(A0, attr(X, "bigM")))

    family <- prep[["family"]] %||% ctx[["family"]] %||% "gaussian"

    if (is.null(archetype_names))
        archetype_names <- paste0("A", seq_len(nrow(A)))
    rownames(A) <- rownames(fit[["B"]]) <- colnames(fit[["S"]]) <- archetype_names
    rownames(A_feature) <- archetype_names
    if (!is.null(A0))
        rownames(A0) <- archetype_names
    colnames(fit[["B"]]) <- rownames(fit[["S"]]) <- rownames(X)

    if (!fit[["converged"]]) {
        fmt <- "Algorithm did not converge after %d iterations"
        warning(sprintf(fmt, fit[["i"]]), call. = FALSE)
    }

    if (ctx[["verbose"]]) {
        fmt <- ifelse(fit[["converged"]],
                      "Converged after %d iterations:",
                      "Final iteration %d:")
        fmt <- paste(fmt, "loss = %.4g, R2 = %.3f")
        message(sprintf(fmt, fit[["i"]], loss[j, "loss"], loss[j, "r2"]))
    }

    init <- ctx[["init"]]
    if (!is.null(init) && !is.character(init))
        init <- if (is.function(init)) "function" else "matrix"
    fit_info <- fit_info %|p|% list(
        family = prep[["family"]] %||% ctx[["family"]] %||% "gaussian",
        robust = !identical(ctx[["robust"]], FALSE),
        robust_psi = .aa_robust_label(ctx[["robust"]]) %||% NA_character_,
        robust_args = I(list(ctx[["robust_args"]])),
        missing = isTRUE(ctx[["missing"]]),
        delta = fit[["delta"]] %||% 0,
        init = init,
        scaling = .aa_scale_label(ctx[["scale"]]),
        sample_weights = !is.null(ctx[["weights"]])
    )

    archetypes(
        call         = ctx[["call"]],
        data         = ctx[["x"]],
        weights      = if (!is.null(fit[["row_weights"]])) fit[["row_weights"]] else attr(prep[["X"]], "weights"),
        init         = A0,
        coefficients = fit[["B"]],
        compositions = fit[["S"]],
        slack        = fit[["delta"]] %||% 0,
        loss         = loss,
        converged    = fit[["converged"]],
        family       = family,
        fit_info     = fit_info,
        feature_map  = feature_map,
        A            = A_feature
    )
}

.aa_scale_label <- function(scale) {
    if (identical(scale, FALSE))
        return("none")
    if (isTRUE(scale))
        return("z-score")
    if (is.numeric(scale) && is.null(dim(scale)))
        return("custom")
    if (is_matrix(scale))
        return("metric")
    "custom"
}
