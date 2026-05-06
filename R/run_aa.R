#' Run Archetypal Analysis
#'
#' `run_aa()` is the common entry point for fitting archetypal analysis models.
#' It handles shared validation, preprocessing, initialization, and output
#' formatting, then delegates the optimization loop to the selected method.
#'
#' @param data data matrix (rows = samples, columns = dimensions), or an
#'   object with a class-specific `run_aa()` method.
#' @param K number of archetypes
#' @param ... arguments passed to methods.
#'
#' @returns An object of class \code{\link{archetypes}}
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "YAAAP"))
#' run_aa(as.matrix(toy), K = 3)
#' run_aa(as.matrix(toy), K = 3, method = "nnls")
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
#' arch_fd <- coordinates_fd(fit_fd)
#' }
#'
#' @export
run_aa <- function(data, K, ...) {
    UseMethod("run_aa")
}

#' @rdname run_aa
#' @param method fitting method. One of `"pgd"` or `"nnls"` (default: `"pgd"`).
#' @param init function, method string, or numeric coordinate matrix to initialize
#'   archetypes (default: `"furthest_sum"`). When a matrix is supplied it must
#'   have dimension `K x ncol(data)`. Rows outside the allowed data hull are
#'   projected into it with a warning; row names, when present, are used as
#'   archetype names.
#' @param init_args list of additional arguments for the initialization function
#' @param weights optional vector of sample weights (default: NULL)
#' @param scale scaling or metric embedding used before fitting. `TRUE` applies
#'   the default z-score preprocessing, `FALSE` leaves columns on their original
#'   scale, a positive numeric vector divides columns by user-supplied scale
#'   factors, and a symmetric positive-definite matrix applies the corresponding
#'   feature metric embedding in the original data column space.
#' @param robust whether to use Tukey bisquare row reweighting (default: FALSE)
#' @param tukey_c tuning constant for Tukey bisquare weights (default: 4.685)
#' @param sd_threshold threshold for feature standard deviation to filter
#'   low-variance features (default: 1e-6)
#' @param max_iter maximum number of iterations (default: 100)
#' @param tol convergence tolerance based on residual sum of squares (default: 1e-6)
#' @param tol_r2 convergence tolerance based on R^2 (default: 0.9999)
#' @param max_kappa maximum condition number for archetypes (default: 1000)
#' @param eps small positive number to ensure numerical stability
#'   (default: 0 for sparse input 1e-8 for dense)
#' @param verbose whether to print progress messages (default: FALSE)
#' @param ... method-specific arguments. For `"pgd"`, these are `delta`,
#'   `pseudo_pgd`, `step_size`, `max_iter_optimizer`, `step_shrinkage`, and
#'   `max_no_update`. For `"nnls"`, these are `ols_solver`, `bigM`, and
#'   `max_no_update`.
#'
#' @exportS3Method
run_aa.default <- function(data,
                           K,
                           method = c("pgd", "nnls"),
                           init = "furthest_sum",
                           init_args = list(),
                           weights = NULL,
                           scale = TRUE,
                           robust = FALSE,
                           tukey_c = 4.685,
                           sd_threshold = 1e-6,
                           max_iter = 100L,
                           tol = 1e-6,
                           tol_r2 = 0.9999,
                           max_kappa = 1000,
                           eps = ifelse(inherits(data, "sparseMatrix"), 0, 1e-8),
                           verbose = FALSE,
                           ...) {
    call <- match.call()
    call[[1L]] <- quote(run_aa)

    .aa_run_aa_default(
        call = call,
        data = data,
        K = K,
        method = method,
        init = init,
        init_args = init_args,
        weights = weights,
        scale = scale,
        robust = robust,
        tukey_c = tukey_c,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        ...
    )
}

#' @rdname run_aa
#' @details
#' `run_aa.fd()` fits archetypal analysis to the coefficient matrix of an
#' `fda::fd` object. The feature scaling is set to the basis inner-product
#' matrix from `fda::eval.penalty(data$basis, Lfdobj = 0)`, so
#' `scale` cannot be supplied to this method. The fitted
#' `coordinates` remain a coefficient matrix; use [coordinates_fd()] to convert
#' them to an `fda::fd` object on demand.
#'
#' @exportS3Method
run_aa.fd <- function(data, K, ...) {
    if (!requireNamespace("fda", quietly = TRUE))
        stop("Package `fda` is required for `run_aa.fd()`.", call. = FALSE)

    method_args <- list(...)
    if ("scale" %in% names(method_args))
        stop("`scale` is computed from the fd basis and cannot be supplied to `run_aa.fd()`.",
             call. = FALSE)

    coefs <- stats::coef(data)
    if (length(dim(coefs)) != 2L)
        stop("`run_aa.fd()` currently supports fd objects with a 2-D coefficient matrix.",
             call. = FALSE)

    B <- t(coefs)
    G <- fda::eval.penalty(data[["basis"]], Lfdobj = 0)

    call <- match.call()
    call[[1L]] <- quote(run_aa)

    fit <- .aa_run_aa_default(
        call = call,
        data = B,
        K = K,
        scale = G,
        ...
    )

    fit[["data"]] <- data
    fit
}

.aa_run_aa_default <- function(call,
                               data,
                               K,
                               method = c("pgd", "nnls"),
                               init = "furthest_sum",
                               init_args = list(),
                               weights = NULL,
                               scale = TRUE,
                               robust = FALSE,
                               tukey_c = 4.685,
                               sd_threshold = 1e-6,
                               max_iter = 100L,
                               tol = 1e-6,
                               tol_r2 = 0.9999,
                               max_kappa = 1000,
                               eps = ifelse(inherits(data, "sparseMatrix"), 0, 1e-8),
                               verbose = FALSE,
                               ...) {
    .aa_check_inputs( # nolint: object_usage_linter.
        data = data,
        K = K,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        robust = robust,
        tukey_c = tukey_c,
        scale = scale
    )
    method <- match.arg(method, c("pgd", "nnls"))

    method_config <- switch(
        method,
        pgd = {
            args <- .aa_pgd_method_args(...)
            list(
                bigM = 0,
                delta = args[["delta"]],
                fit_fun = .aa_fit_pgd,
                fit_args = list(
                    delta = args[["delta"]],
                    pseudo_pgd = args[["pseudo_pgd"]],
                    step_size = args[["step_size"]],
                    max_iter_optimizer = args[["max_iter_optimizer"]],
                    step_shrinkage = args[["step_shrinkage"]],
                    max_no_update = args[["max_no_update"]]
                )
            )
        },
        nnls = {
            args <- .aa_nnls_method_args(...)
            list(
                bigM = args[["bigM"]],
                delta = 0,
                fit_fun = .aa_fit_nnls,
                fit_args = list(
                    ols_solver = args[["ols_solver"]],
                    max_no_update = args[["max_no_update"]]
                )
            )
        }
    )

    weight_fun <- if (robust) {
        function(row_rss) .aa_bisquare_weights(row_rss, c = tukey_c)
    } else {
        NULL
    }

    pre <- .aa_preprocess(
        data,
        sd_threshold,
        weights,
        verbose,
        bigM = method_config[["bigM"]],
        scale = scale
    )
    X <- pre[["X"]]
    undo_scale <- pre[["undo_scale"]]
    rm(pre)

    out <- .aa_checks_edge_cases(X, K, verbose)
    if (!is.null(out)) {
        return(.aa_prepare_output(
            call = call,
            data = data,
            X = X,
            A0 = out[["init"]],
            A = out[["coordinates"]],
            B = out[["coefficients"]],
            S = out[["compositions"]],
            delta = out[["slack"]],
            i = nrow(out[["loss"]]) - 1L,
            loss = out[["loss"]],
            converged = out[["converged"]],
            undo_scale = undo_scale,
            max_iter = max_iter,
            verbose = verbose
        ))
    }

    common_args <- list(
        X = X,
        weight_fun = weight_fun,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose
    )
    init_vars <- .aa_init_vars( # nolint: object_usage_linter.
        X = X,
        K = K,
        init = init,
        init_args = init_args,
        eps = eps,
        max_iter = max_iter,
        verbose = verbose,
        delta = method_config[["delta"]]
    )
    fit <- do.call(
        what = method_config[["fit_fun"]],
        args = c(common_args, init_vars, method_config[["fit_args"]])
    )

    .aa_prepare_output(
        call = call,
        data = data,
        X = X,
        A0 = fit[["A0"]],
        A = fit[["A"]],
        B = fit[["B"]],
        S = fit[["S"]],
        delta = fit[["delta"]],
        i = fit[["i"]],
        loss = fit[["loss"]],
        converged = fit[["converged"]],
        undo_scale = undo_scale,
        max_iter = max_iter,
        verbose = verbose
    )
}

.aa_pgd_method_args <- function(delta = 0,
                                pseudo_pgd = TRUE,
                                step_size = 1.0,
                                max_iter_optimizer = 10L,
                                step_shrinkage = 0.5,
                                max_no_update = 5L) {
    stopifnot("step_size must be positive" = step_size > 0)
    stopifnot("step_shrinkage must be between (0, 1)" =
                  step_shrinkage > 0 && step_shrinkage < 1)
    stopifnot("delta must be single non-negative number" =
                  length(delta) == 1 && delta >= 0)
    stopifnot("max_no_update must be a positive integer" =
                  max_no_update == as.integer(max_no_update) && max_no_update >= 1L)

    list(
        delta = delta,
        pseudo_pgd = pseudo_pgd,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage,
        max_no_update = as.integer(max_no_update)
    )
}

.aa_nnls_method_args <- function(ols_solver = c("qr", "ginv", "BFGS"),
                                 bigM = NULL,
                                 max_no_update = 5L) {
    ols_solver <- match.arg(ols_solver)
    if (!is.null(bigM)) {
        stopifnot("`bigM` must be NULL or a positive number" =
                      is.numeric(bigM) && length(bigM) == 1L && is.finite(bigM) && bigM > 0)
    }
    stopifnot("max_no_update must be a positive integer" =
                  max_no_update == as.integer(max_no_update) && max_no_update >= 1L)

    list(
        ols_solver = ols_solver,
        bigM = bigM,
        max_no_update = as.integer(max_no_update)
    )
}
