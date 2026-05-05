#' Run Archetypal Analysis
#'
#' `run_aa()` is the common entry point for fitting archetypal analysis models.
#' It handles shared validation, preprocessing, initialization, and output
#' formatting, then delegates the optimization loop to the selected method.
#'
#' @param data data matrix (rows = samples, columns = dimensions)
#' @param K number of archetypes
#' @param method fitting method. One of `"pgd"` or `"nnls"` (default: `"pgd"`).
#' @param init function, method string, or numeric coordinate matrix to initialize
#'   archetypes (default: `"furthest_sum"`). When a matrix is supplied it must
#'   have dimension `K x ncol(data)`. Rows outside the allowed data hull are
#'   projected into it with a warning; row names, when present, are used as
#'   archetype names.
#' @param init_args list of additional arguments for the initialization function
#' @param weights optional vector of sample weights (default: NULL)
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
#'   `pseudo_pgd`, `step_size`, `max_iter_optimizer`, and `step_shrinkage`.
#'   For `"nnls"`, these are `ols_solver` and `bigM`.
#'
#' @returns An object of class \code{\link{archetypes}}
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "YAAAP"))
#' run_aa(as.matrix(toy), K = 3)
#' run_aa(as.matrix(toy), K = 3, method = "nnls")
#'
#' @export
run_aa <- function(data,
                   K,
                   method = c("pgd", "nnls"),
                   init = "furthest_sum",
                   init_args = list(),
                   weights = NULL,
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
    # Keep the shared implementation in `.aa_run()` so compatibility wrappers
    # can pass their own `match.call()` into the fitted object.
    .aa_run(
        call = match.call(),
        data = data,
        K = K,
        method = method,
        init = init,
        init_args = init_args,
        weights = weights,
        robust = robust,
        tukey_c = tukey_c,
        sd_threshold = sd_threshold,
        max_iter = max_iter,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        verbose = verbose,
        method_args = list(...)
    )
}

.aa_run <- function(call,
                    data,
                    K,
                    method,
                    init,
                    init_args,
                    weights,
                    robust,
                    tukey_c,
                    sd_threshold,
                    max_iter,
                    tol,
                    tol_r2,
                    max_kappa,
                    eps,
                    verbose,
                    method_args) {
    .aa_check_inputs( # nolint: object_usage_linter.
        data = data,
        K = K,
        tol = tol,
        tol_r2 = tol_r2,
        max_kappa = max_kappa,
        eps = eps,
        robust = robust,
        tukey_c = tukey_c
    )
    method <- match.arg(method, c("pgd", "nnls"))

    out <- .aa_checks_edge_cases(data, K, verbose)
    if (!is.null(out)) return(out)

    method_config <- switch(
        method,
        pgd = {
            args <- do.call(.aa_pgd_method_args, method_args)
            list(
                bigM = 0,
                delta = args[["delta"]],
                fit_fun = .aa_fit_pgd,
                fit_args = list(
                    robust = robust,
                    delta = args[["delta"]],
                    pseudo_pgd = args[["pseudo_pgd"]],
                    step_size = args[["step_size"]],
                    max_iter_optimizer = args[["max_iter_optimizer"]],
                    step_shrinkage = args[["step_shrinkage"]]
                )
            )
        },
        nnls = {
            args <- do.call(.aa_nnls_method_args, method_args)
            list(
                bigM = args[["bigM"]],
                delta = 0,
                fit_fun = .aa_fit_nnls,
                fit_args = list(
                    ols_solver = args[["ols_solver"]]
                )
            )
        }
    )

    weight_fun <- if (robust) {
        function(row_rss) .aa_bisquare_weights(row_rss, c = tukey_c)
    } else {
        .aa_unit_weights
    }

    pre <- .aa_preprocess(data, sd_threshold, weights, verbose, bigM = method_config[["bigM"]])
    X <- pre[["X"]]
    undo_scale <- pre[["undo_scale"]]
    rm(pre)

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
                                step_shrinkage = 0.5) {
    stopifnot("step_size must be positive" = step_size > 0)
    stopifnot("step_shrinkage must be between (0, 1)" =
                  step_shrinkage > 0 && step_shrinkage < 1)
    stopifnot("delta must be single non-negative number" =
                  length(delta) == 1 && delta >= 0)

    list(
        delta = delta,
        pseudo_pgd = pseudo_pgd,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage
    )
}

.aa_nnls_method_args <- function(ols_solver = c("qr", "ginv", "BFGS"),
                                 bigM = 200) {
    ols_solver <- match.arg(ols_solver)
    stopifnot("`bigM` must be greater than 0" = bigM > 0)

    list(
        ols_solver = ols_solver,
        bigM = bigM
    )
}
