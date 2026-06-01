#' Archetypal Analysis Preprocessing Step for recipes
#'
#' `step_archetypes()` creates a *specification* of a recipe step that projects
#' numeric predictors onto a K-archetype simplex.  During `prep()` the
#' archetypes are fitted once (with optional random-restart via `nrep`); during
#' `bake()` new data are projected onto the learnt simplex.
#'
#' @param recipe A recipe object.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See [recipes::selections()] for more details.
#' @param role Role for the *new* columns produced by this step.
#'   Default `"predictor"`.
#' @param trained Logical. Has the step been trained (i.e. has `prep()` been
#'   called)?  Set by `prep()` — do not change manually.
#' @param num_comp Number of archetypes K (default `3L`). Tunable via
#'   [tune::tune()].
#' @param delta Convexity penalty for the pgd solver (default `0`). Ignored
#'   when `fit_method != "pgd"`. Tunable via [tune::tune()].
#' @param fit_method Solver to use: `"pgd"` (default), `"nnls"`, or `"paa"`.
#' @param options A named list of additional arguments passed to [run_aa()]
#'   (e.g. `nrep`, `max_iter`, `scale`). These are not validated here.
#' @param reconstruct If `TRUE`, `bake()` also appends columns named
#'   `rec_<original_col>` containing the low-rank reconstruction
#'   `S %*% Z` (compositions times archetype coordinates).
#' @param keep_original_cols Should the original predictor columns be retained
#'   after baking? Default `FALSE`.
#' @param res The fitted [archetypes][yaap::archetypes_pgd] object. `NULL`
#'   before `prep()`, set automatically by `prep()`.
#' @param col_names Character vector of selected column names. Set by `prep()`.
#' @param seed Integer seed for reproducible fitting. Evaluated once at
#'   construction time (default: a random integer).
#' @param skip Should the step be skipped during `bake()` on the test set?
#'   Default `FALSE`.
#' @param id A character string that identifies this step in the recipe.
#'
#' @return An updated recipe with the new step appended.
#'
#' @details
#' `step_archetypes()` is **not** supported for `method = "directional"` or
#' `method = "kernel"`.
#'
#' The `num_comp` and `delta` parameters support [tune::tune()] for
#' hyperparameter search.
#'
#' @examples
#' if (requireNamespace("recipes", quietly = TRUE)) {
#'     rec <- recipes::recipe(Species ~ ., data = iris) |>
#'         step_archetypes(
#'             recipes::all_numeric_predictors(),
#'             num_comp = 3L,
#'             options = list(max_iter = 20L, tol_r2 = 0.95)
#'         ) |>
#'         recipes::prep(training = iris)
#'     recipes::bake(rec, new_data = iris)
#' }
#'
#' @export
step_archetypes <- function(
  recipe,
  ...,
  role = "predictor",
  trained = FALSE,
  num_comp = 3L,
  delta = 0,
  fit_method = "pgd",
  options = list(),
  reconstruct = FALSE,
  keep_original_cols = FALSE,
  res = NULL,
  col_names = NULL,
  seed = sample.int(1e5L, 1L),
  skip = FALSE,
  id = recipes::rand_id("archetypes")
) {
    rlang::check_installed(c("recipes", "withr"), reason = "required by step_archetypes()")
    recipes::add_step(
        recipe,
        .aa_step_archetypes_new(
            terms              = rlang::enquos(...),
            role               = role,
            trained            = trained,
            num_comp           = num_comp,
            delta              = delta,
            fit_method         = fit_method,
            options            = options,
            reconstruct        = reconstruct,
            keep_original_cols = keep_original_cols,
            res                = res,
            col_names          = col_names,
            seed               = seed,
            skip               = skip,
            id                 = id
        )
    )
}

.aa_step_archetypes_new <- function(
  terms, role, trained, num_comp, delta, fit_method, options,
  reconstruct, keep_original_cols, res, col_names, seed, skip, id
) {
    recipes::step(
        subclass           = "archetypes",
        terms              = terms,
        role               = role,
        trained            = trained,
        num_comp           = num_comp,
        delta              = delta,
        fit_method         = fit_method,
        options            = options,
        reconstruct        = reconstruct,
        keep_original_cols = keep_original_cols,
        res                = res,
        col_names          = col_names,
        seed               = seed,
        skip               = skip,
        id                 = id
    )
}

#' @exportS3Method recipes::prep
prep.step_archetypes <- function(x, training, info = NULL, ...) {
    col_names <- recipes::recipes_eval_select(x$terms, training, info)

    X <- as.matrix(training[, col_names, drop = FALSE])

    aa_args <- c(
        list(x = X, K = x$num_comp, method = x$fit_method),
        if (identical(x$fit_method, "pgd")) list(delta = x$delta),
        x$options
    )

    res <- withr::with_seed(x$seed, do.call(run_aa, aa_args))

    .aa_step_archetypes_new(
        terms              = x$terms,
        role               = x$role,
        trained            = TRUE,
        num_comp           = x$num_comp,
        delta              = x$delta,
        fit_method         = x$fit_method,
        options            = x$options,
        reconstruct        = x$reconstruct,
        keep_original_cols = recipes::get_keep_original_cols(x),
        res                = res,
        col_names          = col_names,
        seed               = x$seed,
        skip               = x$skip,
        id                 = x$id
    )
}

#' @exportS3Method recipes::bake
bake.step_archetypes <- function(object, new_data, ...) {
    recipes::check_new_data(object$col_names, object, new_data)

    X_new <- as.matrix(new_data[, object$col_names, drop = FALSE])
    S <- predict(object$res, newdata = X_new, type = "compositions") # N x K compositions
    arch_names <- anames(object$res)
    S_df <- as.data.frame(S)
    colnames(S_df) <- arch_names

    # check_name handles name conflicts (e.g. existing columns with same name)
    S_df <- recipes::check_name(S_df, new_data, object, newname = arch_names)
    new_data <- vctrs::vec_cbind(new_data, tibble::as_tibble(S_df), .name_repair = "minimal")

    if (isTRUE(object$reconstruct)) {
        X_hat <- predict(object$res, newdata = X_new, type = "reconstruction")
        rec_names <- paste0("rec_", object$col_names)
        rec_df <- as.data.frame(X_hat)
        colnames(rec_df) <- rec_names
        rec_df <- recipes::check_name(rec_df, new_data, object, newname = rec_names)
        new_data <- vctrs::vec_cbind(new_data, tibble::as_tibble(rec_df), .name_repair = "minimal")
    }

    new_data <- recipes::remove_original_cols(new_data, object, object$col_names)
    new_data
}

#' @exportS3Method base::print
print.step_archetypes <- function(x, width = max(20, options()$width - 29), ...) {
    title <- paste0(
        "Archetypal Analysis Projection (K=", x$num_comp,
        if (identical(x$fit_method, "pgd")) paste0(", delta=", x$delta),
        ") "
    )
    recipes::print_step(x$col_names, x$terms, x$trained, title, width)
    invisible(x)
}

#' @exportS3Method generics::tidy
tidy.step_archetypes <- function(x, matrix = "coordinates", ...) {
    if (!recipes::is_trained(x)) {
        res <- tibble::tibble(
            archetype = character(0L),
            term      = character(0L),
            value     = double(0L),
            id        = x$id
        )
        return(res)
    }
    out <- tidy.archetypes(x$res, matrix = matrix, ...)
    out$id <- x$id
    out
}

#' @exportS3Method tune::tunable
tunable.step_archetypes <- function(x, ...) {
    tibble::tibble(
        name = c("num_comp", "delta"),
        call_info = list(
            list(pkg = "dials", fun = "num_comp", range = c(1L, 10L)),
            list(pkg = "dials", fun = "mixture", range = c(0, 1))
        ),
        source = "recipe",
        component = "step_archetypes",
        component_id = x$id
    )
}

#' @exportS3Method recipes::required_pkgs
required_pkgs.step_archetypes <- function(x, ...) {
    c("yaap", "withr")
}
