# broom tidiers for archetypes and kernel_archetypes objects.
# Registers tidy(), glance(), and augment() generics from the `generics`
# package so dispatch works without requiring the user to attach `broom`.

# -------------------------------------------------------------------------
# tidy
# -------------------------------------------------------------------------

#' Tidy an archetypes model
#'
#' Converts an `archetypes` model into a tidy long-form tibble.
#'
#' @param x An object of class `archetypes` or `kernel_archetypes`.
#' @param matrix Which component to return: `"coordinates"` (K x M archetype
#'   coordinates, default), `"coefficients"` (K x N weights of archetypes over
#'   samples), or `"compositions"` (N x K weights of samples over archetypes).
#' @param ... Ignored.
#'
#' @return
#' A tibble. Column names depend on `matrix`:
#' * `"coordinates"`: `archetype`, `term`, `value` (K * M rows)
#' * `"coefficients"`: `archetype`, `sample`, `value` (K * N rows)
#' * `"compositions"`: `sample`, `archetype`, `value` (N * K rows)
#'
#' @seealso [glance.archetypes()], [augment.archetypes()]
#'
#' @exportS3Method generics::tidy
tidy.archetypes <- function(x, matrix = c("coordinates", "coefficients", "compositions"), ...) {
    matrix <- match.arg(matrix)
    switch(matrix,
        coordinates = {
            A <- x[["coordinates"]]
            .aa_pivot_matrix(
                mat      = A,
                row_name = "archetype",
                col_name = "term",
                row_ids  = rownames(A),
                col_ids  = colnames(A),
                row_fill = paste0("A", seq_len(nrow(A))),
                col_fill = paste0("X", seq_len(ncol(A)))
            )
        },
        coefficients = {
            B <- x[["coefficients"]]
            .aa_pivot_matrix(
                mat      = B,
                row_name = "archetype",
                col_name = "sample",
                row_ids  = rownames(B),
                col_ids  = colnames(B),
                row_fill = paste0("A", seq_len(nrow(B))),
                col_fill = paste0("S", seq_len(ncol(B)))
            )
        },
        compositions = {
            S <- x[["compositions"]]
            .aa_pivot_matrix(
                mat      = S,
                row_name = "sample",
                col_name = "archetype",
                row_ids  = rownames(S),
                col_ids  = colnames(S),
                row_fill = paste0("S", seq_len(nrow(S))),
                col_fill = paste0("A", seq_len(ncol(S)))
            )
        }
    )
}

#' @rdname tidy.archetypes
#' @details
#' For `kernel_archetypes`, `matrix = "coordinates"` returns the
#' `coordinates` matrix (`coefficients %*% data`) when available, and
#' emits a warning with an empty tibble otherwise. `"coefficients"` and
#' `"compositions"` behave identically to the `archetypes` method.
#'
#' @exportS3Method generics::tidy
tidy.kernel_archetypes <- function(x, matrix = c("coordinates", "coefficients", "compositions"), ...) {
    matrix <- match.arg(matrix)
    if (identical(matrix, "coordinates")) {
        A <- x[["coordinates"]]
        if (is.null(A)) {
            warning(paste(
                "`coordinates` is NULL for this kernel_archetypes object.",
                "Supply `data` with `kernel = 'precomputed'` to get input-space coordinates.",
                "Returning an empty tibble."
            ), call. = FALSE)
            return(tibble::tibble(archetype = character(), term = character(), value = numeric()))
        }
        return(.aa_pivot_matrix(
            mat      = A,
            row_name = "archetype",
            col_name = "term",
            row_ids  = rownames(A),
            col_ids  = colnames(A),
            row_fill = paste0("A", seq_len(nrow(A))),
            col_fill = paste0("X", seq_len(ncol(A)))
        ))
    }
    # coefficients and compositions: identical structure to archetypes
    tidy.archetypes(x, matrix = matrix, ...)
}

# -------------------------------------------------------------------------
# glance
# -------------------------------------------------------------------------

#' Glance at an archetypes model
#'
#' Returns a one-row tibble of model-level summary statistics.
#'
#' @param x An object of class `archetypes` or `kernel_archetypes`.
#' @param ... Ignored.
#'
#' @return
#' A one-row tibble with columns:
#' * `K` — number of archetypes
#' * `converged` — did the optimizer converge?
#' * `loss` — final residual sum of squares
#' * `r2` — final R²
#' * `n_iter` — number of iterations run (excluding initialisation row)
#' * `family` — observation family (`archetypes` only)
#'
#' @seealso [tidy.archetypes()], [augment.archetypes()]
#'
#' @exportS3Method generics::glance
glance.archetypes <- function(x, ...) {
    loss_df <- x[["loss"]]
    n       <- nrow(loss_df)
    tibble::tibble(
        K         = nrow(x[["coordinates"]]),
        converged = isTRUE(x[["converged"]]),
        loss      = loss_df[["loss"]][n],
        r2        = loss_df[["r2"]][n],
        n_iter    = n - 1L,
        aic       = tryCatch(AIC(x), error = function(e) NA_real_, warning = function(w) NA_real_),
        family    = x[["family"]] %||% "gaussian"
    )
}

#' @rdname glance.archetypes
#' @exportS3Method generics::glance
glance.kernel_archetypes <- function(x, ...) {
    loss_df <- x[["loss"]]
    n       <- nrow(loss_df)
    tibble::tibble(
        K         = nrow(x[["coefficients"]]),
        converged = isTRUE(x[["converged"]]),
        loss      = loss_df[["loss"]][n],
        r2        = loss_df[["r2"]][n],
        n_iter    = n - 1L
    )
}

# -------------------------------------------------------------------------
# augment
# -------------------------------------------------------------------------

#' Augment data with composition weights from an archetypes model
#'
#' Adds per-sample archetype composition columns to the original data.
#'
#' @param x An object of class `archetypes` or `kernel_archetypes`.
#' @param data Optional data frame or matrix to augment. If `NULL`, uses the
#'   data stored inside `x` (if available). For `archetypes` objects, passing
#'   new data triggers [predict.archetypes()] to compute compositions.
#' @param ... Passed to [predict.archetypes()] when `data` is provided and
#'   `x` is an `archetypes` object.
#'
#' @return
#' A tibble with all columns from `data` plus one column per archetype named
#' `.A1`, `.A2`, etc. (dot-prefixed `anames(x)`).
#'
#' @seealso [tidy.archetypes()], [glance.archetypes()]
#'
#' @exportS3Method generics::augment
augment.archetypes <- function(x, data = NULL, ...) {
    if (is.null(data)) {
        stored <- x[["data"]]
        if (is.null(stored))
            stop(paste(
                "Data must be provided either as an argument to `augment()` or",
                "when constructing the archetypes object."
            ), call. = FALSE)
        out <- tibble::as_tibble(stored)
        S   <- x[["compositions"]]
    } else {
        out <- tibble::as_tibble(data)
        S   <- predict(x, newdata = data, ...)
    }
    comp_names <- paste0(".", anames(x))
    comp_df    <- tibble::as_tibble(S, .name_repair = "minimal")
    colnames(comp_df) <- comp_names
    tibble::add_column(out, comp_df)
}

#' @rdname augment.archetypes
#' @details
#' For `kernel_archetypes`, `data` is used only to convert the stored data to
#' a tibble; compositions always come from the stored `x[["compositions"]]`
#' since projecting new samples requires the original Gram matrix.
#'
#' @exportS3Method generics::augment
augment.kernel_archetypes <- function(x, data = NULL, ...) {
    raw <- if (!is.null(data)) data else x[["data"]]
    if (is.null(raw))
        stop(paste(
            "Data must be provided either as an argument to `augment()` or",
            "when constructing the kernel_archetypes object."
        ), call. = FALSE)
    out        <- tibble::as_tibble(raw)
    S          <- x[["compositions"]]
    comp_names <- paste0(".", anames(x))
    comp_df    <- tibble::as_tibble(S, .name_repair = "minimal")
    colnames(comp_df) <- comp_names
    tibble::add_column(out, comp_df)
}

# -------------------------------------------------------------------------
# Internal helper
# -------------------------------------------------------------------------

# Pivot a matrix to a long tibble with (row_name, col_name, value) columns.
# NULL row/column names are replaced by row_fill / col_fill.
.aa_pivot_matrix <- function(mat, row_name, col_name, row_ids, col_ids,
                             row_fill, col_fill) {
    row_ids <- row_ids %||% row_fill
    col_ids <- col_ids %||% col_fill
    # R stores matrices column-major, so as.vector() iterates columns first.
    # rep(..., times = ncol) cycles the K row labels M times  (A1,A2,..,AK repeated M times)
    # rep(..., each = nrow)  repeats each column label K times (X1,X1,..,X1, X2,..)
    out <- tibble::tibble(
        row = rep(row_ids, times = ncol(mat)),
        col = rep(col_ids, each  = nrow(mat)),
        value = as.vector(mat)
    )
    colnames(out)[1:2] <- c(row_name, col_name)
    out
}
