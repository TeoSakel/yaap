# Directional Archetypes Class -----------------------------------------------

#' Internal directional archetype analysis constructor
#'
#' Builds a `directional_archetypes` object from fitted directional AA matrices.
#' `coordinates` must be K x M, `coefficients` K x N, and `compositions` N x K.
#' The base `archetypes()` constructor enforces stochasticity and shared object
#' invariants; this constructor adds generator-space metadata for directional methods.
#'
#' @param coordinates Raw directional archetype coordinates.
#' @param coefficients Archetype generator coefficients.
#' @param compositions Sample compositions.
#' @param loss Per-iteration directional loss metrics.
#' @param converged Logical convergence flag.
#' @param call Matched call.
#' @param data Original data matrix.
#' @param init Initial archetype coordinates.
#' @param generator_data Row-normalized and possibly hemisphere-flipped data
#'   used to generate archetypes.
#' @param hemisphere_direction Dominant hemisphere direction, or `NULL`.
#' @param row_norms Original data row norms.
#' @param precision Precision mode used for the Watson loss.
#' @param fit_info Internal fit metadata.
#' @param feature_map Internal prediction feature-map metadata.
#'
#' @noRd
.aa_new_directional_archetypes <- function(coordinates,
                                           coefficients,
                                           compositions,
                                           loss = NULL,
                                           converged = TRUE,
                                           call = NULL,
                                           data = NULL,
                                           init = NULL,
                                           generator_data = NULL,
                                           hemisphere_direction = NULL,
                                           row_norms = NULL,
                                           precision = NULL,
                                           fit_info = list(),
                                           feature_map = NULL) {
    out <- archetypes(
        A = coordinates,
        coefficients = coefficients,
        compositions = compositions,
        slack = 0,
        loss = loss,
        converged = converged,
        call = call,
        data = data,
        init = init,
        family = "watson",
        fit_info = fit_info,
        feature_map = feature_map
    )
    out[["generator_data"]] <- generator_data
    out[["hemisphere_direction"]] <- hemisphere_direction
    out[["row_norms"]] <- row_norms
    out[["precision"]] <- precision
    out[["directions"]] <- .aa_unit_rows(coordinates)
    class(out) <- c("directional_archetypes", class(out))
    out
}

#' @exportS3Method
fitted.directional_archetypes <- function(object, ...) {
    Y <- compositions(object) %*% .aa_input_coordinates_matrix(object)
    Y <- .aa_unit_rows(Y)
    X <- object[["data"]]
    if (!is.null(X)) {
        X <- .aa_unit_rows(as.matrix(X))
        Y <- .aa_align_rows(Y, X)
    }
    Y
}

#' @exportS3Method
residuals.directional_archetypes <- function(object, ...) {
    X <- object[["data"]]
    if (is.null(X)) {
        stop("Original data is missing from `object$data`.", call. = FALSE)
    }

    X <- .aa_unit_rows(as.matrix(X))
    Y <- .aa_unit_rows(compositions(object) %*% .aa_input_coordinates_matrix(object))
    Y <- .aa_align_rows(Y, X)
    X - Y
}

#' Predict compositions or reconstructions for directional archetypes
#'
#' Projects new directional samples onto the fitted archetype space, like
#' [predict.archetypes()]. Directional reconstructions are returned as unit-length directions;
#' `type = "compositions"` returns the archetype weights.
#'
#' @param object An object of class `directional_archetypes`.
#' @param newdata New directional data matrix.
#' @param type Prediction output type. `"reconstruction"` (default) or
#'   `"compositions"`.
#' @param max_iter,eps,step_size,max_iter_optimizer,step_shrinkage Optimization
#'   controls for composition fitting.
#' @param ... Ignored.
#'
#' @return A composition matrix (N_new x K) when `type = "compositions"`, or
#' a directional reconstruction matrix (N_new x M) when
#' `type = "reconstruction"`.
#'
#' @exportS3Method
predict.directional_archetypes <- function(object,
                                           newdata,
                                           type = c("reconstruction", "compositions"),
                                           max_iter = 100L,
                                           eps = 1e-8,
                                           step_size = 1.0,
                                           max_iter_optimizer = 10L,
                                           step_shrinkage = 0.5,
                                           ...) {
    type <- match.arg(type)
    # Prediction fixes the learned directional archetypes A and solves only for
    # S, the simplex composition of each new row in that archetype hull. This is
    # the same S subproblem used in the alternating training loop, but without
    # the B/C update because new data must not change the fitted archetypes.
    A <- .aa_input_coordinates_matrix(object)
    X <- as.matrix(newdata)
    if (ncol(X) != ncol(A)) {
        fmt <- "`newdata` has %d columns but `object$coordinates` has %d columns"
        stop(sprintf(fmt, ncol(X), ncol(A)))
    }
    .aa_check_no_zero_rows(X)
    X_loss <- if (identical(object[["precision"]], "unit")) .aa_unit_rows(X) else X
    S <- .aa_directional_fit_S(
        X_loss = X_loss,
        A = A,
        max_iter = max_iter,
        eps = eps,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage
    )
    if (identical(type, "compositions")) {
        return(S)
    }

    Y <- .aa_unit_rows(S %*% A)
    .aa_align_rows(Y, .aa_unit_rows(X))
}

#' @exportS3Method
AIC.directional_archetypes <- function(object, ...) {
    stop("AIC is not defined for Watson-loss directional archetypes.", call. = FALSE)
}
