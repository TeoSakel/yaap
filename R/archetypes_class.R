# archetypes_class.R: Archetypes S3 class definition and methods

# Base Archetypes Class -------------------------------------------------------

# Internal constructor for archetypes S3 objects. Validates dimensions and
# stochasticity constraints, then wraps the result in the `archetypes` class.
archetypes <- function(A = NULL,
                       coefficients,
                       compositions,
                       slack = 0,
                       loss = NULL,
                       converged = TRUE,
                       call = NULL,
                       data = NULL,
                       init = NULL,
                       family = "gaussian",
                       weights = NULL,
                       fit_info = list(),
                       feature_map = NULL) {
    call <- call %||% match.call()
    loss <- loss %||% data.frame(loss = NA_real_, r2 = NA_real_)

    K <- if (!is.null(A)) nrow(A) else nrow(coefficients)
    N <- nrow(compositions)

    if (nrow(coefficients) != K) {
        fmt <- "nrow(coefficients) = %d does not match number of archetypes (%d)"
        stop(sprintf(fmt, nrow(coefficients), K))
    }
    if (ncol(coefficients) != N) {
        fmt <- "Inconsistent number of samples between compositions (%d) and coefficients (%d)"
        stop(sprintf(fmt, N, ncol(coefficients)))
    }
    if (ncol(compositions) != K) {
        fmt <- "ncol(compositions) = %d does not match number of archetypes (%d)"
        stop(sprintf(fmt, ncol(compositions), K))
    }

    if (is.null(feature_map)) {
        if (is.null(A)) {
            stop("`feature_map` must be supplied when `A` is NULL.", call. = FALSE)
        }
        feature_map <- .aa_identity_feature_map(A)
    }
    if (!inherits(feature_map, "aa_feature_map")) {
        stop("`feature_map` must be an `aa_feature_map` object.", call. = FALSE)
    }

    if (!is.null(A) && nrow(A) != K) {
        fmt <- "nrow(A) = %d does not match number of archetypes (%d)"
        stop(sprintf(fmt, nrow(A), K), call. = FALSE)
    }

    archetype_order <- .aa_canonical_archetype_order(A, coefficients)
    if (!identical(archetype_order, seq_along(archetype_order))) {
        if (!is.null(A)) {
            A <- A[archetype_order, , drop = FALSE]
        }
        coefficients <- coefficients[archetype_order, , drop = FALSE]
        compositions <- compositions[, archetype_order, drop = FALSE]
        if (!is.null(init)) {
            init <- init[archetype_order, , drop = FALSE]
        }
        if (length(slack) == length(archetype_order)) {
            slack <- slack[archetype_order]
        }
    }

    archetype_names <- rownames(coefficients) %||% paste0("A", seq_len(K))
    rownames(coefficients) <- archetype_names
    colnames(compositions) <- archetype_names
    if (!is.null(A)) {
        rownames(A) <- archetype_names
    }
    if (!is.null(init)) {
        rownames(init) <- archetype_names
    }

    if (!is_all_non_negative(slack)) {
        stop("All `slack` values must be non-negative.", call. = FALSE)
    }
    if (any(slack > 0)) {
        slack_tol <- sqrt(.Machine$double.eps)
        a <- rowSums(coefficients)
        if (!all(a <= 1 + slack + slack_tol)) {
            stop("Some `coefficients` row sums are above the allowed slack.", call. = FALSE)
        }
        if (!all(a >= 1 - slack - slack_tol)) {
            stop("Some `coefficients` row sums are below the allowed slack.", call. = FALSE)
        }
    } else if (!isTRUE(all.equal(rowSums(coefficients), rep(1, K), check.attributes = FALSE))) {
        stop("Coefficients must be row-stochastic (each row sums to 1) when slack = 0")
    }

    if (!is_row_stochastic(compositions)) {
        stop("Compositions must be row-stochastic (each row sums to 1)")
    }

    if (!inherits(loss, "data.frame")) stop("loss must be compatible with data.frame")
    family <- family %||% "gaussian"
    if (!is_non_empty_string(family)) {
        stop("`family` must be a single non-empty string.", call. = FALSE)
    }
    if (!is.list(fit_info)) {
        stop("`fit_info` must be a list.", call. = FALSE)
    }
    if (!is.null(weights)) {
        if (length(weights) != N) {
            stop("`weights` must have one value per sample.", call. = FALSE)
        }
        if (!is_all_non_negative(weights)) {
            stop("`weights` must be finite and non-negative.", call. = FALSE)
        }
    }

    out <- structure(
        list(
            A            = A,
            coefficients = coefficients,
            compositions = compositions,
            slack        = slack,
            init         = init,
            loss         = loss,
            converged    = converged,
            data         = data,
            call         = call,
            family       = family,
            weights      = weights,
            fit_info     = fit_info,
            feature_map  = feature_map
        ),
        class = "archetypes"
    )
    realized <- .aa_realize_coordinates(out)
    if (!is.null(init) && !is.null(realized) && any(dim(init) != dim(realized))) {
        stop("Initial archetypes must have the same number of columns as coordinates")
    }
    out
}


# Methods for Base Archetypes Class ------------------------------------------------

#' Coefficients for archetypes objects
#'
#' Returns the weights that express each archetype as a combination of samples.
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#'
#' @return
#' A numeric matrix (K x N) where each row contains the sample weights used to
#' form the corresponding archetype.
#'
#' @seealso [fitted.archetypes()], [predict.archetypes()], [residuals.archetypes()]
#'
#' @examples
#' fit <- run_aa(iris[, 1:4], K = 3)
#' coefficients(fit)
#'
#' @exportS3Method
coefficients.archetypes <- function(object, ...) {
    object[["coefficients"]]
}


#' Access archetype coordinates and compositions
#'
#' `coordinates()` returns archetype coordinates in the same dimensions as the original input data.
#' `compositions()` returns the weights that express each sample as a combination of archetypes.
#'
#' @param object An archetype analysis result.
#' @param ... Ignored.
#'
#' @return
#' `coordinates()` returns a numeric matrix with K rows of archetype coordinates.
#' `compositions()` returns an N x K numeric matrix of composition weights.
#'
#' @details
#'
#' For kernel fits `coordinates()` return an input-space proxy computed from
#' the fitted generator weights and stored data, since the real coordinates
#' live in the space where the kernel inner product is defined.
#'
#' For non-Gaussian families , coordinates live in parameter space (logit for
#' binomial, log for Poisson, etc.) and may not be directly comparable to the
#' original data.
#'
#' @rdname archetype_accessors
#' @export
coordinates <- function(object, ...) {
    UseMethod("coordinates")
}

#' @rdname archetype_accessors
#' @exportS3Method
coordinates.archetypes <- function(object, ...) {
    A <- .aa_realize_coordinates(object)
    if (is.null(A)) {
        return(NULL)
    }
    if (inherits(object[["data"]], "fd")) {
        return(.aa_fd_from_rows(A, object[["data"]], rownames(A)))
    }
    A
}


#' @rdname archetype_accessors
#' @export
compositions <- function(object, ...) {
    UseMethod("compositions")
}

#' @rdname archetype_accessors
#' @exportS3Method
compositions.archetypes <- function(object, ...) {
    object[["compositions"]]
}


#' Archetype names
#'
#' Get or set the names of archetypes in an archetype analysis result.
#'
#' @param x An archetype analysis result.
#' @param value Character vector with one name per archetype.
#'
#' @return
#' `anames()` returns a character vector. The replacement method returns `x`
#' with names updated consistently across archetype coordinates, coefficients,
#' compositions, and initial coordinates when present.
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' fit <- run_aa(as.matrix(toy), K = 3, max_iter = 20, tol_r2 = 0.95)
#' anames(fit)
#' anames(fit) <- c("A", "B", "C")
#' anames(fit)
#'
#' @rdname archetypes
#' @export
anames <- function(x) {
    UseMethod("anames")
}

#' @rdname archetypes
#' @export
`anames<-` <- function(x, value) {
    UseMethod("anames<-")
}

#' @rdname archetypes
#' @exportS3Method
anames.archetypes <- function(x) {
    rownames(x[["coefficients"]])
}

#' @rdname archetypes
#' @method anames<- archetypes
#' @export
`anames<-.archetypes` <- function(x, value) {
    K <- nrow(x[["coefficients"]])
    if (length(value) != K) {
        fmt <- "Expected %d archetype names, got %d"
        stop(sprintf(fmt, K, length(value)))
    }
    if (any(is.na(value))) {
        stop("Archetype names must not be missing.", call. = FALSE)
    }
    if (!all(nzchar(value))) {
        stop("Archetype names must not be empty.", call. = FALSE)
    }
    if (anyDuplicated(value)) {
        stop("Archetype names must be unique.", call. = FALSE)
    }

    rownames(x[["coefficients"]]) <- value
    colnames(x[["compositions"]]) <- value
    if (!is.null(x[["init"]])) {
        rownames(x[["init"]]) <- value
    }
    if (!is.null(x[["A"]])) {
        rownames(x[["A"]]) <- value
    }

    x
}

#' Fitted values for archetypes objects
#'
#' Computes the projection of the data on the convex hull
#' of the archetypes fit as `compositions %*% coordinates`.
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#'
#' @return
#' Fitted values in the input feature space for Gaussian fits and in the family parameter space for non-Gaussian fits.
#'
#' @seealso [residuals.archetypes()], [predict.archetypes()], [coefficients.archetypes()]
#'
#' @examples
#' toy <- read.csv(system.file("extdata", "toy.csv", package = "yaap"))
#' fit <- run_aa(as.matrix(toy), K = 3, max_iter = 20, tol_r2 = 0.95)
#' Xhat <- fitted(fit)
#' dim(Xhat)
#'
#' @exportS3Method
fitted.archetypes <- function(object, ...) {
    S <- compositions(object)
    A <- .aa_input_coordinates_matrix(object)
    X_hat <- S %*% A
    if (inherits(object[["data"]], "fd")) {
        return(.aa_fd_from_rows(X_hat, object[["data"]], rownames(X_hat)))
    }
    X_hat
}

#' Residuals for archetypes objects
#'
#' @param object An object of class `archetypes`.
#' @param type Residual type. `"response"` (default) returns observed minus fitted values on the observation scale.
#'   `"pearson"` returns GLM-style Pearson residuals, dividing response
#'   residuals by the square root of the family variance function and applying
#'   stored row weights when present.
#' @param ... Ignored.
#'
#' @return
#' Residuals in the same representation as `data`.
#'
#' @seealso
#' [fitted.archetypes()], [predict.archetypes()].
#' For kernel fits see [residuals.kernel_archetypes()].
#'
#' @exportS3Method
residuals.archetypes <- function(object,
                                 type = c("response", "pearson"),
                                 ...) {
    type <- match.arg(type)
    X <- object[["data"]]
    if (is.null(X)) {
        stop("Original data is missing from `object$data`.", call. = FALSE)
    }

    X_hat <- fitted(object)
    if (inherits(X, "fd")) {
        return(X - X_hat)
    }

    stopifnot(all(dim(X) == dim(X_hat)))
    fitted_param <- X_hat
    if (identical(object[["family"]] %||% "gaussian", "multinomial")) {
        X_hat <- rowSums(as.matrix(X)) * X_hat
    }
    out <- X - X_hat
    if (identical(type, "response")) {
        return(out)
    }

    # Pearson residuals
    variance <- .aa_family_variance(object, fitted_param)
    out <- out / sqrt(variance)
    weights <- object[["weights"]] %||% 1
    sqrt(weights) * out
}

.aa_family_variance <- function(object, fitted) {
    family <- object[["family"]] %||% "gaussian"
    X <- object[["data"]]
    eps <- sqrt(.Machine$double.eps)

    variance <- switch(family,
        gaussian    = 1,
        poisson     = fitted,
        binomial    = fitted * (1 - fitted),
        multinomial = rowSums(as.matrix(X)) * fitted * (1 - fitted)
    )

    if (is.null(variance)) {
        stop(
            sprintf("Pearson residuals are not defined for family '%s'.", family),
            call. = FALSE
        )
    }

    pmax(variance, eps)
}

#' Predict compositions or reconstructions for new data from an archetypes model
#'
#' Projects new samples onto the archetype space by solving for composition.
#' The representation is either in the original coordinates (`type = "reconstruction"`)
#' or as barycentric coordinates (`type = "compositions"`) on the archetype simplex.
#'
#' @param object An object of class `archetypes`.
#' @param newdata New data to fit. Must contain the features (columns) used to
#'   fit `object`; fd-backed fits may pass an `fda::fd` object.
#' @param type Prediction output type.
#'   `"reconstruction"` (default) returns `X_hat = S %*% A`;
#'   `"compositions"` returns `S`.
#' @param ... Passed to \code{\link{fit_simplex}} and family-specific
#'   prediction routines.
#'
#' @return
#' If `type = "compositions"`, a numeric matrix (N_new x K) with non-negative
#' row-stochastic composition weights. If `type = "reconstruction"`, a
#' reconstruction matrix (N_new x M) in the input feature space.
#'
#' @details
#' For non-Gaussian families, the reconstructions live in parameter space.
#' In particular, "binomial" family returns probabilities and "poisson"
#' family returns positive rates. For multinomial family, the reconstructions
#' are category probabilities, so each row sums to one.
#'
#' @seealso
#' [fit_simplex()], [residuals.archetypes()], [fitted.archetypes()].
#' For class-specific overrides see [predict.directional_archetypes()] and
#' [predict.kernel_archetypes()].
#'
#' @exportS3Method
predict.archetypes <- function(object,
                               newdata,
                               type = c("reconstruction", "compositions"),
                               ...) {
    type <- match.arg(type)
    is_fd_newdata <- inherits(newdata, "fd")
    if (is_fd_newdata) {
        newdata <- .aa_fd_to_matrix(newdata)
    }

    # Compositions Branch
    family <- object[["family"]] %||% "gaussian"
    if (identical(family, "gaussian")) {
        info <- object[["fit_info"]]
        # TODO: consider supporting prediction for missing-data fits
        if (isTRUE(info[["missing"]])) {
            stop("predict() is not supported for missing-data Gaussian fits.", call. = FALSE)
        }

        map <- object[["feature_map"]]
        X <- .aa_feature_map_transform(map, newdata)
        A <- object[["A"]]
        fit_args <- list(...)
        fit_args[["eps"]] <- fit_args[["eps"]] %||% map[["eps"]] %||% 0
        S <- do.call(fit_simplex, c(list(A = A, X = X), fit_args))
    } else {
        A <- .aa_input_coordinates_matrix(object)
        X <- if (!is.null(colnames(A)) && !is.null(colnames(newdata))) {
            if (inherits(newdata, "data.table")) {
                newdata[, colnames(A), with = FALSE]
            } else {
                newdata[, colnames(A), drop = FALSE]
            }
        } else {
            newdata
        }
        S <- .aa_paa_predict_S(object, X, ...)
    }

    if (identical(type, "compositions")) {
        return(S)
    }

    # Reconstruction Branch
    A <- .aa_input_coordinates_matrix(object)
    X_hat <- S %*% A
    if (is_fd_newdata && inherits(object[["data"]], "fd")) {
        return(.aa_fd_from_rows(X_hat, object[["data"]], rownames(X_hat)))
    }
    X_hat
}

#' @exportS3Method
print.archetypes <- function(x, ...) {
    call_str <- paste(deparse(x[["call"]]), sep = "\n", collapse = "\n")
    cat("\nCall:\n", call_str, "\n\n", sep = "")
    cat(.aa_print_title(x), ":\n", sep = "")
    K <- nrow(x[["coefficients"]])
    loss <- x[["loss"]]
    cat("Number of Archetypes:", K, "\n")
    conv_info <- sprintf(
        "%s after %d iterations.\n",
        ifelse(x[["converged"]], "Converged", "DID NOT CONVERGE"),
        nrow(loss) - 1L
    )
    cat(conv_info)
    cat("Final Loss Metrics:\n")
    print(loss[nrow(loss), , drop = FALSE], row.names = FALSE)
    cat("\n")
    invisible(x)
}

#' @exportS3Method
summary.archetypes <- function(object, ...) {
    loss <- object[["loss"]]
    n_archetypes <- nrow(coefficients(object))
    n_features <- NA_integer_
    A <- .aa_input_coordinates_matrix(object)
    if (!is.null(A)) {
        n_archetypes <- nrow(A)
        n_features <- ncol(A)
    }
    out <- list(
        call         = object[["call"]],
        title        = .aa_print_title(object),
        fit_info     = object[["fit_info"]],
        n_archetypes = n_archetypes,
        n_samples    = nrow(object[["compositions"]]),
        n_features   = n_features,
        converged    = object[["converged"]],
        n_iter       = nrow(loss) - 1L,
        loss         = loss,
        final_loss   = loss[nrow(loss), , drop = FALSE],
        coordinates  = coordinates(object)
    )
    class(out) <- "summary.archetypes"
    out
}

#' @exportS3Method
print.summary.archetypes <- function(x, ...) {
    cat(x[["title"]], " Summary:\n", sep = "")
    cat("Number of Archetypes:", x[["n_archetypes"]], "\n")
    cat("Number of Samples:", x[["n_samples"]], "\n")
    if (!is.na(x[["n_features"]])) {
        cat("Number of Features:", x[["n_features"]], "\n")
    }

    info <- x[["fit_info"]]
    if (is.list(info) && length(info) > 0L) {
        cat("\nFit Details:\n")
        print(as.data.frame(info, stringsAsFactors = FALSE), row.names = FALSE)
    }

    cat("\nConvergence:\n")
    cat(ifelse(x[["converged"]], "Converged", "DID NOT CONVERGE"),
        " after ", x[["n_iter"]], " iterations.\n",
        sep = ""
    )

    cat("\nFinal Loss Metrics:\n")
    print(x[["final_loss"]], row.names = FALSE)

    if (!is.null(x[["coordinates"]])) {
        cat("\nCoordinates:\n")
        print(x[["coordinates"]])
    }
    invisible(x)
}

#' Plot method for archetypes objects
#'
#' Draws diagnostic and geometric plots for a fitted archetypal analysis model.
#'
#' @param x An object of class `archetypes`
#' @param what Character string naming the plot to draw. Supported values are:
#'   \describe{
#'   \item{`"composition"`, `"composision"`}{Draw a stacked barplot of
#'   fitted samples compositions, with bars representing observations and colors
#'   representing archetypes. Rows and columns are clustered by default.}
#'   \item{`"ternary"`, `"simplex"`}{Compositions are plotted as points in the
#'   2D simplex, with corners representing archetypes. For K > 3, multiple 2D
#'   projections are drawn. This plot requires package `compositions`.}
#'   \item{`"loss"`}{Plot the evolution of the objective value across optimization iterations.}
#'   \item{`"profiles"`}{Plot a barplot of the fitted archetypes.
#'   Each set of bars represents a dimension and archetypes are separated by color.}
#'   \item{`"coordinates"`}{Scatterplot of archetype coordinates, optionally over the
#'   original observations when data are available.}
#'   }
#' @param subset Optional sample subset for plots that display observations.
#'   May be numeric row indices, sample names, or a logical vector.
#'   Subsetting is applied before clustering or projecting.
#' @param plot Logical. Should the plot be drawn? If `FALSE`, the prepared
#'   plotting data is returned without drawing.
#' @param ... Additional graphical parameters passed to the selected plotting helper.
#'
#' @return The prepared data used by the selected plot, returned invisibly
#'   when `plot = TRUE`.
#'
#' @importFrom graphics lines pairs plot points
#' @exportS3Method
plot.archetypes <- function(x,
                            what = c("compositions", "loss", "coordinates", "profiles"),
                            subset = NULL,
                            plot = TRUE,
                            ...) {
    stopifnot(inherits(x, "archetypes"))
    what <- .aa_plot_normalize_what(what)
    dots <- list(...)

    if (what == "compositions") {
        S <- .aa_plot_subset_rows(as.matrix(compositions(x)), subset)
        args <- list(plot = plot, cluster_rows = TRUE, cluster_cols = TRUE) %|p|% dots
        args[["compositions"]] <- S
        return(do.call(plot_archetypes_compositions, args))
    }

    if (what %in% c("ternary", "simplex")) {
        S <- .aa_plot_subset_rows(as.matrix(compositions(x)), subset)
        colnames(S) <- colnames(S) %||% paste0("A", seq_len(ncol(S)))
        if (isTRUE(plot)) {
            if (!requireNamespace("compositions", quietly = TRUE)) {
                stop(
                    "Ternary/simplex plots require package `compositions`. ",
                    "Install it to use `plot(..., what = 'ternary'/'simplex')`.",
                    call. = FALSE
                )
            }
            args <- list(x = compositions::acomp(S), axes = TRUE) %|p|% dots
            do.call(graphics::plot, args)
            # main is ignored by compositions::plot.acomp() so add it separately if provided
            graphics::title(main = dots[["main"]])
        }
        # target package ggtern uses wide format no "pca"
        out <- as.data.frame(S)
        if (isTRUE(plot)) {
            return(invisible(out))
        }
        return(out)
    }

    if (what == "loss") {
        args <- list(plot = plot) %|p|% dots
        args[["loss"]] <- x[["loss"]]
        return(do.call(plot_archetypes_loss, args))
    }

    if (what == "profiles") {
        coordinates <- coordinates(x)
        args <- list(
            family = x[["family"]] %||% "gaussian",
            archetype_names = anames(x),
            plot = plot
        ) %|p|% dots
        args[["coordinates"]] <- coordinates
        return(do.call(plot_archetypes_profiles, args))
    }

    family <- x[["family"]] %||% "gaussian"
    if (!(family %in% c("gaussian", "watson"))) {
        msg <- paste(
            "Coordinate plots are not defined when archetype coordinates live",
            "in parameter space but `data` lives in observation space;",
            sprintf("family = '%s'", family)
        )
        stop(msg, call. = FALSE)
    }

    args <- list(data = x[["data"]], archetype_names = anames(x), plot = plot) %|p|% dots
    if (!is.null(args[["data"]])) {
        X <- if (inherits(args[["data"]], "fd")) {
            .aa_fd_to_matrix(args[["data"]])
        } else {
            as.matrix(args[["data"]])
        }
        args[["data"]] <- .aa_plot_subset_rows(X, subset)
    }
    coord <- coordinates(x)
    args[["coordinates"]] <- if (inherits(coord, "fd")) .aa_fd_to_matrix(coord) else coord
    do.call(plot_archetypes_coordinates, args)
}

#' AIC for archetypes objects
#'
#' Computes the AIC-like validity criterion for an `archetypes` object.
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#' @return Numeric scalar with the adapted AIC-like criterion.
#'
#' @details
#' This is not the classical likelihood-based Akaike Information Criterion
#' \eqn{-2 \log L + 2k}. Instead, it implements the adapted archetypal-analysis
#' criterion proposed by Suleman (2017), using the reconstruction error
#' and an efficiency-adjusted complexity penalty.
#' The value is intended for comparing Euclidean Gaussian archetype fits on
#' the same data, especially across different numbers of archetypes `K`.
#' For `K = 1`, the efficiency adjustment is undefined and `AIC()` returns
#' `NA_real_`.
#'
#' @references
#' A. Suleman, "Validation of archetypal analysis"
#' 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE),
#' Naples, Italy, 2017, pp. 1-6, \doi{10.1109/FUZZ-IEEE.2017.8015385}
#'
#' @exportS3Method
AIC.archetypes <- function(object, ...) {
    if (any(object[["slack"]] > 0)) {
        warning(paste(
            "AIC computation assumes coefficients are row-stochastic;",
            "slack > 0 may violate this assumption."
        ))
    }

    X <- object[["data"]]
    if (is.null(X)) {
        stop(paste(
            "AIC requires original data `X`;",
            "provide it when constructing the archetypes object."
        ))
    }
    X <- if (inherits(X, "fd")) .aa_fd_to_matrix(X) else X

    # Compute AIC
    S <- compositions(object)
    K <- ncol(S) # Number of Archetypes
    N <- nrow(X)
    M <- ncol(X)
    nelem <- prod(dim(X)) # number of elements in X
    family <- object[["family"]] %||% "gaussian"
    if (!identical(family, "gaussian")) {
        stop("AIC is not defined for non-Gaussian archetypes objects.", call. = FALSE)
    }
    if (N <= M) {
        warning(
            paste(
                "Adapted AIC is undefined when the number of samples",
                "is not larger than the number of features; returning NA."
            ),
            call. = FALSE
        )
        return(NA_real_)
    }

    # K = 1 is the mean-model baseline. Its reconstruction covariance is zero,
    # so the efficiency-adjusted adapted AIC is undefined.
    if (K == 1L) {
        return(NA_real_)
    }

    A <- .aa_input_coordinates_matrix(object)
    X_hat <- S %*% A
    npar <- N * (K - 1) + K * (N - 1) + 1 # K_mu + K_beta + 1

    eta <- tryCatch(.aa_effic(X, X_hat), error = function(e) NA_real_)
    if (!is.finite(eta) || eta <= 0) {
        warning(
            paste(
                "Adapted AIC is undefined because the efficiency term",
                "is non-positive or non-finite; returning NA."
            ),
            call. = FALSE
        )
        return(NA_real_)
    }
    rss <- norm(X - X_hat, "F")^2 # "loss" is not necessarily the RSS on the original data
    log(rss / nelem) + 2 * npar / (N * eta)
}


# Helper functions ------------------------------------------------------------

.aa_canonical_archetype_order <- function(coordinates, coefficients) {
    key <- if (!is.null(coordinates)) {
        cbind(as.matrix(coordinates), as.matrix(coefficients))
    } else {
        as.matrix(coefficients)
    }
    K <- nrow(key)
    if (K <= 1L) {
        return(seq_len(K))
    }

    key <- as.data.frame(key, optional = TRUE)
    key[[".index"]] <- seq_len(K)
    do.call(order, c(key, list(na.last = TRUE)))
}

.aa_realize_coordinates <- function(object) {
    .aa_feature_map_inverse(
        object[["feature_map"]],
        object[["A"]],
        object = object
    )
}

.aa_input_coordinates_matrix <- function(object) {
    A <- coordinates(object)
    if (is.null(A)) {
        return(NULL)
    }
    if (inherits(A, "fd")) {
        return(.aa_fd_to_matrix(A))
    }
    as.matrix(A)
}

.aa_fd_to_matrix <- function(x) {
    if (!requireNamespace("fda", quietly = TRUE)) {
        stop("Package `fda` is required for fd data.", call. = FALSE)
    }
    t(stats::coef(x))
}

.aa_fd_from_rows <- function(X, template, curve_names = NULL) {
    if (!requireNamespace("fda", quietly = TRUE)) {
        stop("Package `fda` is required for fd data.", call. = FALSE)
    }

    fdnames <- template[["fdnames"]]
    curve_names <- curve_names %||% rownames(X)
    if (!is.null(curve_names)) {
        fdnames[["reps"]] <- curve_names
    }

    fda::fd(t(X), template[["basis"]], fdnames = fdnames)
}

# Normalises the `what` argument for archetype plot methods.
.aa_plot_normalize_what <- function(what) {
    what <- match.arg(
        tolower(what[1L]),
        c(
            "compositions", "composition", "composision", "composisions",
            "ternary", "simplex",
            "loss", "coordinates", "profiles"
        )
    )
    if (what %in% c("composition", "composision", "composisions")) {
        what <- "compositions"
    }
    what
}

# Subsets rows of a matrix by index, name, or logical vector.
.aa_plot_subset_rows <- function(z, subset) {
    if (is.null(subset)) {
        return(z)
    }
    if (is.logical(subset) && length(subset) != nrow(z)) {
        stop("Logical `subset` must have one value per sample", call. = FALSE)
    }
    if (is.character(subset) && anyNA(match(subset, rownames(z)))) {
        stop("Some `subset` values are not sample names", call. = FALSE)
    }
    out <- z[subset, , drop = FALSE]
    if (nrow(out) == 0L) {
        stop("`subset` selects no samples", call. = FALSE)
    }
    out
}

.aa_print_title <- function(x) {
    info <- x[["fit_info"]]
    if (!is.list(info) || length(info) == 0L) {
        return("Archetypal Analysis")
    }

    secondary <- if (isTRUE(info[["robust"]])) "Robust" else character()
    method <- switch(info[["method"]] %||% "",
        pgd = "PGD",
        nnls = "NNLS",
        paa = "PAA",
        kernel = "Kernel",
        directional = "Directional",
        info[["method"]] %||% ""
    )
    family <- info[["family"]]
    suffix <- if (!is.null(info[["kernel"]])) {
        sprintf(" (%s kernel)", info[["kernel"]])
    } else if (identical(info[["method"]], "paa") && !is.null(family)) {
        sprintf(" (%s family)", family)
    } else {
        ""
    }

    paste(c(secondary, method, "Archetypal Analysis"), collapse = " ") |>
        paste0(suffix)
}
