#' Archetype Analysis Result Object
#'
#' Creates an S3 object of class `archetypes` holding the fitted archetypes,
#' their decomposition of the data, and diagnostic information from the fitting
#' process.
#'
#' @param coordinates Numeric matrix (K x M) giving the archetype coordinates
#'   in the original feature space. `K` is the number of archetypes, `M` is the
#'   number of features.
#' @param coefficients Numeric matrix (K x N) giving the weights that express
#'   each archetype as a linear combination of the N samples in the input data.
#' @param compositions Numeric matrix (N x K) giving the barycentric
#'   coordinates of each sample in the archetype space (how much each archetype
#'   contributes to each sample).
#' @param loss Data frame containing per-iteration metrics
#' @param converged Logical. Whether the optimization converged.
#' @param data Numeric matrix (N x M) with the original data used for
#'   fitting. If supplied, residuals are computed against
#'   `compositions %*% coordinates`.
#' @param call The matched function call that created the object (defaults to `match.call()`).
#' @param init Optional numeric matrix (K x M) with the initial archetype
#'   coordinates before optimization.
#'
#' @details
#' Let `X` be the original data matrix (N x M).
#'
#' The fitted values are:
#'
#' ```
#'   X_hat = compositions %*% coordinates
#'   coordinates = coefficients %*% X
#' ```
#' `compositions` is always a row-stochastic matrix (each row sums to 1).
#' `coefficients` is typically a row-stochastic matrix as well unless the
#' constraint is relaxed (`slack` > 0) in which case the archetypes can lay
#' outside the convex hull of the data.
#'
#' ## Invariants and Dimensions
#'
#' * `nrow(coordinates) == ncol(compositions) == nrow(coefficients) == K`
#' * `nrow(compositions) == ncol(coefficients) == N`
#' * `rowSums(compositions) == 1`
#' * If `data` is supplied: `dim(data) == c(N, M)` where `M = ncol(coordinates)`
#'
#' @return
#' An object of class `archetypes`, which is a list with components:
#'
#' - `coordinates` - (K x M) archetype coordinates.
#' - `coefficients` - (K x N) archetype weights on samples.
#' - `compositions` - (N x K) sample weights on archetypes.
#' - `init` - optional initial coordinates.
#' - `loss` - data frame of per-iteration metrics.
#' - `converged` - logical convergence flag.
#' - `residuals` - `data - compositions %*% coordinates`
#' - `call` - the matched call.
#'
#' @seealso
#' [plot.archetypes()], [predict.archetypes()], [residuals.archetypes()], [fitted.archetypes()]
#'
#' @export
archetypes <- function(coordinates,
                       coefficients,
                       compositions,
                       slack = 0,
                       loss = NULL,
                       converged = TRUE,
                       call = NULL,
                       data = NULL,
                       init = NULL) {
    if (is.null(call))
        call <- match.call()
    if (is.null(loss))
        loss <- data.frame(rss = NA_real_, r2  = NA_real_, k_S = NA_real_, k_A = NA_real_)

    return(
        new_archetypes(
            coordinates  = coordinates,
            coefficients = coefficients,
            compositions = compositions,
            slack        = slack,
            loss         = loss,
            AIC          = NA_real_,
            converged    = converged,
            call         = call,
            data         = data,
            init         = init
        )
    )
}

#' Archetypes object constructor
new_archetypes <- function(coordinates,
                           coefficients,
                           compositions,
                           slack,
                           loss,
                           converged,
                           AIC = NA_real_,
                           call = NULL,
                           data = NULL,
                           init = NULL) {


    # Dimension checks
    K <- nrow(coordinates)   # number of archetypes
    # M <- ncol(coordinates)   # number of features
    N <- nrow(compositions)  # number of samples


    # Check Coordinates
    if (!is.null(init) && any(dim(init) != dim(coordinates)))
        stop("Initial archetypes must have the same number of columns as coordinates")

    # Check Coefficients
    if (nrow(coefficients) != K) {
        fmt <- "nrow(coefficients) = %d does not match number of archetypes (%d)"
        stop(sprintf(fmt, nrow(coefficients), K))
    }

    if (ncol(coefficients) != N) {
        fmt <- "Inconsistent number of samples between compositions (%d) and coefficients (%d)"
        stop(sprintf(fmt, N, ncol(coefficients)))
    }

    stopifnot("All slack values must be non-negative" = all(slack >= 0))
    if (any(slack > 0)) {
        a <- rowSums(coefficients)
        ix <- which(a > 1 + slack)
        if (length(ix) > 0) {
            fmt <- "Some rowSums(coefficients) are above allowed slack: %s"
            if (length(ix) > 10) {
                ix <- head(ix, 10)
                fmt <- paste(fmt, "... (truncated)")
            }
            stop(sprintf(fmt, paste(ix, collapse = ", ")))
        }

        ix <- which(a < 1 - slack)
        if (length(ix) > 0) {
            fmt <- "Some rowSums(coefficients) are below allowed slack: %s"
            if (length(ix) > 10) {
                ix <- head(ix, 10)
                fmt <- paste(fmt, "... (truncated)")
            }
            stop(sprintf(fmt, paste(ix, collapse = ", ")))
        }
    } else if (!all.equal(rowSums(coefficients), rep(1, K), check.attributes = FALSE)) {
        stop("Coefficients must be row-stochastic (each row sums to 1) when slack = 0")
    }

    # Check Compositions
    if (ncol(compositions) != K) {
        fmt <- "ncol(compositions) = %d does not match number of archetypes (%d)"
        stop(sprintf(fmt, ncol(compositions), K))
    }
    if (!all.equal(rowSums(compositions), rep(1, N), check.attributes = FALSE))
        stop("Compositions must be row-stochastic (each row sums to 1)")

    # Check loss
    if(!inherits(loss, "data.frame")) stop("loss must be compatible with data.frame")

    structure(
        list(
            coordinates  = coordinates,
            coefficients = coefficients,
            compositions = compositions,
            slack        = slack,
            init         = init,
            loss         = loss,
            AIC          = AIC,
            converged    = converged,
            call         = call
        ),
        class = "archetypes"
    )
}

#' Coefficients for archetypes objects
#'
#' Returns the `coefficients` component of an `archetypes` object,
#' i.e., the weights that express each archetype as a convex (or
#' linear, if constraint is relaxed) combination of samples.
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#'
#' @return
#' A numeric matrix (K x N) where each row contains the weights
#' over samples used to form one archetype.
#'
#' @seealso [archetypes()], [fitted.archetypes()],
#'   [predict.archetypes()], [residuals.archetypes()]
#'
#' @examples
#' # coefficients(fit)
#'
#' @exportS3Method
coefficients.archetypes <- function(object, ...)
    return(object[["coefficients"]])

#' Fitted values for archetypes objects
#'
#' Computes the projection of the data on the convex hull
#' of the archetypes fit as `compositions %*% coordinates`.
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#'
#' @return
#' A numeric matrix (N x M) of fitted values.
#'
#' @seealso [archetypes()], [residuals.archetypes()],
#'   [predict.archetypes()], [coefficients.archetypes()]
#'
#' @examples
#' # Xhat <- fitted(fit)
#'
#' @exportS3Method
fitted.archetypes <- function(object, ...)
    return(with(object, compositions %*% coordinates))

#' Residuals for archetypes objects
#'
#' @param object An object of class `archetypes`.
#' @param data Numeric matrix (N x M) with the original data. If not provided, the
#'   function attempts to use `object$data`.
#' @param ... Ignored.
#'
#' @return
#' A numeric matrix of residuals.
#'
#' @exportS3Method
residuals.archetypes <- function(object, data = NULL, ...) {
    X <- if (is.null(data)) {
        if (is.null(object[["data"]])) {
            msg <- paste("Original data must be provided either when",
                         "constructing the archetypes object or as an argument",
                         "to `residuals.archetypes()`")
            stop(msg)
        }
        object[["data"]]
    } else {
        data
    }
    X_hat <- fitted(object)
    stopifnot(all(dim(X) == dim(X_hat)))
    return(X - X_hat)
}

#' Predict compositions for new data from an archetypes model
#'
#' Projects new samples onto the archetype space by solving
#' non-negative least squares (NNLS) for the composition weights
#' given fixed archetype coordinates.
#'
#' @param object An object of class `archetypes`.
#' @param newdata New data to fit. Must contain the features (columns) used to fit `object`.
#' @param ... Passed to \code{\link{fit_simplex}}.
#'
#' @return
#' A numeric matrix (N_new x K) of non-negative composition weights for `newdata`.
#'
#' @seealso [archetypes()], [fit_simplex()]
#'
#' @exportS3Method
predict.archetypes <- function(object, newdata, ...) {
    A <- object[["coordinates"]]
    X <- if (!is.null(colnames(A))) {  # extract relevant columns
        if (class(newdata) == "data.table") {
            newdata[, colnames(A), with = FALSE]
        } else {
            newdata[, colnames(A), drop = FALSE]
        }
    } else {
        newdata
    }
    S <- fit_simplex(A, X, ...)
    return(S)
}

#' @exportS3Method
print.archetypes <- function(object, ...) {
    call_str <- paste(deparse(object[["call"]]), sep = "\n", collapse = "\n")
    cat("\nCall:\n", call_str, "\n\n", sep = "")
    cat("Archetypes Summary:\n")
    K <- nrow(object[["coordinates"]])
    cat("Number of Archetypes:", K, "\n")
    conv_info <- sprintf(
        "%s after %d iterations.\n",
        ifelse(object[["converged"]], "Converged", "DID NOT CONVERGE"),
        nrow(object[["loss"]]) - 1L
    )
    cat(conv_info)
    cat("Final Loss Metrics:\n")
    print(tail(object[["loss"]], 1), row.names = FALSE)
    cat("\n")
    invisible(object)
}

#' Plot method for archetypes objects
#'
#' @param x An object of class `archetypes`
#'
#' @exportS3Method
plot.archetypes <- function(x, what = c("compositions", "residuals", "coefficients"), ...) {
    stopifnot(inherits(x, "archetypes"))
    stop("Not yet implemented")
}

#' AIC for archetypes objects
#'
#' Computes the Akaike Information Criterion (AIC) for an `archetypes` object
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#' @return Numeric scalar with the AIC value.
#'
#' @references
#' A. Suleman, "Validation of archetypal analysis,"
#' 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE),
#' Naples, Italy, 2017, pp. 1-6, \url{www.doi.org/10.1109/FUZZ-IEEE.2017.8015385}
#'
#' @exportS3Method
AIC.archetypes <- function(object, ...) {
    if (any(object[["slack"]] > 0))
        warning(paste("AIC computation assumes coefficients are row-stochastic;",
                      "slack > 0 may violate this assumption."))

    aic <- object[["AIC"]]
    if (!is.na(aic)) return(aic)  # check for cached value
    X <- object[["data"]]
    if (!is.null(X))
        stop(paste("AIC was not precomputed because original data `X`",
                   "was not provided when constructing the archetypes object."))

    # Compute AIC
    K     <- nrow(object[["coordinates"]])        # Number of Archetypes
    nelem <- prod(dim(X))                         # number of elements in X
    rss   <- tail(object[["loss"]][["rss"]], 1L)  # final RSS
    X_hat <- fitted(object)
    aic   <- log(rss / nelem) + 2 * (2*K - 1) / effic(X, X_hat)

    object[["AIC"]] <- aic  # cache for future calls
    return(aic)
}

