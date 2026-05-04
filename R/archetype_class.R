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
#' @param slack Non-negative numeric scalar or vector giving the allowed
#'   relaxation of coefficient row sums away from 1.
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
#' - `data` - optional original data matrix.
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
#'
#' @param coordinates Numeric matrix (K x M) giving archetype coordinates.
#' @param coefficients Numeric matrix (K x N) giving archetype weights on samples.
#' @param compositions Numeric matrix (N x K) giving sample weights on archetypes.
#' @param slack Non-negative numeric scalar or vector giving allowed coefficient
#'   row-sum relaxation.
#' @param loss Data frame containing per-iteration metrics.
#' @param converged Logical. Whether the optimization converged.
#' @param AIC Numeric scalar with a precomputed AIC value.
#' @param call The matched function call that created the object.
#' @param data Optional original data matrix.
#' @param init Optional initial archetype coordinates.
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
                ix <- utils::head(ix, 10)
                fmt <- paste(fmt, "... (truncated)")
            }
            stop(sprintf(fmt, paste(ix, collapse = ", ")))
        }

        ix <- which(a < 1 - slack)
        if (length(ix) > 0) {
            fmt <- "Some rowSums(coefficients) are below allowed slack: %s"
            if (length(ix) > 10) {
                ix <- utils::head(ix, 10)
                fmt <- paste(fmt, "... (truncated)")
            }
            stop(sprintf(fmt, paste(ix, collapse = ", ")))
        }
    } else if (!isTRUE(all.equal(rowSums(coefficients), rep(1, K), check.attributes = FALSE))) {
        stop("Coefficients must be row-stochastic (each row sums to 1) when slack = 0")
    }

    # Check Compositions
    if (ncol(compositions) != K) {
        fmt <- "ncol(compositions) = %d does not match number of archetypes (%d)"
        stop(sprintf(fmt, ncol(compositions), K))
    }
    if (!isTRUE(all.equal(rowSums(compositions), rep(1, N), check.attributes = FALSE)))
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
            data         = data,
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

#' Archetype names
#'
#' Get or set the names of archetypes in an `archetypes` object.
#'
#' @param x An object of class `archetypes`.
#' @param value Character vector with one name per archetype.
#'
#' @return
#' `names.archetypes()` returns a character vector. The replacement method
#' returns `x` with names updated consistently across archetype coordinates,
#' coefficients, compositions, and initial coordinates when present.
#'
#' @examples
#' # names(fit)
#' # names(fit) <- c("A", "B", "C")
#'
#' @exportS3Method
names.archetypes <- function(x)
    rownames(x[["coordinates"]])

#' @rdname names.archetypes
#' @method names<- archetypes
#' @export
`names<-.archetypes` <- function(x, value) {
    K <- nrow(x[["coordinates"]])
    if (length(value) != K) {
        fmt <- "Expected %d archetype names, got %d"
        stop(sprintf(fmt, K, length(value)))
    }
    stopifnot("Archetype names must not be missing" = !any(is.na(value)))
    stopifnot("Archetype names must not be empty" = all(nzchar(value)))
    stopifnot("Archetype names must be unique" = !anyDuplicated(value))

    rownames(x[["coordinates"]]) <- value
    rownames(x[["coefficients"]]) <- value
    colnames(x[["compositions"]]) <- value
    if (!is.null(x[["init"]]))
        rownames(x[["init"]]) <- value

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
        if (inherits(newdata, "data.table")) {
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
print.archetypes <- function(x, ...) {
    call_str <- paste(deparse(x[["call"]]), sep = "\n", collapse = "\n")
    cat("\nCall:\n", call_str, "\n\n", sep = "")
    cat("Archetypes Summary:\n")
    K <- nrow(x[["coordinates"]])
    loss <- x[["loss"]]
    cat("Number of Archetypes:", K, "\n")
    conv_info <- sprintf(
        "%s after %d iterations.\n",
        ifelse(x[["converged"]], "Converged", "DID NOT CONVERGE"),
        nrow(loss) - 1L
    )
    cat(conv_info)
    cat("Final Loss Metrics:\n")
    print(loss[nrow(loss), ], row.names = FALSE)
    cat("\n")
    invisible(x)
}

#' Plot method for archetypes objects
#'
#' Draws diagnostic and geometric plots for a fitted archetypal analysis model.
#' The `what` argument controls which part of the fit is visualized:
#'
#' * `"compositions"` plots the rows of `x$compositions`, i.e. the weights that
#'   express each observation as a mixture of the archetypes. For three
#'   archetypes this is a ternary plot: points near a corner are dominated by
#'   one archetype, points near an edge mix two archetypes, and points near the
#'   center mix all three.
#' * `"loss"` plots the residual sum of squares stored in `x$loss$rss` across
#'   optimization iterations. This is useful for checking whether the fitting
#'   algorithm reduced reconstruction error and whether the loss curve has
#'   plateaued.
#' * `"coordinates"` plots the original observations together with the fitted
#'   archetype coordinates in feature space. In two dimensions the archetypes
#'   are connected as a closed polygon. In more than two dimensions the method
#'   draws pairwise scatterplots and connects archetypes with dashed closed
#'   paths in each panel. With `projection = "pca"`, higher-dimensional data and
#'   archetypes are first projected to the first two principal components.
#'
#' @param x An object of class `archetypes`
#' @param what Character string naming the plot to draw. Supported values are
#'   `"compositions"`, `"loss"`, and `"coordinates"`. `"composition"` and
#'   `"composision"` are accepted as aliases for `"compositions"`.
#' @param data Optional numeric matrix with the original data. Required for
#'   `what = "coordinates"` if the object does not store its original data.
#' @param projection Projection to use for coordinate plots. Use `"pca"` to
#'   project data with more than two dimensions to the first two principal
#'   component scores before plotting.
#' @param ... Additional graphical parameters passed to the underlying plotting
#'   functions. For `what = "coordinates"`, `col`, `pch`, and `cex` may be
#'   length-two vectors: the first value is used for observations and the second
#'   value is used for archetypes. If only one color is supplied, observations
#'   use that color and archetypes are drawn in red.
#'
#' @importFrom compositions acomp
#' @importFrom graphics lines pairs plot points
#' @exportS3Method
plot.archetypes <- function(x,
                            what = c("compositions", "loss", "coordinates"),
                            data = NULL,
                            projection = c("none", "pca"),
                            ...) {
    stopifnot(inherits(x, "archetypes"))

    what <- match.arg(
        tolower(what[1L]),
        c("compositions", "composition", "composision", "composisions", "loss", "coordinates")
    )

    # Composition Ternary Plot -----------------------------------------------

    if (what %in% c("composition", "composision",  "composisions"))
        what <- "compositions"
    projection <- match.arg(projection)

    dots <- list(...)
    plot_args <- function(defaults, dots) {
        defaults[names(dots)] <- dots
        defaults
    }
    arg_or <- function(name, default) {
        value <- dots[[name]]
        if (is.null(value)) default else value
    }

    if (what == "compositions") {
        S <- as.matrix(x[["compositions"]])
        if (is.null(colnames(S)))
            colnames(S) <- paste0("A", seq_len(ncol(S)))
        args <- plot_args(
            list(x = acomp(S), axes = TRUE),
            dots
        )
        do.call(plot, args)
        return(invisible(x))
    }

    # Loss Plot --------------------------------------------------------------

    if (what == "loss") {
        loss <- x[["loss"]]
        if (is.null(loss[["rss"]]))
            stop("`x$loss` must contain an `rss` column for loss plots")
        args <- plot_args(
            list(
                x = seq_len(nrow(loss)) - 1L,
                y = loss[["rss"]],
                type = "l",
                xlab = "Iteration",
                ylab = "RSS"
            ),
            dots
        )
        do.call(plot, args)
        return(invisible(x))
    }

    # Coordinate Plot -------------------------------------------------------

    X <- if (is.null(data)) x[["data"]] else data
    if (is.null(X)) {
        msg <- paste("Original data must be provided either when constructing",
                     "the archetypes object or through `data` for coordinate plots")
        stop(msg)
    }

    X <- as.matrix(X)
    A <- as.matrix(x[["coordinates"]])
    if (ncol(X) != ncol(A)) {
        fmt <- "`data` has %d columns but `x$coordinates` has %d columns"
        stop(sprintf(fmt, ncol(X), ncol(A)))
    }
    if (ncol(X) < 2L)
        stop("Coordinate plots require at least two dimensions")

    if (projection == "pca" && ncol(X) > 2L) {
        combined <- rbind(X, A)
        pc <- stats::prcomp(combined, center = TRUE, scale. = FALSE)
        scores <- pc[["x"]][, seq_len(2L), drop = FALSE]
        rownames(scores) <- rownames(combined)
        x_projected <- x
        x_projected[["data"]] <- scores[seq_len(nrow(X)), , drop = FALSE]
        x_projected[["coordinates"]] <- scores[nrow(X) + seq_len(nrow(A)), , drop = FALSE]
        colnames(x_projected[["data"]]) <- colnames(x_projected[["coordinates"]]) <- c("PC1", "PC2")
        plot(x_projected, what = "coordinates", projection = "none", ...)
        return(invisible(x))
    }

    data_col <- arg_or("col", "black")
    arch_col <- if (length(data_col) > 1L) data_col[2L] else "red"
    data_col <- data_col[1L]
    data_pch <- arg_or("pch", 1)
    arch_pch <- if (length(data_pch) > 1L) data_pch[2L] else 16
    data_pch <- data_pch[1L]
    cex <- arg_or("cex", 1)
    data_cex <- cex[1L]
    arch_cex <- if (length(cex) > 1L) cex[2L] else 1.3
    lwd <- arg_or("lwd", 1)

    dots[c("col", "pch", "cex", "lwd", "lty")] <- NULL

    if (is.null(colnames(X)))
        colnames(X) <- paste0("V", seq_len(ncol(X)))
    if (is.null(colnames(A)))
        colnames(A) <- colnames(X)

    if (ncol(X) == 2L) {
        args <- plot_args(
            list(
                x = X[, 1L],
                y = X[, 2L],
                xlab = colnames(X)[1L],
                ylab = colnames(X)[2L],
                asp = 1,
                col = data_col,
                pch = data_pch,
                cex = data_cex
            ),
            dots
        )
        do.call(plot, args)
        A_closed <- A[c(seq_len(nrow(A)), 1L), , drop = FALSE]
        lines(A_closed[, 1L], A_closed[, 2L], col = arch_col, lwd = lwd, lty = 1)
        points(A[, 1L], A[, 2L], col = arch_col, pch = arch_pch, cex = arch_cex)
        return(invisible(x))
    }

    combined <- rbind(X, A)
    n_data <- nrow(X)
    n_total <- nrow(combined)
    panel <- function(x, y, ...) {
        data_ix <- seq_len(n_data)
        arch_ix <- (n_data + 1L):n_total
        arch_closed <- c(arch_ix, arch_ix[1L])
        points(x[data_ix], y[data_ix],
               col = data_col, pch = data_pch, cex = data_cex)
        lines(x[arch_closed], y[arch_closed], col = arch_col, lwd = lwd, lty = 2)
        points(x[arch_ix], y[arch_ix],
               col = arch_col, pch = arch_pch, cex = arch_cex)
    }
    args <- plot_args(
        list(x = combined, panel = panel, lower.panel = panel, upper.panel = panel),
        dots
    )
    do.call(pairs, args)
    invisible(x)
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
    if (is.null(X))
        stop(paste("AIC was not precomputed because original data `X`",
                   "was not provided when constructing the archetypes object."))

    # Compute AIC
    K     <- nrow(object[["coordinates"]])        # Number of Archetypes
    nelem <- prod(dim(X))                         # number of elements in X
    rss   <- object[["loss"]][["rss"]]
    X_hat <- fitted(object)
    aic   <- log(rss[length(rss)] / nelem) + 2 * (2*K - 1) / effic(X, X_hat)

    object[["AIC"]] <- aic  # cache for future calls
    return(aic)
}
