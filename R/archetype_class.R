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
#' @param data Original data used for fitting. Usually a numeric matrix (N x M),
#'   but class-specific entry points may store the original input object. If
#'   supplied, residuals are computed against fitted values.
#' @param call The matched function call that created the object (defaults to `match.call()`).
#' @param init Optional numeric matrix (K x M) with the initial archetype
#'   coordinates before optimization.
#' @param family Observation family. Defaults to `"gaussian"`.
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
#' ## Archetype names
#'
#' Archetype labels are stored as row names on `coordinates` and `coefficients`,
#' and as column names on `compositions`. Use `anames()` and `anames<-()` to get
#' or set them consistently:
#'
#' ```
#' anames(fit)
#' anames(fit) <- c("A", "B", "C")
#' ```
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
#' - `data` - optional original data.
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
                       init = NULL,
                       family = "gaussian") {
    if (is.null(call))
        call <- match.call()
    if (is.null(loss))
        loss <- data.frame(loss = NA_real_, r2  = NA_real_, k_S = NA_real_, k_A = NA_real_)

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
            init         = init,
            family       = family
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
#' @param data Optional original data.
#' @param init Optional initial archetype coordinates.
#' @param family Observation family.
new_archetypes <- function(coordinates,
                           coefficients,
                           compositions,
                           slack,
                           loss,
                           converged,
                           AIC = NA_real_,
                           call = NULL,
                           data = NULL,
                           init = NULL,
                           family = "gaussian") {


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
    if (is.null(family))
        family <- "gaussian"
    stopifnot("family must be a single non-empty string" =
                  is.character(family) && length(family) == 1L && nzchar(family))

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
            call         = call,
            family       = family
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
#' # anames(fit)
#' # anames(fit) <- c("A", "B", "C")
#'
#' @rdname archetypes
#' @export
anames <- function(x)
    UseMethod("anames")

#' @rdname archetypes
#' @export
`anames<-` <- function(x, value)
    UseMethod("anames<-")

#' @rdname archetypes
#' @exportS3Method
anames.archetypes <- function(x)
    rownames(x[["coordinates"]])

#' @rdname archetypes
#' @method anames<- archetypes
#' @export
`anames<-.archetypes` <- function(x, value) {
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
#' Fitted values. Usually a numeric matrix (N x M); fd-backed fits return an
#' `fda::fd` object.
#'
#' @seealso [archetypes()], [residuals.archetypes()],
#'   [predict.archetypes()], [coefficients.archetypes()]
#'
#' @examples
#' # Xhat <- fitted(fit)
#'
#' @exportS3Method
fitted.archetypes <- function(object, ...) {
    X_hat <- with(object, compositions %*% coordinates)
    family <- object[["family"]]
    if (is.null(family))
        family <- "gaussian"
    if (identical(family, "multinomial") && !is.null(object[["data"]])) {
        totals <- rowSums(as.matrix(object[["data"]]))
        X_hat <- totals * X_hat
    }
    if (inherits(object[["data"]], "fd"))
        return(.aa_fd_from_rows(X_hat, object[["data"]], rownames(X_hat)))
    X_hat
}

#' Residuals for archetypes objects
#'
#' @param object An object of class `archetypes`.
#' @param data Original data. If not provided, the function attempts to use
#'   `object$data`.
#' @param ... Ignored.
#'
#' @return
#' Residuals in the same representation as `data`.
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
    if (inherits(X, "fd"))
        return(X - X_hat)

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
#' @param newdata New data to fit. Must contain the features (columns) used to
#'   fit `object`; fd-backed fits may pass an `fda::fd` object.
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
    if (inherits(newdata, "fd"))
        newdata <- .aa_fd_to_matrix(newdata)
    X <- if (!is.null(colnames(A))) {  # extract relevant columns
        if (inherits(newdata, "data.table")) {
            newdata[, colnames(A), with = FALSE]
        } else {
            newdata[, colnames(A), drop = FALSE]
        }
    } else {
        newdata
    }
    family <- object[["family"]]
    if (is.null(family))
        family <- "gaussian"
    if (!identical(family, "gaussian"))
        return(.aa_paa_predict_S(object, X, ...))
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

#' Stacked barplot of composition weights
#'
#' Draws a horizontal stacked barplot for a matrix-like set of composition
#' weights, with rows interpreted as samples and columns interpreted as
#' archetypes or other compositional parts.
#'
#' This display is intended for small composition matrices, roughly of tens of
#' samples and/or archetypes. For many samples or many archetypes, a heatmap is
#' usually easier to read.
#'
#' @param x Numeric matrix or data frame. Rows are samples and columns are
#'   archetypes. Rows should contain non-negative composition weights.
#' @param cluster_rows,cluster_cols Logical values or `hclust` objects. When
#'   `TRUE`, rows or columns are reordered by hierarchical clustering. Row
#'   clustering is computed on `compositions::cdt(compositions::acomp(x))`.
#' @param distance Optional distance metric used for both row and column
#'   clustering when `distance_rows` or `distance_cols` are not supplied.
#' @param distance_rows,distance_cols Distance metrics used when clustering
#'   rows or columns. Values may be any method accepted by [stats::dist()],
#'   `"correlation"`, a function that returns a `dist` object, or a precomputed
#'   `dist` object. Row distances are computed after the cdt transform.
#' @param linkage Linkage method passed to [stats::hclust()].
#' @param col Optional vector of colors, one per archetype. Defaults to a
#'   qualitative HCL palette.
#' @param legend Logical. Should an archetype legend be drawn?
#' @param border Border color for the stacked bar segments.
#' @param ... Additional graphical parameters passed to [graphics::barplot()].
#'
#' @return Invisibly returns a list with the reordered matrix, row and column
#'   orders, and clustering objects.
#'
#' @export
composition_barplot <- function(x,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                distance = NULL,
                                distance_rows = "euclidean",
                                distance_cols = "euclidean",
                                linkage = "complete",
                                col = NULL,
                                legend = TRUE,
                                border = NA,
                                ...) {
    S <- as.matrix(x)
    if (!is.numeric(S))
        stop("`x` must be a numeric matrix or data frame", call. = FALSE)
    if (!all(is.finite(S)))
        stop("`x` must contain only finite values", call. = FALSE)
    if (any(S < 0))
        stop("`x` must contain non-negative composition weights", call. = FALSE)
    if (nrow(S) == 0L || ncol(S) == 0L)
        stop("`x` must have at least one row and one column", call. = FALSE)
    if (is.null(rownames(S)))
        rownames(S) <- seq_len(nrow(S))
    if (is.null(colnames(S)))
        colnames(S) <- paste0("A", seq_len(ncol(S)))
    if (!is.null(distance)) {
        if (missing(distance_rows))
            distance_rows <- distance
        if (missing(distance_cols))
            distance_cols <- distance
    }

    make_dist <- function(data, distance) {
        if (inherits(distance, "dist"))
            return(distance)
        if (is.function(distance)) {
            d <- distance(data)
            if (inherits(d, "dist"))
                return(d)
            return(stats::as.dist(d))
        }
        if (!is.character(distance) || length(distance) != 1L || !nzchar(distance))
            stop("Clustering distance must be a non-empty string, function, or dist object", call. = FALSE)
        if (distance == "correlation")
            return(stats::as.dist(1 - stats::cor(t(data))))
        stats::dist(data, method = distance)
    }

    make_hclust <- function(value, data, margin, distance) {
        if (inherits(value, "hclust"))
            return(value)
        if (!isTRUE(value))
            return(NULL)
        if (margin == "rows") {
            if (nrow(data) < 2L)
                return(NULL)
            transformed <- as.matrix(compositions::cdt(compositions::acomp(data)))
            return(stats::hclust(make_dist(transformed, distance), method = linkage))
        }
        if (ncol(data) < 2L)
            return(NULL)
        stats::hclust(make_dist(t(data), distance), method = linkage)
    }

    row_hclust <- make_hclust(cluster_rows, S, "rows", distance_rows)
    col_hclust <- make_hclust(cluster_cols, S, "cols", distance_cols)
    row_order <- if (is.null(row_hclust)) seq_len(nrow(S)) else row_hclust[["order"]]
    col_order <- if (is.null(col_hclust)) seq_len(ncol(S)) else col_hclust[["order"]]
    S_plot <- S[row_order, col_order, drop = FALSE]

    if (is.null(col)) {
        col <- grDevices::hcl.colors(ncol(S_plot), palette = "Dark 3")
    } else if (!is.null(names(col)) && all(colnames(S_plot) %in% names(col))) {
        col <- col[colnames(S_plot)]
    }
    col <- rep_len(col, ncol(S_plot))

    dots <- list(...)
    args <- list(
        height = t(S_plot),
        col = col,
        border = border,
        xlab = "Composition",
        names.arg = rownames(S_plot),
        horiz = TRUE,
        las = 1L
    )
    args[names(dots)] <- dots
    do.call(graphics::barplot, args)
    if (isTRUE(legend)) {
        graphics::legend(
            "topright",
            legend = colnames(S_plot),
            fill = col,
            border = border,
            bty = "n",
            cex = 0.8
        )
    }

    invisible(list(
        compositions = S_plot,
        row_order = row_order,
        col_order = col_order,
        row_hclust = row_hclust,
        col_hclust = col_hclust
    ))
}

#' Plot method for archetypes objects
#'
#' Draws diagnostic and geometric plots for a fitted archetypal analysis model.
#' The `what` argument controls which part of the fit is visualized:
#'
#' * `"composition"` and `"compositions"` draw a stacked barplot of
#'   `x$compositions`, with bars representing observations and colors
#'   representing archetypes. Rows and columns are clustered by default.
#' * `"ternary"` and `"simplex"` plot the rows of `x$compositions` as points in
#'   simplex coordinates. For three archetypes this is a ternary plot: points
#'   near a corner are dominated by one archetype, points near an edge mix two
#'   archetypes, and points near the center mix all three.
#' * `"loss"` plots the objective value stored in `x$loss$loss` across
#'   optimization iterations. This is useful for checking whether the fitting
#'   algorithm reduced the fitted objective and whether the loss curve has plateaued.
#' * `"profiles"` plots the fitted archetypes in their natural representation:
#'   coefficient functions for `fd` fits, and archetype profiles for matrix and
#'   probabilistic fits.
#' * `"coordinates"` plots the original observations together with the fitted
#'   archetype coordinates in feature space. In two dimensions the archetypes
#'   are connected as a closed polygon. In more than two dimensions the method
#'   draws pairwise scatterplots and connects archetypes with dashed closed
#'   paths in each panel. With `projection = "pca"`, higher-dimensional data and
#'   archetypes are first projected to the first two principal components.
#'
#' @param x An object of class `archetypes`
#' @param what Character string naming the plot to draw. Supported values are
#'   `"composition"`, `"compositions"`, `"ternary"`, `"simplex"`, `"loss"`,
#'   `"profiles"`, and `"coordinates"`. `"composision"` and `"composisions"` are
#'   accepted as aliases for `"compositions"`.
#' @param samples Optional sample subset for plots that display observations:
#'   `"composition"`, `"compositions"`, `"ternary"`, `"simplex"`, and
#'   `"coordinates"`. May be numeric row indices, sample names, or a logical
#'   vector. Subsetting is applied before clustering or projection.
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
                            what = c("compositions", "loss", "coordinates", "profiles"),
                            samples = NULL,
                            data = NULL,
                            projection = c("none", "pca"),
                            ...) {
    stopifnot(inherits(x, "archetypes"))

    what <- match.arg(
        tolower(what[1L]),
        c("compositions", "composition", "composision", "composisions",
          "ternary", "simplex",
          "loss", "coordinates", "profiles")
    )

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
    subset_samples <- function(z) {
        if (is.null(samples))
            return(z)
        if (is.logical(samples) && length(samples) != nrow(z))
            stop("Logical `samples` must have one value per sample", call. = FALSE)
        if (is.character(samples) && anyNA(match(samples, rownames(z))))
            stop("Some `samples` are not sample names", call. = FALSE)
        out <- z[samples, , drop = FALSE]
        if (nrow(out) == 0L)
            stop("`samples` selects no samples", call. = FALSE)
        out
    }

    if (what == "compositions") {
        S <- subset_samples(as.matrix(x[["compositions"]]))
        args <- plot_args(
            list(x = S, cluster_rows = TRUE, cluster_cols = TRUE),
            dots
        )
        do.call(composition_barplot, args)
        return(invisible(x))
    }

    # Composition Ternary/Simplex Plot ---------------------------------------

    if (what %in% c("ternary", "simplex")) {
        S <- subset_samples(as.matrix(x[["compositions"]]))
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
        if (is.null(loss[["loss"]]))
            stop("`x$loss` must contain a `loss` column for loss plots")
        args <- plot_args(
            list(
                x = seq_len(nrow(loss)) - 1L,
                y = loss[["loss"]],
                type = "l",
                xlab = "Iteration",
                ylab = "Loss"
            ),
            dots
        )
        do.call(plot, args)
        return(invisible(x))
    }

    # Profile Plot ----------------------------------------------------------

    if (what == "profiles") {
        if (inherits(x[["data"]], "fd")) {
            plot(coordinates_fd(x), ...)
            return(invisible(x))
        }
        .aa_plot_profiles(x, ...)
        return(invisible(x))
    }

    # Coordinate Plot -------------------------------------------------------

    family <- x[["family"]]
    if (is.null(family))
        family <- "gaussian"
    if (!(family %in% c("gaussian", "directional"))) {
        msg <- paste(
            "Coordinate plots are not defined when archetype coordinates live",
            "in parameter space but `data` lives in observation space;",
            sprintf("family = '%s'", family)
        )
        stop(msg, call. = FALSE)
    }

    X <- if (is.null(data)) x[["data"]] else data
    if (is.null(X)) {
        msg <- paste("Original data must be provided either when constructing",
                     "the archetypes object or through `data` for coordinate plots")
        stop(msg)
    }

    X <- if (inherits(X, "fd")) .aa_fd_to_matrix(X) else as.matrix(X)
    X <- subset_samples(X)
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

.aa_plot_profiles <- function(x, ...) {
    A <- as.matrix(x[["coordinates"]])
    if (ncol(A) < 1L)
        stop("Profile plots require at least one feature.", call. = FALSE)

    dots <- list(...)
    arg_or <- function(name, default) {
        value <- dots[[name]]
        if (is.null(value)) default else value
    }
    dots[c("x", "y", "type", "lty", "col", "xlab", "ylab", "xaxt")] <- NULL

    family <- x[["family"]]
    if (is.null(family))
        family <- "gaussian"
    ylab <- arg_or(
        "ylab",
        if (identical(family, "gaussian")) "Value" else sprintf("%s parameter", family)
    )
    feature_names <- colnames(A)
    feature_x <- seq_len(ncol(A))
    xlab <- arg_or("xlab", "Feature")
    col <- arg_or("col", seq_len(nrow(A)))
    lty <- arg_or("lty", 1)

    args <- c(
        list(
            x = feature_x,
            y = t(A),
            type = "l",
            lty = lty,
            col = col,
            xlab = xlab,
            ylab = ylab,
            xaxt = if (is.null(feature_names)) "s" else "n"
        ),
        dots
    )
    do.call(graphics::matplot, args)
    if (!is.null(feature_names)) {
        graphics::axis(1, at = feature_x, labels = feature_names, las = 2)
    }
    invisible(x)
}

# Directional Archetypes Class -----------------------------------------------

#' Directional archetype analysis result object
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
#'
#' @export
directional_archetypes <- function(coordinates,
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
                                   precision = NULL) {
    out <- archetypes(
        coordinates = coordinates,
        coefficients = coefficients,
        compositions = compositions,
        slack = 0,
        loss = loss,
        converged = converged,
        call = call,
        data = data,
        init = init,
        family = "directional"
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
    Y <- object[["compositions"]] %*% object[["coordinates"]]
    Y <- .aa_unit_rows(Y)
    X <- object[["data"]]
    if (!is.null(X)) {
        X <- .aa_unit_rows(as.matrix(X))
        Y <- .aa_align_rows(Y, X)
    }
    Y
}

#' @exportS3Method
residuals.directional_archetypes <- function(object, data = NULL, ...) {
    X <- if (is.null(data)) object[["data"]] else data
    if (is.null(X)) {
        msg <- paste("Original data must be provided either when constructing",
                     "the directional archetypes object or as an argument",
                     "to `residuals.directional_archetypes()`")
        stop(msg)
    }
    X <- .aa_unit_rows(as.matrix(X))
    Y <- .aa_unit_rows(object[["compositions"]] %*% object[["coordinates"]])
    Y <- .aa_align_rows(Y, X)
    X - Y
}

#' @exportS3Method
predict.directional_archetypes <- function(object,
                                           newdata,
                                           max_iter = 100L,
                                           eps = 1e-8,
                                           step_size = 1.0,
                                           max_iter_optimizer = 10L,
                                           step_shrinkage = 0.5,
                                           ...) {
    # Prediction fixes the learned directional archetypes A and solves only for
    # S, the simplex composition of each new row in that archetype hull. This is
    # the same S subproblem used in the alternating training loop, but without
    # the B/C update because new data must not change the fitted archetypes.
    A <- object[["coordinates"]]
    X <- as.matrix(newdata)
    if (ncol(X) != ncol(A)) {
        fmt <- "`newdata` has %d columns but `object$coordinates` has %d columns"
        stop(sprintf(fmt, ncol(X), ncol(A)))
    }
    .aa_check_no_zero_rows(X)
    X_loss <- if (identical(object[["precision"]], "unit")) .aa_unit_rows(X) else X
    .aa_directional_fit_S(
        X_loss = X_loss,
        A = A,
        max_iter = max_iter,
        eps = eps,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage
    )
}

#' @exportS3Method
AIC.directional_archetypes <- function(object, ...) {
    stop("AIC is not defined for Watson-loss directional archetypes.", call. = FALSE)
}

# Kernel Archetypes Class -----------------------------------------------------

#' Kernel archetype analysis result object
#'
#' @param coefficients Numeric matrix (`K x N`) giving each Hilbert-space
#'   archetype as a weighted combination of training samples.
#' @param compositions Numeric matrix (`N x K`) giving each sample as a weighted
#'   combination of kernel archetypes.
#' @param gram Training Gram matrix.
#' @param coordinates_proxy Optional input-space proxy coordinates
#'   `coefficients %*% data`.
#' @param slack Non-negative coefficient row-sum relaxation.
#' @param loss Data frame containing per-iteration metrics.
#' @param converged Logical convergence flag.
#' @param call Matched call.
#' @param data Optional original data matrix.
#' @param init Optional initial coefficient matrix.
#' @param kernel Kernel specification.
#' @param kernel_args Kernel arguments.
#'
#' @export
kernel_archetypes <- function(coefficients,
                              compositions,
                              gram,
                              coordinates_proxy = NULL,
                              slack = 0,
                              loss = NULL,
                              converged = TRUE,
                              call = NULL,
                              data = NULL,
                              init = NULL,
                              kernel = NULL,
                              kernel_args = list()) {
    if (is.null(call))
        call <- match.call()
    if (is.null(loss))
        loss <- data.frame(loss = NA_real_, r2 = NA_real_, k_S = NA_real_, k_A = NA_real_)
    K <- nrow(coefficients)
    N <- ncol(coefficients)
    stopifnot("nrow(compositions) must match number of samples" = nrow(compositions) == N)
    stopifnot("ncol(compositions) must match number of archetypes" = ncol(compositions) == K)
    stopifnot("Gram matrix dimensions must match number of samples" =
                  identical(dim(gram), c(N, N)))
    stopifnot("All slack values must be non-negative" = all(slack >= 0))
    if (any(slack > 0)) {
        a <- rowSums(coefficients)
        stopifnot("Some rowSums(coefficients) are above allowed slack" = all(a <= 1 + slack))
        stopifnot("Some rowSums(coefficients) are below allowed slack" = all(a >= 1 - slack))
    } else if (!isTRUE(all.equal(rowSums(coefficients), rep(1, K), check.attributes = FALSE))) {
        stop("Coefficients must be row-stochastic (each row sums to 1) when slack = 0")
    }
    if (!isTRUE(all.equal(rowSums(compositions), rep(1, N), check.attributes = FALSE)))
        stop("Compositions must be row-stochastic (each row sums to 1)")
    if (!is.null(coordinates_proxy)) {
        stopifnot("coordinates_proxy must have one row per archetype" =
                      nrow(coordinates_proxy) == K)
    }

    structure(
        list(
            coefficients = coefficients,
            compositions = compositions,
            gram = gram,
            coordinates_proxy = coordinates_proxy,
            slack = slack,
            init = init,
            loss = loss,
            converged = converged,
            data = data,
            kernel = kernel,
            kernel_args = kernel_args,
            call = call
        ),
        class = "kernel_archetypes"
    )
}

#' @exportS3Method
coefficients.kernel_archetypes <- function(object, ...)
    object[["coefficients"]]

#' @rdname archetypes
#' @exportS3Method
anames.kernel_archetypes <- function(x)
    rownames(x[["coefficients"]])

#' @rdname archetypes
#' @method anames<- kernel_archetypes
#' @export
`anames<-.kernel_archetypes` <- function(x, value) {
    K <- nrow(x[["coefficients"]])
    if (length(value) != K) {
        fmt <- "Expected %d archetype names, got %d"
        stop(sprintf(fmt, K, length(value)))
    }
    stopifnot("Archetype names must not be missing" = !any(is.na(value)))
    stopifnot("Archetype names must not be empty" = all(nzchar(value)))
    stopifnot("Archetype names must be unique" = !anyDuplicated(value))

    rownames(x[["coefficients"]]) <- value
    colnames(x[["compositions"]]) <- value
    if (!is.null(x[["coordinates_proxy"]]))
        rownames(x[["coordinates_proxy"]]) <- value
    if (!is.null(x[["init"]]))
        rownames(x[["init"]]) <- value
    x
}

#' @exportS3Method
fitted.kernel_archetypes <- function(object, ...) {
    stop(
        "`fitted()` is not defined for nonlinear kernel archetypes; use ",
        "`residuals()` for Hilbert-space residual norms or `coordinates_proxy` ",
        "for input-space visualization.",
        call. = FALSE
    )
}

#' @exportS3Method
residuals.kernel_archetypes <- function(object, ...) {
    G <- object[["gram"]]
    H <- object[["coefficients"]]
    S <- object[["compositions"]]
    AAt <- H %*% G %*% t(H)
    XAt <- G %*% t(H)
    out <- pmax(diag(G) - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt)), 0)
    names(out) <- rownames(S)
    out
}

#' @exportS3Method
print.kernel_archetypes <- function(x, ...) {
    call_str <- paste(deparse(x[["call"]]), sep = "\n", collapse = "\n")
    cat("\nCall:\n", call_str, "\n\n", sep = "")
    cat("Kernel Archetypes Summary:\n")
    cat("Number of Archetypes:", nrow(x[["coefficients"]]), "\n")
    cat("Number of Samples:", ncol(x[["coefficients"]]), "\n")
    if (!is.null(x[["coordinates_proxy"]]))
        cat("Input-space proxy coordinates: available\n")
    loss <- x[["loss"]]
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

#' @exportS3Method
plot.kernel_archetypes <- function(x,
                                   what = c("compositions", "loss", "coordinates", "profiles"),
                                   data = NULL,
                                   projection = c("none", "pca"),
                                   ...) {
    what <- match.arg(
        tolower(what[1L]),
        c("compositions", "composition", "composision", "composisions",
          "ternary", "simplex",
          "loss", "coordinates", "profiles")
    )
    if (what %in% c("composition", "composision", "composisions"))
        what <- "compositions"
    if (what == "compositions") {
        S <- as.matrix(x[["compositions"]])
        dots <- list(...)
        args <- list(x = S, cluster_rows = TRUE, cluster_cols = TRUE)
        args[names(dots)] <- dots
        do.call(composition_barplot, args)
        return(invisible(x))
    }
    if (what %in% c("ternary", "simplex")) {
        S <- as.matrix(x[["compositions"]])
        if (is.null(colnames(S)))
            colnames(S) <- paste0("A", seq_len(ncol(S)))
        graphics::plot(compositions::acomp(S), axes = TRUE, ...)
        return(invisible(x))
    }
    if (what == "loss") {
        loss <- x[["loss"]]
        graphics::plot(
            seq_len(nrow(loss)) - 1L,
            loss[["loss"]],
            type = "l",
            xlab = "Iteration",
            ylab = "Loss",
            ...
        )
        return(invisible(x))
    }
    if (what == "profiles") {
        stop(
            "Profile plots are not defined for kernel archetypes because ",
            "their natural coordinates live in implicit Hilbert space.",
            call. = FALSE
        )
    }

    X <- if (is.null(data)) x[["data"]] else data
    A <- x[["coordinates_proxy"]]
    if (is.null(X) || is.null(A)) {
        stop(
            "Coordinate plots for kernel archetypes require original `data` ",
            "and available `coordinates_proxy`.",
            call. = FALSE
        )
    }
    proxy <- archetypes(
        coordinates = A,
        coefficients = x[["coefficients"]],
        compositions = x[["compositions"]],
        slack = x[["slack"]],
        loss = x[["loss"]],
        converged = x[["converged"]],
        call = x[["call"]],
        data = X,
        init = if (!is.null(x[["init"]]) && !is.null(x[["data"]])) {
            x[["init"]] %*% as.matrix(x[["data"]])
        } else {
            NULL
        }
    )
    plot(proxy, what = "coordinates", data = X, projection = projection, ...)
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
    X <- if (inherits(X, "fd")) .aa_fd_to_matrix(X) else X

    # Compute AIC
    K     <- nrow(object[["coordinates"]])        # Number of Archetypes
    nelem <- prod(dim(X))                         # number of elements in X
    family <- object[["family"]]
    if (is.null(family))
        family <- "gaussian"
    if (!identical(family, "gaussian"))
        stop("AIC is not defined for non-Gaussian archetypes objects.", call. = FALSE)
    rss   <- object[["loss"]][["loss"]]
    X_hat <- with(object, compositions %*% coordinates)
    aic   <- log(rss[length(rss)] / nelem) + 2 * (2*K - 1) / effic(X, X_hat)

    object[["AIC"]] <- aic  # cache for future calls
    return(aic)
}

.aa_fd_to_matrix <- function(x) {
    if (!requireNamespace("fda", quietly = TRUE))
        stop("Package `fda` is required for fd data.", call. = FALSE)
    t(stats::coef(x))
}

.aa_fd_from_rows <- function(X, template, curve_names = NULL) {
    if (!requireNamespace("fda", quietly = TRUE))
        stop("Package `fda` is required for fd data.", call. = FALSE)

    fdnames <- template[["fdnames"]]
    if (is.null(curve_names))
        curve_names <- rownames(X)
    if (!is.null(curve_names))
        fdnames[["reps"]] <- curve_names

    fda::fd(t(X), template[["basis"]], fdnames = fdnames)
}

#' Convert archetype coordinates to functional data
#'
#' Reconstructs the fitted archetype coordinates as an `fda::fd` object. This
#' is mainly useful for models fitted with [run_aa()] on `fda::fd` input, where
#' the optimizer stores coordinates as a coefficient matrix.
#'
#' @param object An object of class \code{\link{archetypes}}.
#' @param data Optional `fda::fd` object to use as the source for the basis and
#'   functional data names. Defaults to `object$data`.
#' @param basis Optional `fda` basis object. Used when `data` is not available.
#' @param fdnames Optional `fdnames` list passed to `fda::fd()`. Defaults to
#'   the names from `data`, with archetype names used as `reps`.
#' @param ... Ignored.
#'
#' @return An `fda::fd` object whose columns are the archetype functions.
#'
#' @export
coordinates_fd <- function(object, data = object[["data"]], basis = NULL, fdnames = NULL, ...) {
    stopifnot(inherits(object, "archetypes"))
    if (!requireNamespace("fda", quietly = TRUE))
        stop("Package `fda` is required for `coordinates_fd()`.", call. = FALSE)

    if (!is.null(data) && !inherits(data, "fd"))
        stop("`data` must be an `fda::fd` object.", call. = FALSE)

    if (is.null(basis))
        basis <- data[["basis"]]
    if (is.null(basis))
        stop("No fd basis found. Supply `data` or `basis` explicitly.", call. = FALSE)

    A <- object[["coordinates"]]
    if (is.null(fdnames))
        fdnames <- data[["fdnames"]]
    if (is.null(fdnames))
        fdnames <- list(args = "time", reps = rownames(A), funs = "values")
    fdnames[["reps"]] <- rownames(A)

    fda::fd(t(A), basis, fdnames = fdnames)
}
