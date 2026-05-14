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
#' @param weights Optional numeric vector of sample weights used during
#'   fitting. When present, it is used by `residuals(..., type = "pearson")`.
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
                       family = "gaussian",
                       weights = NULL) {
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
            converged    = converged,
            call         = call,
            data         = data,
            init         = init,
            family       = family,
            weights      = weights
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
#' @param call The matched function call that created the object.
#' @param data Optional original data.
#' @param init Optional initial archetype coordinates.
#' @param family Observation family.
#' @param weights Optional numeric vector of sample weights used during
#'   fitting.
new_archetypes <- function(coordinates,
                           coefficients,
                           compositions,
                           slack,
                           loss,
                           converged,
                           call = NULL,
                           data = NULL,
                           init = NULL,
                           family = "gaussian",
                           weights = NULL) {


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
    if (!is.null(weights)) {
        stopifnot("weights must be a numeric vector" = is.numeric(weights))
        stopifnot("weights must have one value per sample" = length(weights) == N)
        stopifnot("weights must contain no missing values" = !any(is.na(weights)))
        stopifnot("weights must be non-negative" = all(weights >= 0))
    }

    structure(
        list(
            coordinates  = coordinates,
            coefficients = coefficients,
            compositions = compositions,
            slack        = slack,
            init         = init,
            loss         = loss,
            converged    = converged,
            data         = data,
            call         = call,
            family       = family,
            weights      = weights
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
#' @param type Residual type. `"response"` (default) returns `data - fitted`.
#'   `"pearson"` returns row-wise weighted residuals,
#'   `sqrt(weights) * (data - fitted)`.
#' @param ... Ignored.
#'
#' @return
#' Residuals in the same representation as `data`.
#'
#' @exportS3Method
residuals.archetypes <- function(object,
                                 type = c("response", "pearson"),
                                 ...) {
    type <- match.arg(type)
    X <- object[["data"]]
    if (is.null(X))
        stop("Original data is missing from `object$data`.", call. = FALSE)

    X_hat <- fitted(object)
    if (inherits(X, "fd"))
        return(X - X_hat)

    stopifnot(all(dim(X) == dim(X_hat)))
    out <- X - X_hat
    weights <- object[["weights"]]

    if (identical(type, "response") || is.null(weights))
        return(out)
    sqrt(weights) * out
}

#' Predict compositions or reconstructions for new data from an archetypes model
#'
#' Projects new samples onto the archetype space by solving for composition
#' weights with fixed archetype coordinates. By default the method returns
#' compositions (`S`). With `type = "reconstruction"` it returns
#' `X_hat = S %*% A`, where `A = object$coordinates`.
#'
#' For non-Gaussian families (`bernoulli`, `poisson`, `multinomial`),
#' coordinates live in parameter space. In particular, multinomial
#' reconstructions preserve existing package behavior by scaling each row of
#' `S %*% A` by the corresponding row total of `newdata`.
#'
#' @param object An object of class `archetypes`.
#' @param newdata New data to fit. Must contain the features (columns) used to
#'   fit `object`; fd-backed fits may pass an `fda::fd` object.
#' @param type Prediction output type.
#'   `"compositions"` (default) returns `S`; `"reconstruction"` returns
#'   `X_hat = S %*% A`.
#' @param ... Passed to \code{\link{fit_simplex}} and family-specific
#'   prediction routines.
#'
#' @return
#' If `type = "compositions"`, a numeric matrix (N_new x K) with non-negative
#' row-stochastic composition weights. If `type = "reconstruction"`, a
#' reconstruction matrix (N_new x M) in the same feature space as
#' `object$coordinates`.
#'
#' @seealso [archetypes()], [fit_simplex()]
#'
#' @exportS3Method
predict.archetypes <- function(object,
                               newdata,
                               type = c("compositions", "reconstruction"),
                               ...) {
    type <- match.arg(type)
    A <- object[["coordinates"]]
    is_fd_newdata <- inherits(newdata, "fd")
    if (is_fd_newdata)
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

    S <- if (identical(family, "gaussian")) {
        fit_simplex(A, X, ...)
    } else {
        .aa_paa_predict_S(object, X, ...)
    }

    if (identical(type, "compositions"))
        return(S)

    X_hat <- S %*% A
    if (identical(family, "multinomial")) {
        totals <- rowSums(as.matrix(X))
        X_hat <- totals * X_hat
    }
    if (is_fd_newdata && inherits(object[["data"]], "fd"))
        return(.aa_fd_from_rows(X_hat, object[["data"]], rownames(X_hat)))
    X_hat
}

#' @exportS3Method
print.archetypes <- function(x, ...) {
    # TODO: add more info about method used (method, missing, robust, etc.)
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
#' @param cluster_rows,cluster_cols Logical values, one of `"PC1"` or
#'   `"AOP"`, or `hclust` objects. When `TRUE`, rows or columns are reordered
#'   by hierarchical clustering. `"PC1"` orders rows by their first principal
#'   component score and columns by their first principal component loading.
#'   `"AOP"` (angle of projection) orders rows by the angle of their PC1/PC2
#'   scores and columns by the angle of their PC1/PC2 loadings. Row clustering
#'   is computed on an internal centered log-ratio (`clr`) transform.
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
        if (distance == "spearman")
            return(stats::as.dist(1 - stats::cor(t(data), method = "spearman")))
        stats::dist(data, method = distance)
    }

    make_component_order <- function(components, mode) {
        if (mode == "PC1")
            return(order(components[, 1L]))
        if (ncol(components) < 2L)
            return(NULL)
        order(atan2(components[, 2L], components[, 1L]))
    }

    make_hclust <- function(value, data, margin, distance) {
        cluster_arg <- if (margin == "rows") "cluster_rows" else "cluster_cols"
        if (inherits(value, "hclust"))
            return(value)
        if (is.character(value) && length(value) == 1L) {
            mode <- toupper(value)
            if (!(mode %in% c("PC1", "AOP")))
                stop(
                    sprintf("`%s` must be logical, 'PC1', 'AOP', or an hclust object", cluster_arg),
                    call. = FALSE
                )
            pca <- stats::prcomp(
                data,
                center = TRUE,
                scale. = FALSE,
                rank. = if (mode == "AOP") 2L else 1L
            )
            components <- if (margin == "rows") pca[["x"]] else pca[["rotation"]]
            order <- make_component_order(components, mode)
            if (is.null(order))
                return(NULL)
            return(list(order = order))
        }
        if (!is.logical(value))
            stop(
                sprintf("`%s` must be logical, 'PC1', 'AOP', or an hclust object", cluster_arg),
                call. = FALSE
            )
        if (!isTRUE(value))
            return(NULL)
        if (margin == "rows") {
            if (nrow(data) < 2L)
                return(NULL)
            transformed <- .aa_clr(data)
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
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    if (isTRUE(legend))
        graphics::par(mar = old_par[["mar"]] + c(0, 0, 0, 4))

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
        usr <- graphics::par("usr")
        legend_x <- usr[2L] + 0.03 * diff(usr[1:2])
        graphics::legend(
            x = legend_x,
            y = usr[4L],
            legend = colnames(S_plot),
            fill = col,
            border = border,
            bty = "n",
            cex = 0.8,
            xpd = NA,
            xjust = 0,
            yjust = 1
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
#' @param show_anames Logical; for `what = "coordinates"`, whether
#'   to draw archetype names next to archetype points.
#' @param args.archetypes Optional named list of graphical arguments used for
#'   archetype points/paths in `what = "coordinates"`. Useful for styling
#'   archetypes separately when `...` aesthetics (e.g., `col`) are mapped to
#'   observation-level vectors.
#' @param ... Additional graphical parameters passed to the underlying plotting
#'   functions.
#'
#' @importFrom graphics lines pairs plot points
#' @exportS3Method
plot.archetypes <- function(x,
                            what = c("compositions", "loss", "coordinates", "profiles"),
                            samples = NULL,
                            data = NULL,
                            projection = c("none", "pca"),
                            show_anames = TRUE,
                            args.archetypes = list(),
                            ...) {
    stopifnot(inherits(x, "archetypes"))
    if (!is.list(args.archetypes))
        stop("`args.archetypes` must be a list", call. = FALSE)

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
        if (!requireNamespace("compositions", quietly = TRUE)) {
            stop(
                "Ternary/simplex plots require package `compositions`. ",
                "Install it to use `plot(..., what = 'ternary'/'simplex')`.",
                call. = FALSE
            )
        }
        S <- subset_samples(as.matrix(x[["compositions"]]))
        if (is.null(colnames(S)))
            colnames(S) <- paste0("A", seq_len(ncol(S)))
        main <- dots[["main"]]
        args <- plot_args(
            list(x = compositions::acomp(S), axes = TRUE),
            dots[setdiff(names(dots), "main")]
        )
        do.call(plot, args)
        if (!is.null(main)) graphics::title(main = main)
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
        plot(
            x_projected,
            what = "coordinates",
            projection = "none",
            show_anames = show_anames,
            args.archetypes = args.archetypes,
            ...
        )
        return(invisible(x))
    }

    split_style <- function(name, data_default, archetype_default) {
        data_value <- dots[[name]]
        archetype_value <- args.archetypes[[name]]
        if (is.null(data_value))
            data_value <- data_default
        if (is.null(archetype_value))
            archetype_value <- archetype_default
        list(data = data_value, archetype = archetype_value)
    }

    col_style <- split_style("col", "black", "red")
    pch_style <- split_style("pch", 1, 16)
    cex_style <- split_style("cex", 1, 1.3)
    lwd_style <- split_style("lwd", 1, 1)

    data_col <- col_style[["data"]]
    arch_col <- col_style[["archetype"]]
    data_pch <- pch_style[["data"]]
    arch_pch <- pch_style[["archetype"]]
    data_cex <- cex_style[["data"]]
    arch_cex <- cex_style[["archetype"]]
    arch_lwd <- lwd_style[["archetype"]]
    arch_lty_2d <- if (!is.null(args.archetypes[["lty"]])) args.archetypes[["lty"]] else 1
    arch_lty_pairs <- if (!is.null(args.archetypes[["lty"]])) args.archetypes[["lty"]] else 2

    dots[["col"]] <- data_col
    dots[["pch"]] <- data_pch
    dots[["cex"]] <- data_cex

    if (is.null(colnames(X)))
        colnames(X) <- paste0("V", seq_len(ncol(X)))
    if (is.null(colnames(A)))
        colnames(A) <- colnames(X)
    archetype_names <- anames(x)

    if (ncol(X) == 2L) {
        args <- plot_args(
            list(
                x = X[, 1L],
                y = X[, 2L],
                xlab = colnames(X)[1L],
                ylab = colnames(X)[2L],
                asp = 1
            ),
            dots
        )
        do.call(plot, args)
        A_closed <- A[c(seq_len(nrow(A)), 1L), , drop = FALSE]
        lines(A_closed[, 1L], A_closed[, 2L], col = arch_col, lwd = arch_lwd, lty = arch_lty_2d)
        points(A[, 1L], A[, 2L], col = arch_col, pch = arch_pch, cex = arch_cex)
        if (isTRUE(show_anames) && !is.null(archetype_names)) {
            graphics::text(
                x = A[, 1L],
                y = A[, 2L],
                labels = archetype_names,
                pos = 4,
                cex = 0.8,
                col = arch_col,
                xpd = NA
            )
        }
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
         lines(x[arch_closed], y[arch_closed], col = arch_col, lwd = arch_lwd, lty = arch_lty_pairs)
        points(x[arch_ix], y[arch_ix],
               col = arch_col, pch = arch_pch, cex = arch_cex)
        if (isTRUE(show_anames) && !is.null(archetype_names)) {
            graphics::text(
                x = x[arch_ix],
                y = y[arch_ix],
                labels = archetype_names,
                pos = 4,
                cex = 0.7,
                col = arch_col,
                xpd = NA
            )
        }
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
    family <- x[["family"]]
    if (is.null(family))
        family <- "gaussian"
    ylab <- arg_or(
        "ylab",
        if (identical(family, "gaussian")) "Value" else sprintf("%s parameter", family)
    )
    xlab <- arg_or("xlab", "Feature")
    col <- arg_or("col", NULL)
    legend_text <- arg_or("legend.text", anames(x))
    args_legend <- arg_or("args.legend", list(bty = "n"))

    args <- list(
        height = A,
        beside = TRUE,
        col = col,
        legend.text = legend_text,
        args.legend = args_legend,
        xlab = xlab,
        ylab = ylab
    )
    dots[["height"]] <- NULL
    args[names(dots)] <- dots
    do.call(graphics::barplot, args)
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
residuals.directional_archetypes <- function(object, ...) {
    X <- object[["data"]]
    if (is.null(X))
        stop("Original data is missing from `object$data`.", call. = FALSE)

    X <- .aa_unit_rows(as.matrix(X))
    Y <- .aa_unit_rows(object[["compositions"]] %*% object[["coordinates"]])
    Y <- .aa_align_rows(Y, X)
    X - Y
}

#' Predict compositions or reconstructions for directional archetypes
#'
#' Solves for new sample compositions with fixed directional archetype
#' coordinates. By default returns compositions (`S`).
#' With `type = "reconstruction"`, returns directional reconstructions obtained
#' from `S %*% A`, then row-normalized and hemisphere-aligned to `newdata`.
#'
#' @param object An object of class `directional_archetypes`.
#' @param newdata New directional data matrix.
#' @param type Prediction output type. `"compositions"` (default) or
#'   `"reconstruction"`.
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
                                           type = c("compositions", "reconstruction"),
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
    A <- object[["coordinates"]]
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
    if (identical(type, "compositions"))
        return(S)

    Y <- .aa_unit_rows(S %*% A)
    .aa_align_rows(Y, .aa_unit_rows(X))
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
#' @param coordinates Optional input-space coordinates `coefficients %*% data`.
#'   For nonlinear kernels these are visualization coordinates, not exact
#'   Hilbert-space archetypes.
#' @param slack Non-negative coefficient row-sum relaxation.
#' @param loss Data frame containing per-iteration metrics.
#' @param converged Logical convergence flag.
#' @param call Matched call.
#' @param data Optional original data matrix.
#' @param init Optional initial coefficient matrix.
#' @param kernel Kernel specification.
#' @param kernel_args Kernel arguments.
#' @param weights Optional numeric vector of sample weights used during
#'   fitting.
#'
#' @export
kernel_archetypes <- function(coefficients,
                              compositions,
                              gram,
                              coordinates = NULL,
                              slack = 0,
                              loss = NULL,
                              converged = TRUE,
                              call = NULL,
                              data = NULL,
                              init = NULL,
                              kernel = NULL,
                              kernel_args = list(),
                              weights = NULL) {
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
    if (!is.null(coordinates)) {
        stopifnot("coordinates must have one row per archetype" =
                      nrow(coordinates) == K)
    }
    if (!is.null(weights)) {
        stopifnot("weights must be a numeric vector" = is.numeric(weights))
        stopifnot("weights must have one value per sample" = length(weights) == N)
        stopifnot("weights must contain no missing values" = !any(is.na(weights)))
        stopifnot("weights must be non-negative" = all(weights >= 0))
    }

    structure(
        list(
            coefficients = coefficients,
            compositions = compositions,
            gram = gram,
            coordinates = coordinates,
            slack = slack,
            init = init,
            loss = loss,
            converged = converged,
            data = data,
            kernel = kernel,
            kernel_args = kernel_args,
            weights = weights,
            call = call
        ),
        class = "kernel_archetypes"
    )
}

#' @exportS3Method
coefficients.kernel_archetypes <- function(object, ...)
    object[["coefficients"]]

#' Predict method for kernel archetypes
#'
#' Out-of-sample prediction is not currently defined for `kernel_archetypes`.
#' The fitted object stores training compositions and optional input-space
#' coordinates for visualization, but projecting `newdata` requires kernel-
#' specific cross-Gram evaluation that is not part of this API.
#'
#' @param object An object of class `kernel_archetypes`.
#' @param newdata New data to project (currently unsupported).
#' @param type Prediction output type. `"compositions"` or `"reconstruction"`.
#' @param ... Ignored.
#'
#' @return No return value. Always errors with an explanatory message.
#'
#' @exportS3Method
predict.kernel_archetypes <- function(object,
                                      newdata,
                                      type = c("compositions", "reconstruction"),
                                      ...) {
    type <- match.arg(type)
    stop(
        "predict() is not currently defined for kernel_archetypes; ",
        "out-of-sample projection requires kernel-specific cross-Gram ",
        "evaluation. Use object[['compositions']] for training compositions.",
        call. = FALSE
    )
}

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
    if (!is.null(x[["coordinates"]]))
        rownames(x[["coordinates"]]) <- value
    if (!is.null(x[["init"]]))
        rownames(x[["init"]]) <- value
    x
}

#' @exportS3Method
fitted.kernel_archetypes <- function(object, ...) {
    stop(
        "`fitted()` is not defined for nonlinear kernel archetypes; use ",
        "`residuals()` for Hilbert-space residual norms or `coordinates` ",
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
    if (!is.null(x[["coordinates"]]))
        cat("Input-space coordinates: available\n")
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
        args <- list(x = S, cluster_rows = TRUE, cluster_cols = TRUE, linkage = "ward.D2")
        args[names(dots)] <- dots
        do.call(composition_barplot, args)
        return(invisible(x))
    }
    if (what %in% c("ternary", "simplex")) {
        if (!requireNamespace("compositions", quietly = TRUE)) {
            stop(
                "Ternary/simplex plots require package `compositions`. ",
                "Install it to use `plot(..., what = 'ternary'/'simplex')`.",
                call. = FALSE
            )
        }
        S <- as.matrix(x[["compositions"]])
        if (is.null(colnames(S)))
            colnames(S) <- paste0("A", seq_len(ncol(S)))
        dots2 <- list(...)
        main <- dots2[["main"]]
        do.call(
            graphics::plot,
            c(list(compositions::acomp(S), axes = TRUE),
              dots2[setdiff(names(dots2), "main")])
        )
        if (!is.null(main)) graphics::title(main = main)
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
    A <- x[["coordinates"]]
    if (is.null(X) || is.null(A)) {
        stop(
            "Coordinate plots for kernel archetypes require original `data` ",
            "and available `coordinates`.",
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
#' Computes the AIC-like validity criterion for an `archetypes` object.
#'
#' This is not the classical likelihood-based Akaike Information Criterion
#' \eqn{-2 \log L + 2k}. Instead, it implements the adapted archetypal-analysis
#' criterion proposed by Suleman (2017), using the reconstruction variance
#' \eqn{\|X - \hat X\|_F^2 / (N M)} and an efficiency-adjusted complexity
#' penalty based on the full parameter count
#' \eqn{K_\mu + K_\beta + 1 = N(K - 1) + K(N - 1) + 1}. The value is intended
#' for comparing Euclidean Gaussian archetype fits on the same data, especially
#' across different numbers of archetypes.
#' Following the assumptions in Suleman (2017), the criterion is undefined when
#' the number of samples is not larger than the number of features. If the
#' covariance matrix of `X` is singular otherwise, the efficiency term is
#' computed with a Moore-Penrose pseudo-inverse and a warning is emitted.
#'
#' @param object An object of class `archetypes`.
#' @param ... Ignored.
#' @return Numeric scalar with the adapted AIC-like criterion.
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

    X <- object[["data"]]
    if (is.null(X))
        stop(paste("AIC requires original data `X`;",
                   "provide it when constructing the archetypes object."))
    X <- if (inherits(X, "fd")) .aa_fd_to_matrix(X) else X

    # Compute AIC
    K <- nrow(object[["coordinates"]]) # Number of Archetypes
    N <- nrow(X)
    M <- ncol(X)
    nelem <- prod(dim(X)) # number of elements in X
    family <- object[["family"]]
    if (is.null(family))
        family <- "gaussian"
    if (!identical(family, "gaussian"))
        stop("AIC is not defined for non-Gaussian archetypes objects.", call. = FALSE)
    if (N <= M) {
        warning(paste("Adapted AIC is undefined when the number of samples",
                      "is not larger than the number of features; returning NA."),
                call. = FALSE)
        return(NA_real_)
    }
    # rss   <- object[["loss"]][["loss"]]
    X_hat <- with(object, compositions %*% coordinates)
    eta <- tryCatch(effic(X, X_hat), error = function(e) NA_real_)
    if (!is.finite(eta) || eta <= 0) {
        warning(paste("Adapted AIC is undefined because the efficiency term",
                      "is non-positive or non-finite; returning NA."),
                call. = FALSE)
        return(NA_real_)
    }
    npar <- N * (K - 1) + K * (N - 1) + 1  # K_mu + K_beta + 1
    rss  <- norm(X - X_hat, "F")^2  # "loss" is not necessarily the RSS on the original data
    log(rss / nelem) + 2 * npar / (N * eta)
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
