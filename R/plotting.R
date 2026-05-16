#' Composition Plot For Archetypes
#'
#' Draws a horizontal stacked barplot for a matrix-like set of composition
#' weights, with rows interpreted as samples and columns interpreted as
#' archetypes or other compositional parts.
#'
#' @param compositions Numeric matrix or data frame. Rows are samples and
#'   columns are archetypes. Rows should contain non-negative composition
#'   weights.
#' @param plot Logical. Should the plot be drawn?
#' @param cluster_rows,cluster_cols Logical values, one of `"PC1"` or
#'   `"AOP"`, or `hclust` objects. When `TRUE`, rows or columns are reordered
#'   by hierarchical clustering. `"PC1"` orders rows by their first principal
#'   component score and columns by their first principal component loading.
#'   `"AOP"` orders rows by the angle of their PC1/PC2 scores and columns by
#'   the angle of their PC1/PC2 loadings.
#' @param distance Optional distance metric used for both row and column
#'   clustering when `distance_rows` or `distance_cols` are not supplied.
#' @param distance_rows,distance_cols Distance metrics used when clustering
#'   rows or columns. Values may be any method accepted by [stats::dist()],
#'   `"correlation"`, a function that returns a `dist` object, or a precomputed
#'   `dist` object. Row distances are computed after the centered log-ratio
#'   transform.
#' @param linkage Linkage method passed to [stats::hclust()].
#' @param col Optional vector of colors, one per archetype. Defaults to a
#'   qualitative HCL palette.
#' @param legend Logical. Should an archetype legend be drawn?
#' @param border Border color for the stacked bar segments.
#' @param ... Additional graphical parameters passed to [graphics::barplot()].
#'
#' @return Invisibly returns a data frame in long format with one row per
#'   sample/archetype pair. Columns are `sample` (factor ordered by the applied
#'   clustering), `archetype` (factor ordered by the applied clustering), and
#'   `weight` (numeric).
#'
#' @export
plot_archetypes_compositions <- function(compositions,
                                         plot = TRUE,
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
    S <- as.matrix(compositions)
    # Sanity checks
    if (!is.numeric(S))
        stop("`compositions` must be a numeric matrix or data frame", call. = FALSE)
    if (!all(is.finite(S)))
        stop("`compositions` must contain only finite values", call. = FALSE)
    if (any(S < 0))
        stop("`compositions` must contain non-negative composition weights", call. = FALSE)
    if (nrow(S) == 0L || ncol(S) == 0L)
        stop("`compositions` must have at least one row and one column", call. = FALSE)

    # Set defaults
    rownames(S) <- rownames(S) %||% seq_len(nrow(S))
    colnames(S) <- colnames(S) %||% paste0("A", seq_len(ncol(S)))
    if (!is.null(distance)) {
        if (missing(distance_rows))
            distance_rows <- distance
        if (missing(distance_cols))
            distance_cols <- distance
    }

    row_hclust <- .aa_composition_hclust(cluster_rows, S, "rows", distance_rows, linkage)
    col_hclust <- .aa_composition_hclust(cluster_cols, S, "cols", distance_cols, linkage)
    row_order <- row_hclust[["order"]] %||% seq_len(nrow(S))
    col_order <- col_hclust[["order"]] %||% seq_len(ncol(S))
    S_plot <- S[row_order, col_order, drop = FALSE]

    if (is.null(col)) {
        col <- grDevices::hcl.colors(ncol(S_plot), palette = "Dark 3")
    } else if (!is.null(names(col)) && all(colnames(S_plot) %in% names(col))) {
        col <- col[colnames(S_plot)]
    }
    col <- rep_len(col, ncol(S_plot))

    dots <- list(...)
    barplot_args <- list(
        height = t(S_plot),
        col = col,
        border = border,
        xlab = "Composition",
        names.arg = rownames(S_plot),
        horiz = TRUE,
        las = 1L
    ) %|p|% dots

    if (isTRUE(plot)) {
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        if (isTRUE(legend))
            graphics::par(mar = old_par[["mar"]] + c(0, 0, 0, 4))
        do.call(graphics::barplot, barplot_args)
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
    }

    out <- data.frame(
        sample    = rep(rownames(S_plot), times = ncol(S_plot)),
        archetype = rep(colnames(S_plot), each  = nrow(S_plot)),
        weight    = as.vector(S_plot),
        stringsAsFactors = FALSE
    )
    out[["sample"]]    <- factor(out[["sample"]],    levels = rownames(S_plot))
    out[["archetype"]] <- factor(out[["archetype"]], levels = colnames(S_plot))
    invisible(out)
}

#' Loss Plot For Archetypes
#'
#' @param loss Data frame containing a `loss` column.
#' @param plot Logical. Should the plot be drawn?
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns a data frame with columns `iteration` (integer,
#'   0-based) and `loss` (numeric).
#'
#' @export
plot_archetypes_loss <- function(loss, plot = TRUE, ...) {
    if (is.null(loss[["loss"]]))
        stop("`loss` must contain a `loss` column for loss plots", call. = FALSE)
    iterations <- seq_len(nrow(loss)) - 1L
    values <- loss[["loss"]]
    plot_args <- list(
        x = iterations,
        y = values,
        type = "l",
        xlab = "Iteration",
        ylab = "Loss"
    ) %|p|% list(...)
    if (isTRUE(plot))
        do.call(graphics::plot, plot_args)
    invisible(data.frame(iteration = iterations, loss = values))
}

#' Profile Plot For Archetypes
#'
#' @param coordinates Archetype coordinates or an `fda::fd` object.
#' @param family Observation family used to choose a default y-axis label.
#' @param archetype_names Optional archetype labels.
#' @param plot Logical. Should the plot be drawn?
#' @param ... Additional graphical parameters passed to [graphics::barplot()]
#'   for matrix coordinates, or [graphics::plot()] for `fd` coordinates.
#'
#' @return Invisibly returns a data frame in long format with columns
#'   `archetype`, `feature`, and `value`. Returns `NULL` invisibly for `fd`
#'   coordinates.
#'
#' @export
plot_archetypes_profiles <- function(coordinates,
                                     family = "gaussian",
                                     archetype_names = NULL,
                                     plot = TRUE,
                                     ...) {
    if (inherits(coordinates, "fd")) {
        plot_args <- list(...)
        if (isTRUE(plot))
            do.call(graphics::plot, c(list(x = coordinates), plot_args))
        return(invisible(NULL))
    }

    A <- as.matrix(coordinates)
    if (!is.numeric(A))
        stop("`coordinates` must be numeric.", call. = FALSE)
    if (ncol(A) < 1L)
        stop("Profile plots require at least one feature.", call. = FALSE)
    family <- family %||% "gaussian"

    dots <- list(...)
    ylab <- dots[["ylab"]] %||%
        if (identical(family, "gaussian")) "Value" else sprintf("%s parameter", family)
    xlab         <- dots[["xlab"]]          %||% "Feature"
    col          <- dots[["col"]]
    legend_text  <- dots[["legend.text"]]   %||% archetype_names
    args_legend  <- dots[["args.legend"]]   %||% list(bty = "n")

    barplot_args <- list(
        height = A,
        beside = TRUE,
        col = col,
        legend.text = legend_text,
        args.legend = args_legend,
        xlab = xlab,
        ylab = ylab
    )
    dots[["height"]] <- NULL
    barplot_args[names(dots)] <- dots

    if (isTRUE(plot))
        do.call(graphics::barplot, barplot_args)
    arch_names <- archetype_names %||% rownames(A) %||% seq_len(nrow(A))
    feat_names <- colnames(A) %||% seq_len(ncol(A))
    out <- data.frame(
        archetype = rep(arch_names, times = ncol(A)),
        feature   = rep(feat_names, each  = nrow(A)),
        value     = as.vector(A),
        stringsAsFactors = FALSE
    )
    invisible(out)
}

#' Coordinate Plot For Archetypes
#'
#' Plots archetype coordinates, optionally overlaid on observation data. When
#' `data` is supplied and `projection = "pca"`, observations and coordinates are
#' projected together onto the first two principal components before plotting.
#'
#' @param coordinates Numeric matrix of archetype coordinates. This is the only
#'   required argument.
#' @param data Optional numeric matrix of observations to draw behind the
#'   archetypes.
#' @param projection Projection to use. `"pca"` is only valid when `data` is
#'   supplied.
#' @param archetype_names Optional labels to draw next to archetype points.
#' @param show_anames Logical. Should archetype labels be drawn?
#' @param args.data.scatter Named list of graphical arguments for observation
#'   points.
#' @param plot Logical. Should the plot be drawn?
#' @param ... Graphical arguments for archetype points and paths. General plot
#'   window arguments such as `main`, `xlab`, `ylab`, `xlim`, `ylim`, and `asp`
#'   are also honored when drawing the plot.
#'
#' @return Invisibly returns a data frame in long format with one row per
#'   point. Columns are the coordinate dimensions, `name` (character label),
#'   and `archetype` (logical). Data rows come first, archetype rows last.
#'
#' @export
plot_archetypes_coordinates <- function(coordinates,
                                        data = NULL,
                                        projection = c("none", "pca"),
                                        archetype_names = NULL,
                                        show_anames = TRUE,
                                        args.data.scatter = list(),
                                        plot = TRUE,
                                        ...) {
    if (!is.list(args.data.scatter))
        stop("`args.data.scatter` must be a list", call. = FALSE)
    projection <- match.arg(projection)
    if (is.null(data) && !identical(projection, "none"))
        stop("`projection` can only be 'none' when `data` is NULL.", call. = FALSE)

    A <- as.matrix(coordinates)
    if (!is.numeric(A))
        stop("`coordinates` must be numeric.", call. = FALSE)
    if (nrow(A) == 0L || ncol(A) == 0L)
        stop("`coordinates` must have at least one row and one column.", call. = FALSE)
    X <- NULL
    if (!is.null(data)) {
        X <- as.matrix(data)
        if (!is.numeric(X))
            stop("`data` must be numeric.", call. = FALSE)
        if (ncol(X) != ncol(A)) {
            fmt <- "`data` has %d columns but `coordinates` has %d columns"
            stop(sprintf(fmt, ncol(X), ncol(A)), call. = FALSE)
        }
        if (nrow(X) == 0L)
            stop("`data` must contain at least one row.", call. = FALSE)
    }

    if (is.null(colnames(A))) {
        if (!is.null(X) && !is.null(colnames(X))) {
            colnames(A) <- colnames(X)
        } else {
            colnames(A) <- paste0("V", seq_len(ncol(A)))
        }
    }
    if (!is.null(X) && is.null(colnames(X)))
        colnames(X) <- colnames(A)

    pc <- NULL
    if (identical(projection, "pca")) {
        pc <- stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = 2L)
        X  <- pc[["x"]]
        A  <- predict(pc, newdata = A)
        colnames(X) <- colnames(A) <- c("PC1", "PC2")
    }

    if (ncol(A) < 2L)
        stop("Coordinate plots require at least two dimensions.", call. = FALSE)

    dots <- list(...)
    coord_args <- .aa_coordinate_args(dots, args.data.scatter)  # internal helper; uses %|p|%
    canvas_args <- coord_args[["canvas"]]
    archetype_args <- coord_args[["archetype"]]
    data_args <- coord_args[["data"]]

    if (isTRUE(plot)) {
        if (ncol(A) == 2L) {
            .aa_plot_coordinates_2d(
                coordinates = A,
                data = X,
                archetype_names = archetype_names,
                show_anames = show_anames,
                canvas_args = canvas_args,
                archetype_args = archetype_args,
                data_args = data_args
            )
        } else {
            .aa_plot_coordinates_pairs(
                coordinates = A,
                data = X,
                archetype_names = archetype_names,
                show_anames = show_anames,
                canvas_args = canvas_args,
                archetype_args = archetype_args,
                data_args = data_args
            )
        }
    }

    df_arch <- data.frame(
        name = archetype_names %||% rownames(A) %||% seq_len(nrow(A)),
        archetype = TRUE,
        stringsAsFactors = FALSE
    )
    df_arch <- cbind(as.data.frame(A), df_arch)

    if (is.null(X))
        return(invisible(df_arch))

    df_X <- data.frame(
        name = rownames(X) %||% seq_len(nrow(X)),
        archetype = FALSE,
        stringsAsFactors = FALSE
    )
    df_X <- cbind(as.data.frame(X), df_X)
    invisible(rbind(df_X, df_arch))
}

.aa_composition_dist <- function(data, distance) {
    if (inherits(distance, "dist"))
        return(distance)
    if (is.function(distance)) {
        d <- distance(data)
        if (inherits(d, "dist"))
            return(d)
        return(stats::as.dist(d))
    }

    stopifnot("Clustering distance must be a non-empty string, function, or dist object" = is_non_empty_string(distance))
    distance <- ifelse(distance == "correlation", "pearson", distance)
    if (distance %in% c("pearson", "spearman", "kendall"))
        return(stats::as.dist(1 - stats::cor(t(data), method = distance)))
    stats::dist(data, method = distance)
}

.aa_component_order <- function(components, mode) {
    if (mode == "PC1")
        return(order(components[, 1L]))
    if (ncol(components) < 2L)
        return(NULL)
    order(atan2(components[, 2L], components[, 1L]))
}

.aa_composition_hclust <- function(value, data, margin, distance, linkage) {
    cluster_arg <- if (margin == "rows") "cluster_rows" else "cluster_cols"
    if (inherits(value, "hclust"))
        return(value)
    if (is_single_string(value)) {
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
        order <- .aa_component_order(components, mode)
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
        return(stats::hclust(.aa_composition_dist(transformed, distance), method = linkage))
    }
    if (ncol(data) < 2L)
        return(NULL)
    stats::hclust(.aa_composition_dist(t(data), distance), method = linkage)
}

.aa_coordinate_args <- function(dots, args.data.scatter) {
    canvas_names <- c(
        "main", "sub", "xlab", "ylab", "xlim", "ylim", "asp", "frame.plot",
        "axes", "ann", "type", "las", "cex.axis", "cex.lab", "cex.main",
        "col.axis", "col.lab", "col.main", "font.axis", "font.lab",
        "font.main", "labels"
    )
    canvas <- dots[intersect(names(dots), canvas_names)]
    archetype_dots <- dots[setdiff(names(dots), canvas_names)]
    list(
        canvas = canvas,
        archetype = list(col = "red", pch = 16, cex = 1.3, lwd = 1) %|p|% archetype_dots,
        data = list(col = "black", pch = 1, cex = 1) %|p|% args.data.scatter
    )
}

.aa_plot_coordinates_2d <- function(coordinates,
                                    data,
                                    archetype_names,
                                    show_anames,
                                    canvas_args,
                                    archetype_args,
                                    data_args) {
    A <- coordinates
    X <- data
    canvas_args[["xlab"]] <- canvas_args[["xlab"]] %||% colnames(A)[1L]
    canvas_args[["ylab"]] <- canvas_args[["ylab"]] %||% colnames(A)[2L]
    canvas_args <- list(asp = 1) %|p|% canvas_args

    if (is.null(X)) {
        plot_args <- list(x = A[, 1L], y = A[, 2L], type = "n") %|p|% canvas_args
        do.call(graphics::plot, plot_args)
    } else {
        plot_args <- c(list(x = X[, 1L], y = X[, 2L]), data_args) %|p|% canvas_args
        do.call(graphics::plot, plot_args)
    }

    A_closed <- A[c(seq_len(nrow(A)), 1L), , drop = FALSE]
    graphics::lines(
        A_closed[, 1L],
        A_closed[, 2L],
        col = archetype_args[["col"]],
        lwd = archetype_args[["lwd"]],
        lty = archetype_args[["lty"]] %||% 1
    )
    graphics::points(
        A[, 1L],
        A[, 2L],
        col = archetype_args[["col"]],
        pch = archetype_args[["pch"]],
        cex = archetype_args[["cex"]]
    )
    if (isTRUE(show_anames) && !is.null(archetype_names)) {
        graphics::text(
            x = A[, 1L],
            y = A[, 2L],
            labels = archetype_names,
            pos = 4,
            cex = 0.8,
            col = archetype_args[["col"]],
            xpd = NA
        )
    }
}

.aa_plot_coordinates_pairs <- function(coordinates,
                                       data,
                                       archetype_names,
                                       show_anames,
                                       canvas_args,
                                       archetype_args,
                                       data_args) {
    combined <- rbind(data, coordinates)
    n_data <- nrow(data) %||% 0L
    n_total <- nrow(combined)
    panel <- function(x, y, ...) {
        if (n_data > 0L) {
            data_ix <- seq_len(n_data)
            graphics::points(
                x[data_ix],
                y[data_ix],
                col = data_args[["col"]],
                pch = data_args[["pch"]],
                cex = data_args[["cex"]]
            )
        }
        arch_ix <- (n_data + 1L):n_total
        arch_closed <- c(arch_ix, arch_ix[1L])
        graphics::lines(
            x[arch_closed],
            y[arch_closed],
            col = archetype_args[["col"]],
            lwd = archetype_args[["lwd"]],
            lty = archetype_args[["lty"]] %||% 2
        )
        graphics::points(
            x[arch_ix],
            y[arch_ix],
            col = archetype_args[["col"]],
            pch = archetype_args[["pch"]],
            cex = archetype_args[["cex"]]
        )
        if (isTRUE(show_anames) && !is.null(archetype_names)) {
            graphics::text(
                x = x[arch_ix],
                y = y[arch_ix],
                labels = archetype_names,
                pos = 4,
                cex = 0.7,
                col = archetype_args[["col"]],
                xpd = NA
            )
        }
    }
    pair_args <- list(x = combined, panel = panel, lower.panel = panel, upper.panel = panel)
    do.call(graphics::pairs, pair_args %|p|% canvas_args)
}
