# feature_map.R: Internal feature-space transforms for prediction

.aa_feature_map_transform <- function(map, newdata, ...)
    UseMethod(".aa_feature_map_transform")

.aa_feature_map_inverse <- function(map, z, ...)
    UseMethod(".aa_feature_map_inverse")

.aa_identity_feature_map <- function(coordinates) {
    A <- as.matrix(coordinates)
    p <- ncol(A)
    original_names <- colnames(A) %||% paste0("V", seq_len(p))
    center <- rep(0, p)
    names(center) <- original_names
    structure(
        list(
            type = "euclidean",
            scale_mode = "none",
            original_names = original_names,
            retained = rep(TRUE, p),
            center = center,
            restore_center = center,
            scale_factor = NULL
        ),
        class = c("euclidean_feature_map", "aa_feature_map")
    )
}

.aa_unsupported_feature_map <- function(method) {
    structure(
        list(type = "unsupported", method = method),
        class = c("unsupported_feature_map", "aa_feature_map")
    )
}

.aa_feature_map_transform.unsupported_feature_map <- function(map, newdata, ...) {
    msg <- sprintf("Feature-map prediction is not defined for %s fits.",
                   map[["method"]])
    stop(msg, call. = FALSE)
}

.aa_feature_map_inverse.unsupported_feature_map <- function(map, z, ...)
    NULL

.aa_kernel_feature_map <- function() {
    structure(
        list(type = "kernel"),
        class = c("kernel_feature_map", "aa_feature_map")
    )
}

.aa_feature_map_transform.kernel_feature_map <- function(map, newdata, ...)
    stop("Feature-map prediction is not defined for kernel fits.", call. = FALSE)

.aa_feature_map_inverse.kernel_feature_map <- function(map, z, object = NULL, ...) {
    if (is.null(object) || is.null(object[["data"]]))
        return(NULL)

    X <- object[["data"]]
    if (inherits(X, "fd"))
        X <- .aa_fd_to_matrix(X)

    B <- coefficients(object)
    A <- B %*% X
    rownames(A) <- rownames(B)
    colnames(A) <- colnames(X)
    A
}

.aa_drop_bigM_column <- function(mat, bigM = NULL) {
    if (is.null(bigM))
        return(as.matrix(mat))
    as.matrix(mat)[, -bigM, drop = FALSE]
}

.aa_euclidean_feature_map <- function(X) {
    iM <- attr(X, "bigM")
    X_no_bigM <- if (is.null(iM)) X else X[, -iM, drop = FALSE]
    x_names <- colnames(X_no_bigM)

    restore_center <- attr(X, "restore:center")
    if (is.null(restore_center)) {
        restore_center <- rep(0, ncol(X_no_bigM))
        names(restore_center) <- x_names %||% paste0("V", seq_along(restore_center))
    }
    original_names <- names(restore_center) %||% colnames(X_no_bigM)
    if (is.null(original_names))
        original_names <- paste0("V", seq_along(restore_center))
    names(restore_center) <- original_names

    retained <- attr(X, "mask")
    if (is.null(retained))
        retained <- rep(TRUE, length(original_names))
    names(retained) <- original_names

    center <- attr(X, "scaled:center")
    if (!is.null(center))
        names(center) <- names(center) %||% original_names

    structure(
        list(
            type = "euclidean",
            scale_mode = attr(X, "scale:mode") %||% "none",
            original_names = original_names,
            retained = retained,
            center = center,
            restore_center = restore_center,
            scale_factor = attr(X, "scale:factor")
        ),
        class = c("euclidean_feature_map", "aa_feature_map")
    )
}

.aa_feature_map_select_columns <- function(map, newdata) {
    X <- if (inherits(newdata, "data.table")) {
        newdata[, map[["original_names"]], with = FALSE]
    } else if (!is.null(colnames(newdata))) {
        missing <- setdiff(map[["original_names"]], colnames(newdata))
        if (length(missing) > 0L)
            stop(sprintf("`newdata` is missing required columns: %s",
                         paste(missing, collapse = ", ")), call. = FALSE)
        newdata[, map[["original_names"]], drop = FALSE]
    } else {
        newdata
    }
    X <- if (inherits(X, "sparseMatrix")) X else as.matrix(X)
    if (ncol(X) != length(map[["original_names"]])) {
        fmt <- "`newdata` has %d columns but the feature map expects %d columns"
        stop(sprintf(fmt, ncol(X), length(map[["original_names"]])), call. = FALSE)
    }
    colnames(X) <- map[["original_names"]]
    X
}

.aa_feature_map_transform.euclidean_feature_map <- function(map, newdata, ...) {
    X <- .aa_feature_map_select_columns(map, newdata)
    retained <- map[["retained"]]
    X <- X[, retained, drop = FALSE]

    scale_mode <- map[["scale_mode"]]
    if (identical(scale_mode, "vector")) {
        center <- map[["center"]]
        if (!is.null(center))
            X <- sweep(as.matrix(X), 2L, center[retained], "-")

        scale_factor <- map[["scale_factor"]] %||% rep(1, ncol(X))
        X <- if (inherits(X, "sparseMatrix")) {
            Matrix::colScale(X, 1 / scale_factor)
        } else {
            sweep(X, 2L, scale_factor, "/")
        }
    } else if (identical(scale_mode, "matrix")) {
        X <- as.matrix(X) %*% .aa_unpack_lower_tri(map[["scale_factor"]])
    } else {
        X <- if (inherits(X, "sparseMatrix")) X else as.matrix(X)
    }

    colnames(X) <- map[["original_names"]][retained]
    X
}

.aa_feature_map_inverse.euclidean_feature_map <- function(map, z, ...) {
    mat <- as.matrix(z)
    mat_names <- rownames(mat)
    retained <- map[["retained"]]
    if (ncol(mat) != sum(retained)) {
        fmt <- "transformed data has %d columns but the feature map expects %d retained columns"
        stop(sprintf(fmt, ncol(mat), sum(retained)), call. = FALSE)
    }

    scale_mode <- map[["scale_mode"]]
    if (identical(scale_mode, "matrix")) {
        L <- .aa_unpack_lower_tri(map[["scale_factor"]]) # we have stored t(chol(scale))
        mat <- t(backsolve(t(L), t(mat)))
    } else if (identical(scale_mode, "vector")) {
        scale_factor <- map[["scale_factor"]] %||% rep(1, ncol(mat))
        mat <- sweep(mat, 2L, scale_factor, "*")
        center <- map[["center"]]
        if (is.null(center)) {
            center <- rep(0, length(map[["restore_center"]]))
            names(center) <- map[["original_names"]]
        }
        mat <- sweep(mat, 2L, center[retained], "+")
    }
    rownames(mat) <- mat_names
    colnames(mat) <- map[["original_names"]][retained]

    x_mean <- map[["restore_center"]]
    out <- matrix(
        x_mean,
        nrow = nrow(mat),
        ncol = length(x_mean),
        byrow = TRUE,
        dimnames = list(rownames(mat), map[["original_names"]])
    )

    out[, retained] <- mat
    out
}
