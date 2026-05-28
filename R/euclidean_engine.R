# euclidean_engine.R: Shared Euclidean fitting lifecycle helpers

.aa_check_missing_route <- function(ctx) {
    if (!(ctx[["method"]] %in% c("pgd", "fw")) && ctx[["missing"]]) {
        stop("`missing = TRUE` is only supported for `method = 'pgd'` or `method = 'fw'`.", call. = FALSE)
    }
    if (ctx[["missing"]] && !identical(ctx[["robust"]], FALSE)) {
        stop("`robust` is not supported with `missing = TRUE`.", call. = FALSE)
    }
    if (ctx[["missing"]] && !is.null(ctx[["weights"]])) {
        stop("`weights` are not supported with `missing = TRUE`.", call. = FALSE)
    }
    if (ctx[["missing"]] && is_matrix(ctx[["scale"]])) {
        stop("matrix `scale` is not supported with `missing = TRUE`.", call. = FALSE)
    }
    invisible(TRUE)
}

.aa_euclidean_check <- function(ctx) {
    .aa_check_missing_route(ctx)
    .aa_check_fit_controls(ctx)
    .aa_check_scale(ctx[["scale"]], ncol(ctx[["x"]]))
}

.aa_euclidean_preprocess <- function(ctx, bigM = 0) {
    data <- ctx[["x"]]
    sd_threshold <- ctx[["sd_threshold"]]
    weights <- ctx[["weights"]]
    scale <- ctx[["scale"]]

    if (ctx[["verbose"]]) message("Preprocessing data...")

    if (ctx[["missing"]]) {
        return(.aa_preprocess_missing(data, sd_threshold, ctx[["verbose"]], scale = scale))
    }

    if (inherits(data, "sparseMatrix")) {
        data <- Matrix::drop0(data)
    }

    original_center <- colMeans(data)
    names(original_center) <- colnames(data)
    scale_mode <- if (identical(scale, FALSE)) {
        "none"
    } else if (isTRUE(scale) || (is.numeric(scale) && is.null(dim(scale)))) {
        "vector"
    } else {
        "matrix"
    }

    X <- if (inherits(data, "sparseMatrix")) {
        data
    } else {
        as.matrix(data)
    }
    if (isTRUE(scale)) {
        n <- nrow(data)
        x_mean <- colMeans(data)
        x2 <- colSums(data * data)
        x_var <- pmax((x2 - n * x_mean * x_mean) / max(n - 1L, 1L), 0)
        scale <- sqrt(x_var)
        names(scale) <- colnames(data)

        if (!inherits(data, "sparseMatrix")) {
            X <- sweep(X, 2L, x_mean, "-")
            attr(X, "scaled:center") <- x_mean
        }
        attr(X, "scaled:scale") <- scale
    }

    X <- .aa_filter_low_variance(X, sd_threshold)
    mask <- attr(X, "mask")

    if (identical(scale_mode, "matrix")) {
        retained_names <- names(original_center)
        if (!is.null(mask)) {
            retained_names <- retained_names[mask]
            scale <- scale[mask, mask, drop = FALSE]
        }
        scale_factor <- t(chol(as.matrix(scale))) # lower triangular
        X <- as.matrix(X) %*% scale_factor
        colnames(X) <- retained_names
        attr(X, "mask") <- mask
        # Compressed storage
        attr(X, "scale:factor") <- .aa_pack_lower_tri(scale_factor)
    } else if (identical(scale_mode, "vector")) {
        if (!is.null(mask)) {
            scale <- scale[mask]
        }
        scale_factor <- ifelse(scale > 0, scale, 1)
        x_attrs <- attributes(X)
        X <- if (inherits(X, "sparseMatrix")) {
            Matrix::colScale(X, 1 / scale_factor)
        } else {
            sweep(X, 2L, scale_factor, "/")
        }
        attributes(X) <- utils::modifyList(attributes(X), x_attrs)
        attr(X, "scale:factor") <- scale_factor
    }
    attr(X, "scale:mode") <- scale_mode
    attr(X, "restore:center") <- original_center
    N <- nrow(X)
    if (is.null(bigM)) {
        bigM <- .aa_auto_bigM(X)
    }

    if (!is.null(weights)) {
        weights <- .aa_check_sample_weights(weights, N)
        weights <- weights / mean(weights)
        x_attrs <- attributes(X)
        X <- X * sqrt(weights) # sqrt will be undone during square loss computation
        attributes(X) <- utils::modifyList(attributes(X), x_attrs)
        attr(X, "weights") <- weights
    }

    if (bigM > 0) {
        x_attrs <- attributes(X)
        bigM_col <- matrix(bigM, nrow = N, ncol = 1L, dimnames = list(rownames(X), "bigM"))
        if (inherits(X, "sparseMatrix")) {
            bigM_col <- as(bigM_col, "sparseMatrix")
        }
        X <- cbind(bigM_col, X)
        attr(X, "scaled:center") <- x_attrs[["scaled:center"]]
        attr(X, "scaled:scale") <- x_attrs[["scaled:scale"]]
        attr(X, "scale:mode") <- x_attrs[["scale:mode"]]
        attr(X, "scale:factor") <- x_attrs[["scale:factor"]]
        attr(X, "restore:center") <- x_attrs[["restore:center"]]
        attr(X, "mask") <- x_attrs[["mask"]]
        attr(X, "bigM") <- 1L
        attr(X, "bigM.value") <- bigM
    }

    list(X = X)
}

.aa_euclidean_edge_case <- function(ctx, prep) {
    out <- .aa_checks_edge_cases(
        prep[["X"]],
        ctx[["K"]],
        ctx[["verbose"]],
        M = prep[["M"]]
    )
    if (is.null(out)) {
        return(NULL)
    }
    .aa_euclidean_output(
        ctx,
        prep,
        list(
            A0 = out[["init"]],
            A = coordinates(out),
            B = out[["coefficients"]],
            S = out[["compositions"]],
            delta = out[["slack"]],
            i = nrow(out[["loss"]]) - 1L,
            loss = out[["loss"]],
            converged = out[["converged"]]
        ),
        fit_info = list(method = ctx[["method"]])
    )
}

.aa_euclidean_init <- function(ctx, prep, delta = 0) {
    X <- prep[["X"]]
    init <- ctx[["init"]]
    init_args <- ctx[["init_args"]]
    if (ctx[["verbose"]]) message("Initializing archetypes...")
    L <- ctx[["max_iter"]] + 1L

    if (is_tabular(init)) {
        init <- .aa_preprocess_init(init, X)
        if (length(init_args) > 0L) {
            warning("`init_args` are ignored when `init` is a matrix", call. = FALSE)
            init_args <- list()
        }
        return(.aa_matrix_init(X, ctx[["K"]], init, ctx[["eps"]], L, delta))
    }

    if (is_non_empty_string(init)) {
        init_args <- c(list(method = init), init_args)
        init <- aa_init
    } else if (!is.function(init)) {
        stop("`init` must be a function, a single non-empty string, or archetypes coordinate matrix")
    }

    init_vars <- do.call(init, args = c(list(X = X, K = ctx[["K"]]), init_args))
    rownames(init_vars[["A"]]) <- .aa_init_names(init_vars[["A"]])
    rownames(init_vars[["B"]]) <- rownames(init_vars[["A"]])
    init_vars[["S"]] <- .aa_init_S(X, init_vars[["A"]], eps = ctx[["eps"]])
    init_vars[["loss"]] <- list(loss = rep(NA_real_, L), r2 = rep(NA_real_, L))
    init_vars
}

.aa_euclidean_output <- function(ctx, prep, fit, fit_info = list()) {
    X <- if (!is.null(prep[["output_X"]])) prep[["output_X"]] else prep[["X"]]
    feature_map <- .aa_euclidean_feature_map(X)
    feature_map[["eps"]] <- ctx[["eps"]]

    j <- fit[["i"]] + 1L
    loss <- as.data.frame(fit[["loss"]])[seq_len(j), , drop = FALSE]
    rownames(loss) <- NULL

    A_feature <- .aa_drop_bigM_column(fit[["A"]], attr(X, "bigM"))
    A <- fit[["A"]]
    A0 <- fit[["A0"]]
    init_arg <- ctx[["init"]]
    archetype_names <- NULL
    if (is_tabular(init_arg)) {
        archetype_names <- rownames(init_arg)
    } else if (is.function(init_arg)) {
        archetype_names <- if (!is.null(A0)) rownames(A0) else rownames(A)
    }
    A <- .aa_feature_map_inverse(feature_map, A_feature)
    if (!is.null(A0)) {
        A0 <- .aa_feature_map_inverse(feature_map, .aa_drop_bigM_column(A0, attr(X, "bigM")))
    }

    family <- prep[["family"]] %||% ctx[["family"]] %||% "gaussian"

    rownames(A) <- rownames(fit[["B"]]) <- colnames(fit[["S"]]) <- archetype_names
    rownames(A_feature) <- archetype_names
    if (!is.null(A0)) {
        rownames(A0) <- archetype_names
    }
    colnames(fit[["B"]]) <- rownames(fit[["S"]]) <- rownames(X)

    if (!fit[["converged"]]) {
        fmt <- "Algorithm did not converge after %d iterations"
        warning(sprintf(fmt, fit[["i"]]), call. = FALSE)
    }

    if (ctx[["verbose"]]) {
        fmt <- ifelse(fit[["converged"]],
            "Converged after %d iterations:",
            "Final iteration %d:"
        )
        fmt <- paste(fmt, "loss = %.4g, R2 = %.3f")
        message(sprintf(fmt, fit[["i"]], loss[j, "loss"], loss[j, "r2"]))
    }

    init <- ctx[["init"]]
    if (!is.null(init) && !is.character(init)) {
        init <- if (is.function(init)) "function" else "matrix"
    }
    fit_info <- fit_info %|p|% list(
        family = prep[["family"]] %||% ctx[["family"]] %||% "gaussian",
        robust = !identical(ctx[["robust"]], FALSE),
        robust_psi = .aa_robust_label(ctx[["robust"]]) %||% NA_character_,
        robust_args = I(list(ctx[["robust_args"]])),
        missing = isTRUE(ctx[["missing"]]),
        delta = fit[["delta"]] %||% 0,
        init = init,
        scaling = .aa_scale_label(ctx[["scale"]]),
        sample_weights = !is.null(ctx[["weights"]])
    )

    archetypes(
        call         = ctx[["call"]],
        data         = ctx[["x"]],
        weights      = if (!is.null(fit[["row_weights"]])) fit[["row_weights"]] else attr(prep[["X"]], "weights"),
        init         = A0,
        coefficients = fit[["B"]],
        compositions = fit[["S"]],
        slack        = fit[["delta"]] %||% 0,
        loss         = loss,
        converged    = fit[["converged"]],
        family       = family,
        fit_info     = fit_info,
        feature_map  = feature_map,
        A            = A_feature
    )
}

.aa_scale_label <- function(scale) {
    if (identical(scale, FALSE)) {
        return("none")
    }
    if (isTRUE(scale)) {
        return("z-score")
    }
    if (is.numeric(scale) && is.null(dim(scale))) {
        return("custom")
    }
    if (is_matrix(scale)) {
        return("metric")
    }
    "custom"
}
