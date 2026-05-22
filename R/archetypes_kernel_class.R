# Kernel Archetypes Class -----------------------------------------------------

#' Internal kernel archetype analysis constructor
#'
#' Builds a `kernel_archetypes` object from fitted kernel AA matrices.
#' `coefficients` must be K x N, `compositions` N x K, and `gram` N x N.
#' The base `archetypes()` constructor enforces stochasticity and shared object
#' invariants; this constructor adds Gram matrix, kernel, and input-space proxy metadata.
#'
#' @param coefficients Numeric matrix (`K x N`) giving each Hilbert-space
#'   archetype as a weighted combination of training samples.
#' @param compositions Numeric matrix (`N x K`) giving each sample as a weighted
#'   combination of kernel archetypes.
#' @param gram Training Gram matrix.
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
#' @param fit_info Internal fit metadata.
#' @param feature_map Internal prediction feature-map metadata.
#'
#' @noRd
.aa_new_kernel_archetypes <- function(coefficients,
                                      compositions,
                                      gram,
                                      slack = 0,
                                      loss = NULL,
                                      converged = TRUE,
                                      call = NULL,
                                      data = NULL,
                                      init = NULL,
                                      kernel = NULL,
                                      kernel_args = list(),
                                      weights = NULL,
                                      fit_info = list(),
                                      feature_map = .aa_kernel_feature_map()) {
    call <- call %||% match.call()
    loss <- loss %||% data.frame(loss = NA_real_, r2 = NA_real_)
    if (!identical(dim(gram), c(ncol(coefficients), ncol(coefficients)))) {
        stop(
            paste(
                "Gram matrix dimensions must match the number of samples",
                "encoded by `coefficients`."
            ),
            call. = FALSE
        )
    }
    out <- archetypes(
        A = NULL,
        coefficients = coefficients,
        compositions = compositions,
        slack = slack,
        loss = loss,
        converged = converged,
        call = call,
        data = data,
        init = NULL, # kernel init is a coefficient matrix; add after
        family = "gaussian",
        weights = weights,
        fit_info = fit_info,
        feature_map = feature_map
    )
    if (!is.null(init) && !is.null(rownames(init))) {
        archetype_names <- rownames(out[["coefficients"]])
        if (!is.null(archetype_names) && setequal(rownames(init), archetype_names)) {
            init <- init[archetype_names, , drop = FALSE]
        }
    }
    out[["init"]] <- init
    out[["gram"]] <- gram
    out[["kernel"]] <- kernel
    out[["kernel_args"]] <- kernel_args
    class(out) <- c("kernel_archetypes", "archetypes")
    out
}

#' Predict method for kernel archetypes
#'
#' Predict compositions for kernel archetypes
#'
#' Projects new samples onto the fitted kernel archetype simplex.
#' `type = "compositions"` returns the fitted simplex weights.
#' `type = "reconstruction"` errors because inverse reconstruction from the
#' implicit feature space is undefined for kernel archetypes.
#'
#' @param object An object of class `kernel_archetypes`.
#' @param newdata New data to project, or an `N_new x N_train` cross-kernel
#'   matrix when `object` was fit with `kernel = "precomputed"`.
#' @param type Prediction output type. Only `"compositions"` is defined.
#' @param max_iter,eps,step_size,max_iter_optimizer,step_shrinkage Optimization
#'   controls for composition fitting.
#' @param ... Ignored.
#'
#' @return A numeric `N_new x K` row-stochastic composition matrix.
#'
#' @exportS3Method
predict.kernel_archetypes <- function(object,
                                      newdata,
                                      type = c("reconstruction", "compositions"),
                                      max_iter = 100L,
                                      eps = 1e-8,
                                      step_size = 1.0,
                                      max_iter_optimizer = 10L,
                                      step_shrinkage = 0.5,
                                      ...) {
    type <- match.arg(type)
    if (identical(type, "reconstruction")) {
        stop(
            "`type = \"reconstruction\"` is undefined for kernel_archetypes objects.",
            call. = FALSE
        )
    }

    C <- .aa_kernel_cross_gram(object, newdata)
    .aa_kernel_predict_S(
        C = C,
        G = object[["gram"]],
        B = coefficients(object),
        max_iter = max_iter,
        eps = eps,
        step_size = step_size,
        max_iter_optimizer = max_iter_optimizer,
        step_shrinkage = step_shrinkage
    )
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

#' Residuals for kernel archetypes objects
#'
#' Residuals are per-sample squared distances in the implicitly defined feature space.
#' The method is provided for completness and to support residual-based diagnostics as
#' the natural `fitted()` method is not defined for nonlinear kernel archetypes.
#'
#' @param object An object of class `kernel_archetypes`.
#' @param ... Ignored.
#'
#' @return A named numeric vector of non-negative squared residual norms, one
#'   per sample.
#'
#' @seealso [residuals.archetypes()] for Euclidean residuals on standard fits.
#'
#' @exportS3Method
residuals.kernel_archetypes <- function(object, ...) {
    G <- object[["gram"]]
    H <- object[["coefficients"]]
    S <- compositions(object)
    AAt <- tcrossprod(H %*% G, H)
    XAt <- tcrossprod(G, H)
    out <- diag(G) - 2 * rowSums(S * XAt) + rowSums(S * (S %*% AAt))
    names(out) <- rownames(S)
    pmax(out, 0)
}

#' @exportS3Method
print.kernel_archetypes <- function(x, ...) {
    call_str <- paste(deparse(x[["call"]]), sep = "\n", collapse = "\n")
    cat("\nCall:\n", call_str, "\n\n", sep = "")
    cat(.aa_print_title(x), ":\n", sep = "")
    cat("Number of Archetypes:", nrow(x[["coefficients"]]), "\n")
    cat("Number of Samples:", ncol(x[["coefficients"]]), "\n")
    if (!is.null(x[["data"]])) {
        cat("Input-space coordinate proxy: available\n")
    }
    loss <- x[["loss"]]
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
plot.kernel_archetypes <- function(x,
                                   what = c("compositions", "loss", "coordinates", "profiles"),
                                   subset = NULL,
                                   plot = TRUE,
                                   ...) {
    what <- .aa_plot_normalize_what(what)

    if (what == "profiles") {
        stop(
            "Profile plots are not defined for kernel archetypes because ",
            "their natural coordinates live in implicit Hilbert space.",
            call. = FALSE
        )
    }

    dots <- list(...)
    if (what == "coordinates" && identical(dots[["projection"]], "pca")) {
        kpca <- .aa_kernel_kpca(x[["gram"]], x[["coefficients"]])
        dots[["projection"]] <- NULL
        args <- list(
            coordinates     = kpca[["archetypes"]],
            data            = .aa_plot_subset_rows(kpca[["data"]], subset),
            archetype_names = anames(x),
            projection      = "none",
            plot            = plot
        ) %|p|% dots
        return(do.call(plot_archetypes_coordinates, args))
    }

    NextMethod()
}

#' @exportS3Method
AIC.kernel_archetypes <- function(object, ...) {
    stop("AIC is not defined for kernel archetypes.", call. = FALSE)
}
