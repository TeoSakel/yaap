# Information criteria for archetypes objects ---------------------------------

#' Information Criteria for Archetype Fits
#'
#' Computes AIC or BIC for fitted archetype models
#' with selectable degrees-of-freedom estimates.
#'
#' @param object An object of class `archetypes`.
#' @param type Criterion to compute: `"aic"` or `"bic"`.
#' @param df_method Degrees-of-freedom estimate. "suleman" uses a
#'   Suleman-style efficiency-adjusted parameter count; "covariance" uses a
#'   fitted-value sensitivity estimate; "parametric" uses the full
#'   archetype parameter count; "active" uses an active-set simplex-face
#'   count.
#' @param n_eff Effective sample size used for the BIC penalty multiplier:
#'   "samples" uses rows of the training matrix; "entries" uses observed
#'   matrix entries. The Gaussian RSS term is always normalized by observed
#'   matrix entries.
#' @param support_tol Non-negative threshold used by df_method = active to
#'   decide whether simplex weights belong to a support set.
#' @param ... Additional arguments passed to `aa_ic()`.
#'
#' @return A numeric scalar.
#'
#' @details
#' `df_method = "suleman"` corrects the full parametric degree of freedom count
#' based on how well the fitted covariance structure matches the data covariance,
#' as proposed by Suleman (2017). Unlike the normalized criterion in Suleman (2017),
#' `aa_ic()` applies this correction on the same Gaussian RSS scale used by the other
#' `df_method` choices for consistency.
#'
#' For `df_method = "covariance"`, the fitted archetype weights are held fixed and
#' the degrees of freedom are estimated from how strongly the fitted values
#' respond to the data. In linear-model terms, this is the trace of the
#' corresponding hat matrix. Here the fitted-value map is treated as
#' `X_hat = H X`, with `H = compositions(object) %*% coefficients(object)`, so
#' the degrees of freedom are `M * sum(diag(H))`, where `M` is the number of
#' features. This is Stein's frozen-smoother covariance formula applied to the
#' fixed archetype-weight map.
#'
#' For `df_method = "parametric"`, the structural degrees of freedom are the
#' full ambient simplex parameter count: `N * (K - 1)` free composition
#' weights and `K * (N - 1)` free coefficient weights.
#'
#' For `df_method = "active"`, unlike "parametric", only the active simplex-face
#' dimensions (ie with non-zero coefficients) are counted towards the degrees of
#' freedom of the compositions and coefficients simplices. The coefficient-face
#' dimensions are further capped by the local affine rank of the corresponding
#' rows of `X` because the archetype vertices must live in the affine span of
#' the data. This estimate is more accurate for methods that produce sparse
#' simplex weights through regularization or constraints.
#'
#' For criteria with a free scalar residual parameter, such as Gaussian fits,
#' one additional degree of freedom is added outside the structural `df_method`
#' calculation, assuming homoscedastic errors. Gaussian AIC and BIC use the
#' profile log-RSS scale `nelem * log(rss / nelem)`, where `nelem` is the number
#' of observed matrix entries. Missing data are not counted towards `nelem`.
#' `n_eff` affects the BIC penalty multiplier but not this RSS scale. For
#' non-Gaussian PAA families, AIC and BIC use twice the final optimized
#' objective (deviance) plus the selected penalty.
#'
#' @references
#' A. Suleman, "Validation of archetypal analysis"
#' 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE),
#' Naples, Italy, 2017, pp. 1-6, \doi{10.1109/FUZZ-IEEE.2017.8015385}
#'
#' Charles M. Stein (1981). Estimation of the Mean of a Multivariate Normal
#' Distribution. *The Annals of Statistics*, 9(6), 1135-1151.
#'
#' @export
aa_ic <- function(object,
                  type = c("aic", "bic"),
                  df_method = c("suleman", "covariance", "parametric", "active"),
                  n_eff = c("samples", "entries"),
                  support_tol = 1e-8) {
    if (!inherits(object, "archetypes")) {
        stop("`object` must be an `archetypes` object.", call. = FALSE)
    }
    type <- match.arg(type)
    df_method <- match.arg(df_method)
    n_eff <- match.arg(n_eff)

    X <- .aa_object_data_matrix(object, toupper(type))
    family <- object[["family"]] %||% "gaussian"
    df <- .aa_ic_total_df(object, X, df_method, family, support_tol)
    df <- .aa_ic_adjust_df(object, X, df, df_method, family, type)
    if (is.na(df)) {
        return(NA_real_)
    }
    n <- .aa_ic_n_eff(object, X, n_eff)
    penalty <- switch(type,
        aic = 2,
        bic = log(n)
    )
    if (identical(family, "gaussian")) {
        fit_term <- .aa_ic_gaussian_fit_term(object, X)
        return(fit_term + penalty * df)
    }

    if (!(family %in% c("binomial", "poisson", "multinomial"))) {
        stop(sprintf(
            "Information criteria are not defined for family '%s'.",
            family
        ), call. = FALSE)
    }
    loss <- .aa_ic_final_loss(object)
    2 * loss + penalty * df
}

#' @rdname aa_ic
#' @exportS3Method
AIC.archetypes <- function(object, ..., k = 2) {
    aa_ic(object, type = "aic", ...)
}

#' @rdname aa_ic
#' @exportS3Method
BIC.archetypes <- function(object, ...) {
    aa_ic(object, type = "bic", ...)
}

.aa_ic_adjust_df <- function(object, X, df, df_method, family, type) {
    if (!identical(df_method, "suleman")) {
        return(df)
    }
    if (!identical(family, "gaussian")) {
        stop(
            "Suleman-style information criteria are not defined for non-Gaussian archetypes objects.",
            call. = FALSE
        )
    }
    if (any(object[["slack"]] > 0)) {
        warning(paste(
            toupper(type), "computation assumes coefficients are row-stochastic;",
            "slack > 0 may violate this assumption."
        ), call. = FALSE)
    }

    N <- nrow(X)
    M <- ncol(X)
    if (N <= M) {
        warning(
            paste(
                "Suleman-style", toupper(type), "is undefined when the number of samples",
                "is not larger than the number of features; returning NA."
            ),
            call. = FALSE
        )
        return(NA_real_)
    }

    A <- .aa_input_coordinates_matrix(object)
    X_hat <- compositions(object) %*% A
    eta <- tryCatch(.aa_effic(X, X_hat), error = function(e) NA_real_)
    if (!is.finite(eta) || eta <= 0) {
        warning(
            paste(
                "Suleman-style", toupper(type), "is undefined because the efficiency term",
                "is non-positive or non-finite; returning NA."
            ),
            call. = FALSE
        )
        return(NA_real_)
    }

    df / (N * eta)
}

.aa_ic_total_df <- function(object, X, df_method, family, support_tol) {
    .aa_ic_df(object, X, df_method, support_tol) + .aa_ic_residual_parameter_df(family)
}

.aa_ic_df <- function(object, X, df_method, support_tol) {
    K <- nrow(coefficients(object))
    N <- nrow(X)
    switch(df_method,
        suleman = .aa_ic_parametric_df(N, K),
        parametric = .aa_ic_parametric_df(N, K),
        covariance = .aa_ic_covariance_df(object, X),
        active = .aa_ic_active_df(object, X, support_tol)
    )
}

.aa_ic_parametric_df <- function(N, K) {
    # N K-simplices (S) + K N-simplices (B)
    N * (K - 1L) + K * (N - 1L)
}

.aa_ic_covariance_df <- function(object, X) {
    # Stein formula: df = E[ div(X_hat) ] = E[ tr(dX_hat/dX) ]
    # X_hat = S*B*X, so H = S B and df = M * tr(H)
    S <- compositions(object)
    B <- coefficients(object)
    M <- ncol(X)
    M * sum(S * t(B))
}

.aa_ic_active_df <- function(object, X, support_tol) {
    # Active simplex-face dimensions for compositions and coefficients,
    # capped by local affine ranks because A must live in X-space.
    if (!is_non_negative(support_tol)) {
        stop("`support_tol` must be a single non-negative finite number.", call. = FALSE)
    }

    s_df <- sum(pmax(rowSums(compositions(object) > support_tol) - 1L, 0L))
    B_active <- apply(coefficients(object) > support_tol, 1L, which, simplify = FALSE)
    b_df <- sum(vapply(B_active, function(active) {
        .aa_ic_local_affine_rank(X, active)
    }, integer(1L)))
    s_df + b_df
}

.aa_ic_local_affine_rank <- function(X, active = seq_len(nrow(X))) {
    n_active <- length(active)
    if (n_active <= 1L) {
        return(0L)
    }
    centered <- scale(X[active, , drop = FALSE], center = TRUE, scale = FALSE)
    min(qr(centered)[["rank"]], n_active - 1L)
}

.aa_ic_residual_parameter_df <- function(family) {
    if (identical(family, "gaussian")) 1L else 0L
}

.aa_ic_n_eff <- function(object, X, n_eff) {
    out <- switch(n_eff,
        samples = nrow(X),
        entries = .aa_ic_n_entries(object, X)
    )
    if (!is_positive(out)) {
        stop("`n_eff` must resolve to a positive finite value.", call. = FALSE)
    }
    out
}

.aa_ic_gaussian_fit_term <- function(object, X) {
    rss <- .aa_ic_rss(object, X)
    if (!is_non_negative(rss)) {
        stop("Gaussian information criteria require finite non-negative RSS.", call. = FALSE)
    }
    nelem <- .aa_ic_n_entries(object, X)
    nelem * log(rss / nelem)
}

.aa_ic_n_entries <- function(object, X) {
    if (isTRUE(object[["fit_info"]][["missing"]])) {
        if (inherits(X, "sparseMatrix")) {
            return(sum(!is.na(X@x)))
        }
        return(sum(!is.na(X)))
    }

    if (inherits(X, "sparseMatrix")) prod(dim(X)) else sum(!is.na(X))
}

.aa_ic_rss <- function(object, X) {
    A <- .aa_input_coordinates_matrix(object)
    X_hat <- compositions(object) %*% A
    norm(X - X_hat, "F")^2
}

.aa_ic_final_loss <- function(object) {
    loss <- object[["loss"]]
    value <- utils::tail(loss[["loss"]], 1L)
    if (!is_number(value)) {
        stop("Information criteria require a finite final loss.", call. = FALSE)
    }
    value
}
