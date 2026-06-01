# Information criteria for archetypes objects ---------------------------------

#' Information Criteria for Archetype Fits
#'
#' Computes AIC, BIC, or Mallows' Cp-like criteria for fitted archetype models
#' with selectable degrees-of-freedom estimates.
#'
#' @param object An object of class `archetypes`.
#' @param type Criterion to compute: `"aic"`, `"bic"`, or `"cp"`.
#' @param df_method Degrees-of-freedom estimate. `"suleman"` preserves the
#'   adapted efficiency-corrected criterion used historically by `AIC()`;
#'   `"covariance"` uses Stein's frozen-smoother formula; `"parametric"`
#'   uses the full archetype parameter count.
#' @param n_eff Effective sample size used in standard AIC/BIC/Cp penalties:
#'   `"samples"` uses rows of the training matrix; `"entries"` uses finite
#'   matrix entries.
#' @param sigma2 Optional Gaussian variance estimate for `type = "cp"`.
#' @param ... Additional arguments passed to `aa_ic()`.
#'
#' @return A numeric scalar.
#'
#' @details
#' `AIC.archetypes()` defaults to the Suleman adapted criterion for backward
#' compatibility. `BIC.archetypes()` uses the analogous Suleman-style penalty by
#' default. For `df_method = "covariance"`, the fitted archetype weights are held
#' fixed and the fitted-value map is treated as `X_hat = H X`, with
#' `H = compositions(object) %*% coefficients(object)`. The degrees of freedom
#' are then `M * sum(diag(H))`, where `M` is the number of features.
#'
#' For non-Suleman Gaussian criteria, AIC and BIC use the profile log-RSS scale.
#' For non-Gaussian PAA families, AIC and BIC use twice the final optimized
#' objective plus the selected penalty. Cp is only defined for Gaussian fits.
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
                  type = c("aic", "bic", "cp"),
                  df_method = c("suleman", "covariance", "parametric"),
                  n_eff = c("samples", "entries"),
                  sigma2 = NULL) {
    if (!inherits(object, "archetypes")) {
        stop("`object` must be an `archetypes` object.", call. = FALSE)
    }
    type <- match.arg(type)
    df_method <- match.arg(df_method)
    n_eff <- match.arg(n_eff)

    X <- .aa_ic_data_matrix(object, type)
    family <- object[["family"]] %||% "gaussian"
    if (identical(df_method, "suleman")) {
        return(.aa_ic_suleman(object, X, type, n_eff))
    }

    df <- .aa_ic_df(object, X, df_method)
    n <- .aa_ic_n_eff(X, n_eff)
    if (identical(type, "cp")) {
        if (!identical(family, "gaussian")) {
            stop("Cp is only defined for Gaussian archetypes objects.", call. = FALSE)
        }
        rss <- .aa_ic_rss(object, X)
        sigma2 <- .aa_ic_sigma2(object, sigma2)
        return(rss / sigma2 - n + 2 * df)
    }

    penalty <- switch(type,
        aic = 2,
        bic = log(n)
    )
    if (identical(family, "gaussian")) {
        rss <- .aa_ic_rss(object, X)
        if (!is_non_negative(rss)) {
            stop("Gaussian information criteria require finite non-negative RSS.", call. = FALSE)
        }
        return(n * log(rss / n) + penalty * df)
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

.aa_ic_data_matrix <- function(object, type) {
    X <- object[["data"]]
    if (is.null(X)) {
        stop(paste(
            toupper(type), "requires original data `X`;",
            "provide it when constructing the archetypes object."
        ), call. = FALSE)
    }
    if (inherits(X, "fd")) {
        return(.aa_fd_to_matrix(X))
    }
    as.matrix(X)
}

.aa_ic_suleman <- function(object, X, type, n_eff) {
    if (identical(type, "cp")) {
        stop("Cp is not defined for `df_method = 'suleman'`.", call. = FALSE)
    }

    family <- object[["family"]] %||% "gaussian"
    if (!identical(family, "gaussian")) {
        stop(
            "Suleman information criteria are not defined for non-Gaussian archetypes objects.",
            call. = FALSE
        )
    }
    if (any(object[["slack"]] > 0)) {
        warning(paste(
            toupper(type), "computation assumes coefficients are row-stochastic;",
            "slack > 0 may violate this assumption."
        ), call. = FALSE)
    }

    S <- compositions(object)
    K <- ncol(S)
    N <- nrow(X)
    M <- ncol(X)
    if (N <= M) {
        warning(
            paste(
                "Adapted", toupper(type), "is undefined when the number of samples",
                "is not larger than the number of features; returning NA."
            ),
            call. = FALSE
        )
        return(NA_real_)
    }

    A <- .aa_input_coordinates_matrix(object)
    X_hat <- S %*% A
    eta <- tryCatch(.aa_effic(X, X_hat), error = function(e) NA_real_)
    if (!is.finite(eta) || eta <= 0) {
        warning(
            paste(
                "Adapted", toupper(type), "is undefined because the efficiency term",
                "is non-positive or non-finite; returning NA."
            ),
            call. = FALSE
        )
        return(NA_real_)
    }

    rss <- norm(X - X_hat, "F")^2
    nelem <- prod(dim(X))
    npar <- .aa_ic_parametric_df(N, K)
    penalty <- switch(type,
        aic = 2,
        bic = log(.aa_ic_n_eff(X, n_eff))
    )
    log(rss / nelem) + penalty * npar / (N * eta)
}

.aa_ic_df <- function(object, X, df_method) {
    K <- nrow(coefficients(object))
    N <- nrow(X)
    switch(df_method,
        parametric = .aa_ic_parametric_df(N, K),
        covariance = .aa_ic_covariance_df(object, X)
    )
}

.aa_ic_parametric_df <- function(N, K) {
    N * (K - 1L) + K * (N - 1L) + 1L
}

.aa_ic_covariance_df <- function(object, X) {
    S <- compositions(object)
    B <- coefficients(object)
    M <- ncol(X)
    M * sum(tcrossprod(S, B))
}

.aa_ic_n_eff <- function(X, n_eff) {
    out <- switch(n_eff,
        samples = nrow(X),
        entries = if (inherits(X, "sparseMatrix")) length(X@x) else sum(!is.na(X))
    )
    if (!is_positive(out)) {
        stop("`n_eff` must resolve to a positive finite value.", call. = FALSE)
    }
    out
}

.aa_ic_rss <- function(object, X) {
    A <- .aa_input_coordinates_matrix(object)
    X_hat <- compositions(object) %*% A
    norm(X - X_hat, "F")^2
}

.aa_ic_sigma2 <- function(object, sigma2) {
    if (is.null(sigma2)) {
        loss <- object[["loss"]]
        if (!is.null(loss[["sigma2"]])) {
            sigma2 <- utils::tail(loss[["sigma2"]], 1L)
        }
    }
    if (is.null(sigma2)) {
        stop("`sigma2` must be supplied for Cp when no finite `loss$sigma2` is stored.", call. = FALSE)
    }
    if (!is_positive(sigma2)) {
        stop("`sigma2` must be a single positive finite number.", call. = FALSE)
    }
    sigma2
}

.aa_ic_final_loss <- function(object) {
    loss <- object[["loss"]]
    value <- utils::tail(loss[["loss"]], 1L)
    if (!is_number(value)) {
        stop("Information criteria require a finite final loss.", call. = FALSE)
    }
    value
}
