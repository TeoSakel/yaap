test_that("anames gets and sets archetype labels without hiding list components", {
    A <- matrix(
        c(0, 0, 1, 0, 0, 1),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(paste0("A", 1:3), c("x", "y"))
    )
    B <- diag(3L)
    rownames(B) <- rownames(A)
    colnames(B) <- paste0("s", 1:3)
    S <- diag(3L)
    rownames(S) <- colnames(B)
    colnames(S) <- rownames(A)

    fit <- archetypes(
        A = A,
        coefficients = B,
        compositions = S,
        init = A,
        data = A,
        feature_map = .aa_identity_feature_map(A)
    )

    expect_equal(anames(fit), c("A1", "A3", "A2"))
    expect_true(all(c("A", "coefficients", "compositions") %in% names(fit)))

    anames(fit) <- c("low", "mid", "high")

    expect_equal(anames(fit), c("low", "mid", "high"))
    expect_equal(rownames(fit[["coefficients"]]), anames(fit))
    expect_equal(colnames(fit[["compositions"]]), anames(fit))
    expect_equal(rownames(fit[["init"]]), anames(fit))
    expect_true(all(c("A", "coefficients", "compositions") %in% names(fit)))

    expect_error(anames(fit) <- c("a", "b"), "Expected 3 archetype names")
    expect_error(anames(fit) <- c("a", "a", "b"), "unique")
    expect_error(anames(fit) <- c("a", "", "b"), "empty")
    expect_error(anames(fit) <- c("a", NA, "b"), "missing")
})

test_that("archetypes constructor canonicalizes archetype order", {
    X <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1,
            0.25, 0.25
        ),
        ncol = 2L,
        byrow = TRUE,
        dimnames = list(paste0("s", 1:4), c("x", "y"))
    )
    A <- matrix(
        c(
            1, 0,
            0, 1,
            0, 0
        ),
        ncol = 2L,
        byrow = TRUE,
        dimnames = list(c("right", "top", "origin"), colnames(X))
    )
    B <- matrix(
        c(
            0, 1, 0, 0,
            0, 0, 1, 0,
            1, 0, 0, 0
        ),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(rownames(A), rownames(X))
    )
    S <- matrix(
        c(
            0, 0, 1,
            1, 0, 0,
            0, 1, 0,
            0.25, 0.25, 0.5
        ),
        nrow = 4L,
        byrow = TRUE,
        dimnames = list(rownames(X), rownames(A))
    )

    fit <- archetypes(
        A = A,
        coefficients = B,
        compositions = S,
        init = A,
        data = X,
        feature_map = .aa_identity_feature_map(A)
    )

    expect_equal(anames(fit), c("origin", "top", "right"))
    expect_equal(unname(coordinates(fit)), matrix(c(0, 0, 0, 1, 1, 0), ncol = 2L, byrow = TRUE))
    expect_equal(rownames(fit[["coefficients"]]), anames(fit))
    expect_equal(colnames(fit[["compositions"]]), anames(fit))
    expect_equal(rownames(fit[["init"]]), anames(fit))
    expect_equal(rownames(fit[["A"]]), anames(fit))
    expect_equal(unname(fitted(fit)), unname(S %*% A))
})

test_that("archetypes constructor canonicalizes without coordinates", {
    B <- matrix(
        c(
            0, 1, 0,
            1, 0, 0,
            0, 0, 1
        ),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(c("second", "first", "third"), paste0("s", 1:3))
    )
    S <- diag(3L)
    colnames(S) <- rownames(B)
    rownames(S) <- colnames(B)

    fit <- archetypes(
        A = NULL,
        coefficients = B,
        compositions = S,
        feature_map = .aa_unsupported_feature_map("manual")
    )

    expect_equal(anames(fit), c("third", "second", "first"))
    expect_equal(rownames(fit[["coefficients"]]), anames(fit))
    expect_equal(colnames(fit[["compositions"]]), anames(fit))
})

test_that("anames gets and sets kernel archetype labels", {
    B <- diag(3L)
    rownames(B) <- paste0("A", 1:3)
    colnames(B) <- paste0("s", 1:3)
    S <- diag(3L)
    rownames(S) <- colnames(B)
    colnames(S) <- rownames(B)
    G <- diag(3L)
    rownames(G) <- colnames(G) <- colnames(B)

    fit <- .aa_new_kernel_archetypes(
        coefficients = B,
        compositions = S,
        gram = G,
        init = B
    )

    expect_equal(rownames(fit[["init"]]), anames(fit))

    anames(fit) <- c("red", "green", "blue")

    expect_equal(anames(fit), c("red", "green", "blue"))
    expect_equal(rownames(fit[["coefficients"]]), anames(fit))
    expect_equal(colnames(fit[["compositions"]]), anames(fit))
    expect_true(all(c("coefficients", "compositions", "gram") %in% names(fit)))
})
