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
        coordinates = A,
        coefficients = B,
        compositions = S,
        init = A,
        data = A
    )

    expect_equal(anames(fit), c("A1", "A2", "A3"))
    expect_true(all(c("coordinates", "coefficients", "compositions") %in% names(fit)))

    anames(fit) <- c("low", "mid", "high")

    expect_equal(anames(fit), c("low", "mid", "high"))
    expect_equal(rownames(fit[["coordinates"]]), anames(fit))
    expect_equal(rownames(fit[["coefficients"]]), anames(fit))
    expect_equal(colnames(fit[["compositions"]]), anames(fit))
    expect_equal(rownames(fit[["init"]]), anames(fit))
    expect_true(all(c("coordinates", "coefficients", "compositions") %in% names(fit)))

    expect_error(anames(fit) <- c("a", "b"), "Expected 3 archetype names")
    expect_error(anames(fit) <- c("a", "a", "b"), "unique")
    expect_error(anames(fit) <- c("a", "", "b"), "empty")
    expect_error(anames(fit) <- c("a", NA, "b"), "missing")
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
    P <- B

    fit <- kernel_archetypes(
        coefficients = B,
        compositions = S,
        gram = G,
        coordinates_proxy = P
    )

    anames(fit) <- c("red", "green", "blue")

    expect_equal(anames(fit), c("red", "green", "blue"))
    expect_equal(rownames(fit[["coefficients"]]), anames(fit))
    expect_equal(colnames(fit[["compositions"]]), anames(fit))
    expect_equal(rownames(fit[["coordinates_proxy"]]), anames(fit))
    expect_true(all(c("coefficients", "compositions", "gram") %in% names(fit)))
})
