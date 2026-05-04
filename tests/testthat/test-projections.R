test_that("projection helpers return row-stochastic matrices", {
    mat <- matrix(
        c(
            -1, 2, 0.5,
            0.2, 0.2, 0.2,
            4, 1, -2
        ),
        nrow = 3,
        byrow = TRUE
    )

    simplex <- proj_simplex(mat)
    l1 <- proj_l1(mat)

    expect_matrix_dim(simplex, 3L, 3L)
    expect_matrix_dim(l1, 3L, 3L)
    expect_row_stochastic(simplex)
    expect_row_stochastic(l1)
})

test_that("fit_simplex projects data onto archetype simplex", {
    A <- matrix(
        c(
            0, 0,
            1, 0,
            0, 1
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("A1", "A2", "A3"), c("x", "y"))
    )
    X <- matrix(
        c(
            0.2, 0.2,
            0.8, 0.1
        ),
        ncol = 2,
        byrow = TRUE
    )

    S <- fit_simplex(A, X)

    expect_matrix_dim(S, 2L, 3L)
    expect_named(as.data.frame(S), rownames(A))
    expect_row_stochastic(S)
})
