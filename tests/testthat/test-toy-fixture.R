test_that("toy CSV fixture is available and well formed", {
    toy <- toy_data()

    expect_equal(names(toy), c("x", "y"))
    expect_equal(nrow(toy), 250L)
    expect_true(all(vapply(toy, is.numeric, logical(1))))
    expect_false(anyNA(toy))
})
