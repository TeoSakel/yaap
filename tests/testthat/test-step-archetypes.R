test_that("step_archetypes creates an untrained step", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L)

  expect_s3_class(rec$steps[[1L]], "step_archetypes")
  expect_false(rec$steps[[1L]]$trained)
  expect_equal(rec$steps[[1L]]$num_comp, 3L)
  expect_equal(rec$steps[[1L]]$fit_method, "pgd")
})

test_that("prep.step_archetypes trains the step", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L, seed = 42L) |>
    prep(training = toy_data())

  trained_step <- rec$steps[[1L]]
  expect_true(trained_step$trained)
  expect_s3_class(trained_step$res, "archetypes")
  expect_equal(nrow(trained_step$res[["coordinates"]]), 3L)
  expect_equal(unname(trained_step$col_names), colnames(toy_data()))
})

test_that("bake.step_archetypes returns composition columns", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L, seed = 42L) |>
    prep(training = toy_data())

  baked <- bake(rec, new_data = toy_data())
  expect_s3_class(baked, "tbl_df")
  # keep_original_cols = FALSE by default: original columns removed
  expect_false(any(colnames(toy_data()) %in% colnames(baked)))
  # composition columns present
  expect_true(all(c("A1", "A2", "A3") %in% colnames(baked)))
  # compositions are non-negative and row-stochastic
  S <- as.matrix(baked[, c("A1", "A2", "A3")])
  expect_true(all(S >= -1e-8))
  expect_true(all(abs(rowSums(S) - 1) < 1e-6))
})

test_that("keep_original_cols = TRUE retains original columns", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L,
                    keep_original_cols = TRUE, seed = 42L) |>
    prep(training = toy_data())

  baked <- bake(rec, new_data = toy_data())
  expect_true(all(colnames(toy_data()) %in% colnames(baked)))
  expect_true(all(c("A1", "A2", "A3") %in% colnames(baked)))
})

test_that("reconstruct = TRUE appends rec_ columns", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L,
                    reconstruct = TRUE, seed = 42L) |>
    prep(training = toy_data())

  baked <- bake(rec, new_data = toy_data())
  rec_cols <- paste0("rec_", colnames(toy_data()))
  expect_true(all(rec_cols %in% colnames(baked)))
})

test_that("step_archetypes works with nnls method", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L,
                    fit_method = "nnls", seed = 42L) |>
    prep(training = toy_data())

  baked <- bake(rec, new_data = toy_data())
  expect_true(all(c("A1", "A2", "A3") %in% colnames(baked)))
})

test_that("seed makes fitting reproducible", {
  skip_if_not_installed("recipes")
  library(recipes)

  prep1 <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L, seed = 123L) |>
    prep(training = toy_data())

  prep2 <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L, seed = 123L) |>
    prep(training = toy_data())

  coords1 <- prep1$steps[[1L]]$res[["coordinates"]]
  coords2 <- prep2$steps[[1L]]$res[["coordinates"]]
  expect_equal(coords1, coords2)
})

test_that("tidy.step_archetypes returns empty tibble when untrained", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L)

  out <- tidy(rec$steps[[1L]])
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0L)
  expect_true("id" %in% colnames(out))
})

test_that("tidy.step_archetypes returns archetype coordinates when trained", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L, seed = 42L) |>
    prep(training = toy_data())

  out <- tidy(rec$steps[[1L]])
  expect_s3_class(out, "tbl_df")
  expect_true(all(c("archetype", "term", "value", "id") %in% colnames(out)))
  expect_equal(nrow(out), 3L * ncol(toy_data()))
})

test_that("tunable.step_archetypes declares num_comp and delta", {
  skip_if_not_installed("recipes")
  skip_if_not_installed("tune")
  library(recipes); library(tune)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = tune(), delta = tune())

  out <- tunable(rec$steps[[1L]])
  expect_s3_class(out, "tbl_df")
  expect_setequal(out$name, c("num_comp", "delta"))
})

test_that("required_pkgs.step_archetypes returns correct packages", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L)

  out <- recipes::required_pkgs(rec$steps[[1L]])
  expect_setequal(out, c("YAAAP", "withr"))
})

test_that("print.step_archetypes does not error", {
  skip_if_not_installed("recipes")
  library(recipes)

  rec <- recipe(~., data = toy_data()) |>
    step_archetypes(all_numeric(), num_comp = 3L, seed = 42L)

  expect_no_error(print(rec$steps[[1L]]))

  rec_trained <- prep(rec, training = toy_data())
  expect_no_error(print(rec_trained$steps[[1L]]))
})
