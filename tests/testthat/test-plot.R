test_that("plot.archetypes smoke tests supported plot modes", {
    fit <- manual_fit()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_identical(plot(fit, "loss"), fit)
    expect_identical(plot(fit, "compositions"), fit)
    expect_identical(plot(fit, "composition"), fit)
    expect_identical(plot(fit, "composision"), fit)
    expect_identical(plot(fit, "ternary"), fit)
    expect_identical(plot(fit, "simplex"), fit)
    expect_identical(plot(fit, "profiles"), fit)
    expect_identical(plot(fit, "coordinates"), fit)
})

test_that("composition_barplot supports matrix-like inputs and clustering", {
    fit <- manual_fit()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    out <- composition_barplot(fit[["compositions"]], legend = FALSE)
    expect_equal(out[["row_order"]], seq_len(nrow(fit[["compositions"]])))
    expect_equal(out[["col_order"]], seq_len(ncol(fit[["compositions"]])))
    expect_null(out[["row_hclust"]])
    expect_null(out[["col_hclust"]])

    clustered <- composition_barplot(
        as.data.frame(fit[["compositions"]]),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        legend = FALSE
    )
    expect_s3_class(clustered[["row_hclust"]], "hclust")
    expect_s3_class(clustered[["col_hclust"]], "hclust")
    expect_equal(sort(clustered[["row_order"]]), seq_len(nrow(fit[["compositions"]])))
    expect_equal(sort(clustered[["col_order"]]), seq_len(ncol(fit[["compositions"]])))
})

test_that("composition_barplot exposes clustering distance and linkage", {
    S <- matrix(
        c(
            0.70, 0.20, 0.10,
            0.55, 0.35, 0.10,
            0.15, 0.70, 0.15,
            0.10, 0.25, 0.65
        ),
        ncol = 3L,
        byrow = TRUE
    )

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    clustered <- composition_barplot(
        S,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        distance_rows = "manhattan",
        distance_cols = "correlation",
        linkage = "average",
        legend = FALSE
    )
    expect_identical(clustered[["row_hclust"]][["method"]], "average")
    expect_identical(clustered[["col_hclust"]][["method"]], "average")

    both <- composition_barplot(
        S,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        distance = "maximum",
        linkage = "single",
        legend = FALSE
    )
    expect_identical(both[["row_hclust"]][["method"]], "single")
    expect_identical(both[["col_hclust"]][["method"]], "single")
})

test_that("plot.archetypes handles higher-dimensional coordinate projections", {
    fit <- manual_fit()
    X <- cbind(fit[["data"]], z = c(0, 1, 1, 0.2), w = c(1, 0, 1, 0.4))
    fit[["data"]] <- X
    fit[["coordinates"]] <- X[1:3, , drop = FALSE]
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_identical(plot(fit, "coordinates"), fit)
    expect_identical(plot(fit, "coordinates", projection = "pca"), fit)
})

test_that("plot.archetypes can subset samples in observation-level plots", {
    fit <- manual_fit()
    sample_names <- paste0("s", seq_len(nrow(fit[["compositions"]])))
    rownames(fit[["compositions"]]) <- sample_names
    rownames(fit[["data"]]) <- sample_names

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_identical(plot(fit, "compositions", samples = c(1L, 3L), legend = FALSE), fit)
    expect_identical(plot(fit, "ternary", samples = c("s2", "s4")), fit)
    expect_identical(plot(fit, "coordinates", samples = c(TRUE, FALSE, TRUE, FALSE)), fit)
    expect_identical(plot(fit, "coordinates", projection = "pca", samples = c("s1", "s4")), fit)

    expect_error(
        plot(fit, "coordinates", samples = c(TRUE, FALSE)),
        "Logical `samples`"
    )
    expect_error(
        plot(fit, "ternary", samples = "missing"),
        "Some `samples`"
    )
})
