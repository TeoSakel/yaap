test_that("plot.archetypes smoke tests supported plot modes", {
    fit <- manual_fit()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_named(plot(fit, "loss"), c("iteration", "loss", "plot_args"))
    expect_named(plot(fit, "compositions"), c(
        "compositions", "row_order", "col_order", "row_hclust",
        "col_hclust", "col", "barplot_args"
    ))
    expect_named(plot(fit, "composition"), c(
        "compositions", "row_order", "col_order", "row_hclust",
        "col_hclust", "col", "barplot_args"
    ))
    expect_named(plot(fit, "composision"), c(
        "compositions", "row_order", "col_order", "row_hclust",
        "col_hclust", "col", "barplot_args"
    ))
    expect_named(plot(fit, "ternary"), c("compositions", "plot_args"))
    expect_named(plot(fit, "simplex"), c("compositions", "plot_args"))
    expect_named(plot(fit, "profiles"), c(
        "coordinates", "family", "archetype_names", "barplot_args"
    ))
    expect_named(plot(fit, "coordinates"), c(
        "coordinates", "data", "projection", "pca", "archetype_names",
        "show_anames", "args.coordinates", "args.data.scatter", "plot_args"
    ))
})

test_that("plot methods leave coordinate-helper arguments in dots", {
    helper_args <- c("data", "projection", "show_anames", "args.data.scatter")
    expect_false(any(helper_args %in% names(formals(plot.archetypes))))
    expect_false(any(helper_args %in% names(formals(plot.kernel_archetypes))))
})

test_that("plot.archetypes profiles uses fixed height but accepts other barplot args", {
    fit <- manual_fit()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    bad_height <- matrix(1, nrow = 2L, ncol = 2L)
    out <- plot(
        fit,
        "profiles",
        height = bad_height,
        col = c("red", "blue", "green"),
        horiz = TRUE,
        border = NA
    )
    expect_equal(out[["barplot_args"]][["height"]], fit[["coordinates"]])
    expect_true(out[["barplot_args"]][["horiz"]])
})

test_that("plot_archetypes_compositions supports matrix-like inputs and clustering", {
    fit <- manual_fit()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    out <- plot_archetypes_compositions(fit[["compositions"]], legend = FALSE)
    expect_equal(out[["row_order"]], seq_len(nrow(fit[["compositions"]])))
    expect_equal(out[["col_order"]], seq_len(ncol(fit[["compositions"]])))
    expect_null(out[["row_hclust"]])
    expect_null(out[["col_hclust"]])

    clustered <- plot_archetypes_compositions(
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

test_that("plot_archetypes_compositions exposes clustering distance and linkage", {
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

    clustered <- plot_archetypes_compositions(
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

    both <- plot_archetypes_compositions(
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

test_that("plot_archetypes_compositions supports PC1 and AOP ordering", {
    S <- matrix(
        c(
            0.70, 0.20, 0.10,
            0.55, 0.35, 0.10,
            0.15, 0.70, 0.15,
            0.10, 0.25, 0.65,
            0.25, 0.50, 0.25
        ),
        ncol = 3L,
        byrow = TRUE
    )

    pca <- stats::prcomp(S, center = TRUE, scale. = FALSE)

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    pc1 <- plot_archetypes_compositions(
        S,
        cluster_rows = "PC1",
        cluster_cols = "PC1",
        legend = FALSE
    )
    expect_equal(pc1[["row_order"]], order(pca[["x"]][, 1L]))
    expect_equal(pc1[["col_order"]], order(pca[["rotation"]][, 1L]))
    expect_equal(pc1[["row_hclust"]][["order"]], pc1[["row_order"]])
    expect_equal(pc1[["col_hclust"]][["order"]], pc1[["col_order"]])

    aop <- plot_archetypes_compositions(
        S,
        cluster_rows = "AOP",
        cluster_cols = "AOP",
        legend = FALSE
    )
    expect_equal(aop[["row_order"]], order(atan2(pca[["x"]][, 2L], pca[["x"]][, 1L])))
    expect_equal(
        aop[["col_order"]],
        order(atan2(pca[["rotation"]][, 2L], pca[["rotation"]][, 1L]))
    )
    expect_equal(aop[["row_hclust"]][["order"]], aop[["row_order"]])
    expect_equal(aop[["col_hclust"]][["order"]], aop[["col_order"]])

    expect_error(
        plot_archetypes_compositions(S, cluster_rows = "bad", legend = FALSE),
        "cluster_rows"
    )
    expect_error(
        plot_archetypes_compositions(S, cluster_cols = "bad", legend = FALSE),
        "cluster_cols"
    )
})

test_that("plot helpers support plot = FALSE", {
    fit <- manual_fit()

    comp <- plot_archetypes_compositions(fit[["compositions"]], plot = FALSE)
    expect_equal(unname(comp[["compositions"]]), unname(fit[["compositions"]]))

    loss <- plot_archetypes_loss(fit[["loss"]], plot = FALSE)
    expect_equal(loss[["loss"]], fit[["loss"]][["loss"]])

    profiles <- plot_archetypes_profiles(
        fit[["coordinates"]],
        archetype_names = anames(fit),
        plot = FALSE
    )
    expect_equal(profiles[["coordinates"]], fit[["coordinates"]])

    coords <- plot_archetypes_coordinates(
        fit[["coordinates"]],
        data = fit[["data"]],
        archetype_names = anames(fit),
        plot = FALSE
    )
    expect_equal(coords[["coordinates"]], fit[["coordinates"]])
    expect_equal(coords[["data"]], fit[["data"]])

    expect_named(plot(fit, "coordinates", plot = FALSE), c(
        "coordinates", "data", "projection", "pca", "archetype_names",
        "show_anames", "args.coordinates", "args.data.scatter", "plot_args"
    ))
})

test_that("plot_archetypes_coordinates handles coordinates-only and argument routing", {
    fit <- manual_fit()
    group <- c("g1", "g2", "g1", "g2")
    data_col <- c(g1 = "#1b9e77", g2 = "#d95f02")[group]

    coords_only <- plot_archetypes_coordinates(fit[["coordinates"]], plot = FALSE)
    expect_null(coords_only[["data"]])
    expect_equal(coords_only[["projection"]], "none")

    expect_error(
        plot_archetypes_coordinates(fit[["coordinates"]], projection = "pca", plot = FALSE),
        "projection"
    )

    routed <- plot_archetypes_coordinates(
        fit[["coordinates"]],
        data = fit[["data"]],
        col = "black",
        pch = 17,
        cex = 1.4,
        args.data.scatter = list(col = data_col, pch = c(1, 2, 3, 4), cex = rep(0.8, 4)),
        plot = FALSE
    )
    expect_identical(routed[["args.coordinates"]][["col"]], "black")
    expect_identical(routed[["args.coordinates"]][["pch"]], 17)
    expect_identical(routed[["args.data.scatter"]][["col"]], data_col)
    expect_identical(routed[["args.data.scatter"]][["pch"]], c(1, 2, 3, 4))
})

test_that("plot.archetypes handles higher-dimensional coordinate projections", {
    fit <- manual_fit()
    X <- cbind(fit[["data"]], z = c(0, 1, 1, 0.2), w = c(1, 0, 1, 0.4))
    fit[["data"]] <- X
    fit[["coordinates"]] <- X[1:3, , drop = FALSE]
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_equal(plot(fit, "coordinates")[["projection"]], "none")
    pca <- plot(fit, "coordinates", projection = "pca")
    expect_equal(pca[["projection"]], "pca")
    expect_s3_class(pca[["pca"]], "prcomp")
    expect_false(plot(fit, "coordinates", show_anames = FALSE)[["show_anames"]])
    expect_false(
        plot(fit, "coordinates", projection = "pca", show_anames = FALSE)[["show_anames"]]
    )
})

test_that("plot.archetypes coordinates supports data vectors with args.data.scatter", {
    fit <- manual_fit()
    group <- c("g1", "g2", "g1", "g2")
    data_col <- c(g1 = "#1b9e77", g2 = "#d95f02")[group]

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    out <- plot(
        fit,
        "coordinates",
        col = "black",
        pch = 17,
        cex = 1.4,
        args.data.scatter = list(col = data_col, pch = c(1, 2, 3, 4), cex = rep(0.8, 4))
    )
    expect_identical(out[["args.coordinates"]][["col"]], "black")
    expect_identical(out[["args.coordinates"]][["pch"]], 17)
    expect_identical(out[["args.data.scatter"]][["col"]], data_col)
})

test_that("plot.archetypes coordinate defaults can be overridden by helper args", {
    fit <- manual_fit()

    no_names <- plot(fit, "coordinates", archetype_names = NULL, plot = FALSE)
    expect_null(no_names[["archetype_names"]])

    custom_names <- c("left", "middle", "right")
    named <- plot(fit, "coordinates", archetype_names = custom_names, plot = FALSE)
    expect_identical(named[["archetype_names"]], custom_names)

    no_profile_names <- plot(fit, "profiles", archetype_names = NULL, plot = FALSE)
    expect_null(no_profile_names[["archetype_names"]])
})

test_that("plot.archetypes can subset samples in observation-level plots", {
    fit <- manual_fit()
    sample_names <- paste0("s", seq_len(nrow(fit[["compositions"]])))
    rownames(fit[["compositions"]]) <- sample_names
    rownames(fit[["data"]]) <- sample_names

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    comp <- plot(fit, "compositions", subset = c(1L, 3L), legend = FALSE)
    expect_equal(rownames(comp[["compositions"]]), sample_names[c(1L, 3L)])
    simplex <- plot(fit, "ternary", subset = c("s2", "s4"))
    expect_equal(rownames(simplex[["compositions"]]), c("s2", "s4"))
    coords <- plot(fit, "coordinates", subset = c(TRUE, FALSE, TRUE, FALSE))
    expect_equal(rownames(coords[["data"]]), sample_names[c(1L, 3L)])
    pca <- plot(fit, "coordinates", projection = "pca", subset = c("s1", "s4"))
    expect_equal(nrow(pca[["data"]]), 2L)

    expect_error(
        plot(fit, "coordinates", subset = c(TRUE, FALSE)),
        "Logical `subset`"
    )
    expect_error(
        plot(fit, "ternary", subset = "missing"),
        "Some `subset`"
    )
})
