test_that("plot.archetypes smoke tests supported plot modes", {
    fit <- manual_fit()
    local_test_pdf()

    expect_named(plot(fit, "loss"), c("iteration", "loss"))
    expect_named(plot(fit, "compositions"), c("sample", "archetype", "weight"))
    expect_named(plot(fit, "composition"), c("sample", "archetype", "weight"))
    expect_named(plot(fit, "composision"), c("sample", "archetype", "weight"))
    expect_s3_class(plot(fit, "ternary"), "data.frame")
    expect_s3_class(plot(fit, "simplex"), "data.frame")
    expect_named(plot(fit, "profiles"), c("archetype", "feature", "value"))
    expect_named(plot(fit, "coordinates"), c("x", "y", "name", "archetype"))
})

test_that("plot methods leave coordinate-helper arguments in dots", {
    helper_args <- c("data", "projection", "show_anames", "args.data.scatter")
    expect_false(any(helper_args %in% names(formals(plot.archetypes))))
    expect_false(any(helper_args %in% names(formals(plot.kernel_archetypes))))
})

test_that("plot.archetypes profiles uses fixed height but accepts other barplot args", {
    fit <- manual_fit()
    local_test_pdf()

    bad_height <- matrix(1, nrow = 2L, ncol = 2L)
    out <- plot(
        fit,
        "profiles",
        height = bad_height,
        col = c("red", "blue", "green"),
        horiz = TRUE,
        border = NA
    )
    expect_named(out, c("archetype", "feature", "value"))
    expect_equal(nrow(out), length(coordinates(fit)))
    expect_equal(out[["value"]], as.vector(coordinates(fit)))
})

test_that("plot_archetypes_compositions supports matrix-like inputs and clustering", {
    fit <- manual_fit()
    local_test_pdf()

    out <- plot_archetypes_compositions(fit[["compositions"]], legend = FALSE)
    expect_named(out, c("sample", "archetype", "weight"))
    expect_equal(levels(out[["sample"]]), as.character(seq_len(nrow(fit[["compositions"]]))))
    expect_equal(levels(out[["archetype"]]), colnames(fit[["compositions"]]))
    expect_equal(out[["weight"]], as.vector(fit[["compositions"]]))

    clustered <- plot_archetypes_compositions(
        as.data.frame(fit[["compositions"]]),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        legend = FALSE
    )
    expect_named(clustered, c("sample", "archetype", "weight"))
    expect_equal(sort(levels(clustered[["sample"]])), as.character(seq_len(nrow(fit[["compositions"]]))))
    expect_equal(sort(levels(clustered[["archetype"]])), sort(colnames(fit[["compositions"]])))
    expect_equal(nrow(clustered), length(fit[["compositions"]]))
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

    local_test_pdf()

    clustered <- plot_archetypes_compositions(
        S,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        distance_rows = "manhattan",
        distance_cols = "correlation",
        linkage = "average",
        legend = FALSE
    )
    expect_named(clustered, c("sample", "archetype", "weight"))
    expect_equal(sort(levels(clustered[["sample"]])), as.character(seq_len(nrow(S))))
    expect_equal(sort(levels(clustered[["archetype"]])), paste0("A", seq_len(ncol(S))))

    both <- plot_archetypes_compositions(
        S,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        distance = "maximum",
        linkage = "single",
        legend = FALSE
    )
    expect_named(both, c("sample", "archetype", "weight"))
    expect_equal(sort(levels(both[["sample"]])), as.character(seq_len(nrow(S))))
    expect_equal(sort(levels(both[["archetype"]])), paste0("A", seq_len(ncol(S))))
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

    local_test_pdf()

    pc1 <- plot_archetypes_compositions(
        S,
        cluster_rows = "PC1",
        cluster_cols = "PC1",
        legend = FALSE
    )
    expect_equal(levels(pc1[["sample"]]), as.character(order(pca[["x"]][, 1L])))
    expect_equal(levels(pc1[["archetype"]]), paste0("A", order(pca[["rotation"]][, 1L])))

    aop <- plot_archetypes_compositions(
        S,
        cluster_rows = "AOP",
        cluster_cols = "AOP",
        legend = FALSE
    )
    expect_equal(
        levels(aop[["sample"]]),
        as.character(order(atan2(pca[["x"]][, 2L], pca[["x"]][, 1L])))
    )
    expect_equal(
        levels(aop[["archetype"]]),
        paste0("A", order(atan2(pca[["rotation"]][, 2L], pca[["rotation"]][, 1L])))
    )

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
    local_test_pdf()

    comp_visible <- withVisible(plot_archetypes_compositions(fit[["compositions"]], plot = FALSE))
    expect_true(comp_visible[["visible"]])
    comp <- comp_visible[["value"]]
    expect_named(comp, c("sample", "archetype", "weight"))
    expect_equal(comp[["weight"]], as.vector(fit[["compositions"]]))

    loss_visible <- withVisible(plot_archetypes_loss(fit[["loss"]], plot = FALSE))
    expect_true(loss_visible[["visible"]])
    loss <- loss_visible[["value"]]
    expect_equal(loss[["loss"]], fit[["loss"]][["loss"]])

    profiles <- plot_archetypes_profiles(
        coordinates(fit),
        archetype_names = anames(fit),
        plot = FALSE
    )
    expect_named(profiles, c("archetype", "feature", "value"))
    expect_equal(profiles[["value"]], as.vector(coordinates(fit)))

    coords_visible <- withVisible(plot_archetypes_coordinates(
        coordinates(fit),
        data = fit[["data"]],
        archetype_names = anames(fit),
        plot = FALSE
    ))
    expect_true(coords_visible[["visible"]])
    coords <- coords_visible[["value"]]
    expect_named(coords, c("x", "y", "name", "archetype"))
    expect_equal(unname(as.matrix(coords[!coords[["archetype"]], c("x", "y")])), unname(fit[["data"]]))
    expect_equal(unname(as.matrix(coords[coords[["archetype"]], c("x", "y")])), unname(coordinates(fit)))

    visible <- withVisible(plot(fit, "coordinates", plot = FALSE))
    expect_true(visible[["visible"]])
    expect_named(visible[["value"]], c("x", "y", "name", "archetype"))

    drawn <- withVisible(plot(fit, "coordinates"))
    expect_false(drawn[["visible"]])
    expect_named(drawn[["value"]], c("x", "y", "name", "archetype"))
})

test_that("plot_archetypes_coordinates handles coordinates-only and argument routing", {
    fit <- manual_fit()
    group <- c("g1", "g2", "g1", "g2")
    data_col <- c(g1 = "#1b9e77", g2 = "#d95f02")[group]

    coords_only <- plot_archetypes_coordinates(coordinates(fit), plot = FALSE)
    expect_named(coords_only, c("x", "y", "name", "archetype"))
    expect_true(all(coords_only[["archetype"]]))

    expect_error(
        plot_archetypes_coordinates(coordinates(fit), projection = "pca", plot = FALSE),
        "projection"
    )

    routed <- plot_archetypes_coordinates(
        coordinates(fit),
        data = fit[["data"]],
        col = "black",
        pch = 17,
        cex = 1.4,
        args.data.scatter = list(col = data_col, pch = c(1, 2, 3, 4), cex = rep(0.8, 4)),
        plot = FALSE
    )
    expect_named(routed, c("x", "y", "name", "archetype"))
    expect_equal(sum(!routed[["archetype"]]), nrow(fit[["data"]]))
    expect_equal(sum(routed[["archetype"]]), nrow(coordinates(fit)))
})

test_that("plot.archetypes handles higher-dimensional coordinate projections", {
    fit <- manual_fit()
    X <- cbind(fit[["data"]], z = c(0, 1, 1, 0.2), w = c(1, 0, 1, 0.4))
    fit[["data"]] <- X
    fit[["A"]] <- X[1:3, , drop = FALSE]
    fit[["feature_map"]] <- .aa_identity_feature_map(fit[["A"]])
    local_test_pdf()

    expect_named(plot(fit, "coordinates"), c("x", "y", "z", "w", "name", "archetype"))
    pca <- plot(fit, "coordinates", projection = "pca")
    expect_named(pca, c("PC1", "PC2", "name", "archetype"))
    expect_equal(nrow(pca), nrow(fit[["data"]]) + nrow(coordinates(fit)))
    expect_named(plot(fit, "coordinates", show_anames = FALSE), c("x", "y", "z", "w", "name", "archetype"))
    expect_named(plot(fit, "coordinates", projection = "pca", show_anames = FALSE), c("PC1", "PC2", "name", "archetype"))
})

test_that("plot.archetypes coordinates supports data vectors with args.data.scatter", {
    fit <- manual_fit()
    group <- c("g1", "g2", "g1", "g2")
    data_col <- c(g1 = "#1b9e77", g2 = "#d95f02")[group]

    local_test_pdf()

    out <- plot(
        fit,
        "coordinates",
        col = "black",
        pch = 17,
        cex = 1.4,
        args.data.scatter = list(col = data_col, pch = c(1, 2, 3, 4), cex = rep(0.8, 4))
    )
    expect_named(out, c("x", "y", "name", "archetype"))
    expect_equal(sum(!out[["archetype"]]), nrow(fit[["data"]]))
    expect_equal(sum(out[["archetype"]]), nrow(coordinates(fit)))
})

test_that("plot.archetypes coordinate defaults can be overridden by helper args", {
    fit <- manual_fit()

    no_names <- plot(fit, "coordinates", archetype_names = NULL, plot = FALSE)
    expect_equal(no_names[["name"]][no_names[["archetype"]]], rownames(coordinates(fit)))

    custom_names <- c("left", "middle", "right")
    named <- plot(fit, "coordinates", archetype_names = custom_names, plot = FALSE)
    expect_identical(named[["name"]][named[["archetype"]]], custom_names)

    no_profile_names <- plot(fit, "profiles", archetype_names = NULL, plot = FALSE)
    expect_equal(
        no_profile_names[["archetype"]],
        rep(rownames(coordinates(fit)), times = ncol(coordinates(fit)))
    )
})

test_that("plot.archetypes can subset samples in observation-level plots", {
    fit <- manual_fit()
    sample_names <- paste0("s", seq_len(nrow(fit[["compositions"]])))
    rownames(fit[["compositions"]]) <- sample_names
    rownames(fit[["data"]]) <- sample_names

    local_test_pdf()

    comp <- plot(fit, "compositions", subset = c(1L, 3L), legend = FALSE)
    expect_equal(sort(levels(comp[["sample"]])), sample_names[c(1L, 3L)])
    simplex <- plot(fit, "ternary", subset = c("s2", "s4"))
    expect_equal(rownames(simplex), c("s2", "s4"))
    coords <- plot(fit, "coordinates", subset = c(TRUE, FALSE, TRUE, FALSE))
    expect_equal(coords[["name"]][!coords[["archetype"]]], sample_names[c(1L, 3L)])
    pca <- plot(fit, "coordinates", projection = "pca", subset = c("s1", "s4"))
    expect_equal(sum(!pca[["archetype"]]), 2L)

    expect_error(
        plot(fit, "coordinates", subset = c(TRUE, FALSE)),
        "Logical `subset`"
    )
    expect_error(
        plot(fit, "ternary", subset = "missing"),
        "Some `subset`"
    )
})
