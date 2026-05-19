test_that("internal helper names follow package conventions", {
    r_files <- list.files("../../R", pattern = "[.]R$", full.names = TRUE)
    if (!length(r_files))
        r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)

    lines <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
    defs <- sub("^([.]?[A-Za-z][A-Za-z0-9_.]*)[[:space:]]*(<-|=)[[:space:]]*function.*", "\\1", lines)
    defs <- defs[defs != lines]
    defs <- unique(defs)

    exports <- c(
        "aa_init", "anames", "anames<-", "archetypes_directional",
        "archetypes_kernel_pgd", "archetypes_nnls", "archetypes_paa",
        "archetypes_pgd", "best", "compositions", "consistency", "coordinates",
        "directional_archetypes", "fit_simplex", "kernel_archetypes",
        "onehot", "plot_archetypes_compositions", "plot_archetypes_coordinates",
        "plot_archetypes_loss", "plot_archetypes_profiles", "proj_l1",
        "proj_simplex", "run_aa", "step_archetypes", "tune_archetypes"
    )

    allowed_bare <- c(
        exports,
        "archetypes",
        "uniform_archetypes", "furthest_first", "kmeans_pp", "furthest_sum",
        "aa_pp", "hull_outmost", "dirichlet"
    )

    bare <- defs[!startsWith(defs, ".")]
    bare_internal <- setdiff(bare, allowed_bare)
    bare_internal <- bare_internal[!grepl("^[A-Za-z][A-Za-z0-9_]*[.].+", bare_internal)]
    bare_internal <- bare_internal[!grepl("^(is|has|grad)_[A-Za-z0-9_]+$", bare_internal)]

    dotted_internal <- defs[startsWith(defs, ".")]
    dotted_internal <- dotted_internal[!startsWith(dotted_internal, ".aa_")]

    expect_equal(sort(bare_internal), character())
    expect_equal(sort(dotted_internal), character())

    bad_loss_names <- defs[grepl("_loss_terms$|_terms$|_update_loss$", defs, perl = TRUE)]

    expect_equal(sort(bad_loss_names), character())
})
