toy_benchmark_baseline_path <- function() {
    candidates <- c(
        "benchmark-toy-baseline.csv",
        file.path("tests", "testthat", "benchmark-toy-baseline.csv")
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0)
        return(existing[[1L]])

    candidates[[2L]]
}

toy_benchmark_iterations <- function() {
    20L
}

toy_benchmark_slowdown_threshold <- function() {
    1.10
}

toy_benchmark_data <- function(N = 2000L, noise = 0) {
    K <- 3L
    M <- 2L

    set.seed(20250804)
    A <- c(
        cos(0), cos(2 * pi / 3), cos(2 * pi / 3 * 2),
        sin(0), sin(2 * pi / 3), sin(2 * pi / 3 * 2)
    )
    A <- matrix(A, nrow = K, ncol = M)
    S <- -log(matrix(runif(K * N), nrow = N, ncol = K))
    S <- S / rowSums(S)

    S %*% A + noise * rnorm(nrow(S) * M)
}

run_toy_benchmarks <- function(iterations = toy_benchmark_iterations()) {
    X <- toy_benchmark_data()

    set.seed(1)
    invisible(archetypes_pgd(X, K = 3L, max_iter = 20L, tol_r2 = 0.95))
    set.seed(1)
    invisible(archetypes_nnls(X, K = 3L, max_iter = 5L, tol_r2 = 0.95))
    invisible(gc())

    results <- bench::mark(
        pgd = {
            set.seed(1)
            archetypes_pgd(X, K = 3L, max_iter = 20L, tol_r2 = 0.95)
        },
        nnls = {
            set.seed(1)
            archetypes_nnls(X, K = 3L, max_iter = 5L, tol_r2 = 0.95)
        },
        iterations = iterations,
        check = FALSE,
        memory = FALSE,
        filter_gc = FALSE
    )

    data.frame(
        name = as.character(results[["expression"]]),
        median_sec = as.numeric(results[["median"]]),
        iterations = iterations,
        stringsAsFactors = FALSE
    )
}

read_toy_benchmark_baseline <- function(path = toy_benchmark_baseline_path()) {
    if (!file.exists(path))
        stop("Benchmark baseline file not found: ", path, call. = FALSE)

    baseline <- utils::read.csv(path, stringsAsFactors = FALSE)
    required <- c("name", "median_sec", "iterations", "package_version", "updated_at")
    missing <- setdiff(required, names(baseline))
    if (length(missing) > 0)
        stop("Benchmark baseline is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

    baseline
}

compare_toy_benchmarks <- function(current,
                                   baseline,
                                   threshold = toy_benchmark_slowdown_threshold()) {
    merged <- merge(
        current,
        baseline[c("name", "median_sec")],
        by = "name",
        suffixes = c("_current", "_baseline"),
        all.x = TRUE
    )

    if (anyNA(merged[["median_sec_baseline"]])) {
        missing <- merged[["name"]][is.na(merged[["median_sec_baseline"]])]
        stop("Benchmark baseline is missing rows: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    merged[["ratio"]] <- merged[["median_sec_current"]] / merged[["median_sec_baseline"]]
    merged[["passed"]] <- merged[["ratio"]] <= threshold
    merged
}

format_toy_benchmark_failures <- function(comparison,
                                          threshold = toy_benchmark_slowdown_threshold()) {
    failures <- comparison[!comparison[["passed"]], , drop = FALSE]
    if (nrow(failures) == 0L)
        return("")

    lines <- sprintf(
        "%s: current %.3fs, baseline %.3fs, ratio %.2fx",
        failures[["name"]],
        failures[["median_sec_current"]],
        failures[["median_sec_baseline"]],
        failures[["ratio"]]
    )

    paste(
        sprintf("Toy benchmarks exceeded the %.0f%% slowdown threshold.", (threshold - 1) * 100),
        paste(lines, collapse = "\n"),
        sep = "\n"
    )
}

write_toy_benchmark_baseline <- function(path = toy_benchmark_baseline_path(),
                                         iterations = toy_benchmark_iterations()) {
    baseline <- run_toy_benchmarks(iterations = iterations)
    baseline[["package_version"]] <- as.character(utils::packageVersion("yaap"))
    baseline[["updated_at"]] <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(baseline, path, row.names = FALSE)
    baseline
}
