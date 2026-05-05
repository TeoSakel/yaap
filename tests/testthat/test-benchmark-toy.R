test_that("toy data benchmarks do not regress", {
    skip_if_not(
        identical(Sys.getenv("YAAAP_BENCHMARKS"), "true"),
        "toy benchmarks are opt-in; set YAAAP_BENCHMARKS=true to run"
    )
    skip_if_not_installed("bench")

    current <- run_toy_benchmarks()
    baseline <- read_toy_benchmark_baseline()
    comparison <- compare_toy_benchmarks(current, baseline)

    expect_true(
        all(comparison[["passed"]]),
        info = format_toy_benchmark_failures(comparison)
    )
})
