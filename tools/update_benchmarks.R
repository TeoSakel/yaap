#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg) > 0L) {
    sub("^--file=", "", file_arg[[1L]])
} else {
    "tools/update_benchmarks.R"
}
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
setwd(root)

if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
} else {
    library(yaap)
}

source(file.path(root, "tests", "testthat", "helper-toy-data.R"))
source(file.path(root, "tests", "testthat", "helper-benchmark-toy.R"))

if (!requireNamespace("bench", quietly = TRUE))
    stop("Package 'bench' is required to update benchmarks.", call. = FALSE)

path <- file.path(root, "tests", "testthat", "benchmark-toy-baseline.csv")
baseline <- write_toy_benchmark_baseline(path)
print(baseline)
