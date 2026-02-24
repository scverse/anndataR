#!/usr/bin/env Rscript

# =============================================================================
# anndataR Benchmark Runner
# =============================================================================
# Runs bench::mark benchmarks and outputs Bencher Metric Format (BMF) JSON.
#
# Usage:
#   Rscript benchmarks/run_benchmarks.R [OPTIONS]
#
# Options:
#   --n-obs N        Number of observations (default: 2000)
#   --n-vars N       Number of variables (default: 1000)
#   --iterations N   Iterations per benchmark (default: 3)
#   --suite NAME     Suite to run: all, read, write, get, set,
#                    convert, subset (default: all)
#   --output FILE    Output JSON path (default: benchmarks/results.json)
#   --help           Show help
# =============================================================================

# Load helpers
source("benchmarks/lib/helpers.R")

# Load suite scripts
source("benchmarks/suites/bench_read.R")
source("benchmarks/suites/bench_write.R")
source("benchmarks/suites/bench_get.R")
source("benchmarks/suites/bench_set.R")
source("benchmarks/suites/bench_convert.R")
source("benchmarks/suites/bench_subset.R")

# Parse CLI arguments
opts <- parse_bench_args()

cat("=== anndataR Benchmark Suite ===\n")
cat(sprintf("  n_obs:       %d\n", opts$n_obs))
cat(sprintf("  n_vars:      %d\n", opts$n_vars))
cat(sprintf("  iterations:  %d\n", opts$iterations))
cat(sprintf("  suite:       %s\n", opts$suite))
cat(sprintf("  output:      %s\n", opts$output))
cat("\n")

# Check Python dependencies (dummy_anndata used for data generation)
check_bench_python_deps()

# ---------------------------------------------------------------------------
# Generate / cache test H5AD files (written by Python anndata)
# ---------------------------------------------------------------------------
# x_types use dummy_anndata naming: float_*, integer_*
x_types <- c(
  "float_matrix",
  "float_csparse",
  "float_rsparse",
  "integer_csparse"
)

cache_dir <- file.path(tempdir(), "anndatar_bench")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating test data (via Python dummy_anndata)...\n")
h5ad_paths <- setNames(
  vapply(
    x_types,
    function(xt) {
      cat(sprintf("  %s... ", xt))
      path <- generate_bench_h5ad(xt, opts$n_obs, opts$n_vars, cache_dir)
      cat("done\n")
      path
    },
    character(1)
  ),
  x_types
)
cat("\n")

# ---------------------------------------------------------------------------
# Run selected suites
# ---------------------------------------------------------------------------
all_results <- list()
suites_to_run <- if (opts$suite == "all") {
  c("read", "write", "get", "set", "convert", "subset")
} else {
  opts$suite
}

for (suite in suites_to_run) {
  cat(sprintf("Running suite: %s\n", suite))

  suite_results <- switch(
    suite,
    read = bench_read(h5ad_paths, opts$iterations, x_types),
    write = bench_write(h5ad_paths, opts$iterations, x_types),
    get = bench_get(h5ad_paths, opts$iterations),
    set = bench_set(h5ad_paths, opts$iterations),
    convert = bench_convert(h5ad_paths, opts$iterations, x_types),
    subset = bench_subset(h5ad_paths, opts$iterations),
    {
      warning("Unknown suite: ", suite)
      list()
    }
  )

  n_ok <- sum(!vapply(suite_results, is.null, logical(1)))
  cat(sprintf("  %d benchmarks completed\n\n", n_ok))
  all_results <- c(all_results, suite_results)
}

# ---------------------------------------------------------------------------
# Write BMF JSON output
# ---------------------------------------------------------------------------
# Filter out NULLs
all_results <- Filter(Negate(is.null), all_results)

cat(sprintf("Total benchmarks: %d\n", length(all_results)))

# Ensure output directory exists
dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
write_bmf_json(all_results, opts$output)

cat(sprintf("Results written to: %s\n", opts$output))
