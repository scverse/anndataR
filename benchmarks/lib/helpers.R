# =============================================================================
# Benchmark helpers: data generation, bench→BMF conversion, JSON output
# =============================================================================

# Load anndataR with internal functions available.
# In development: devtools::load_all() exposes unexported functions.
# In CI: the package is installed, so we use ::: for unexported functions.
if (
  file.exists("DESCRIPTION") &&
    grepl("anndataR", readLines("DESCRIPTION", n = 1))
) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(anndataR)
}
library(Matrix)

# ---------------------------------------------------------------------------
# Python dependency check
# ---------------------------------------------------------------------------

#' Check that reticulate, anndata, and dummy_anndata are available
check_bench_python_deps <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "reticulate is required for benchmarks.\n",
      "Install with: install.packages('reticulate')"
    )
  }
  if (!reticulate::py_module_available("anndata")) {
    stop(
      "Python 'anndata' module is required for benchmarks.\n",
      "Install with: reticulate::py_install('anndata')"
    )
  }
  if (!reticulate::py_module_available("dummy_anndata")) {
    stop(
      "Python 'dummy_anndata' module is required for benchmarks.\n",
      "Install with: reticulate::py_install('dummy-anndata')"
    )
  }
}

# ---------------------------------------------------------------------------
# Data generation (via Python dummy_anndata)
# ---------------------------------------------------------------------------

#' Generate and cache an H5AD file with all slots populated
#'
#' Uses Python's dummy_anndata to generate and write the H5AD file,
#' ensuring canonical encoding independent of anndataR's own writer.
#'
#' @param x_type Matrix type for X (e.g. "float_csparse").
#'   Must be a key in `dummy_anndata.matrix_generators`.
#' @param n_obs Number of observations
#' @param n_vars Number of variables
#' @param cache_dir Directory to cache generated files
#' @return Path to the generated H5AD file
generate_bench_h5ad <- function(x_type, n_obs, n_vars, cache_dir) {
  path <- file.path(cache_dir, paste0("bench_", x_type, ".h5ad"))
  if (file.exists(path)) {
    return(path)
  }

  da <- reticulate::import("dummy_anndata", convert = FALSE)

  adata_py <- da$generate_dataset(
    n_obs = as.integer(n_obs),
    n_vars = as.integer(n_vars),
    x_type = x_type,
    layer_types = list("float_csparse", "integer_csparse"),
    obs_types = list("categorical", "dense_array", "string_array"),
    var_types = list("string_array", "boolean_array", "dense_array"),
    obsm_types = list("float_matrix", "float_csparse"),
    varm_types = list("float_matrix"),
    obsp_types = list("float_csparse"),
    varp_types = list("float_csparse"),
    uns_types = list("string", "integer", "float"),
    nested_uns_types = list("string", "float")
  )

  adata_py$write_h5ad(path)
  path
}

# ---------------------------------------------------------------------------
# bench::mark → BMF JSON conversion
# ---------------------------------------------------------------------------

#' Convert a bench::mark result into a BMF JSON entry
#'
#' @param bm A bench::mark result (tibble with one row)
#' @param name Benchmark name
#' @return A named list with BMF structure
bench_to_bmf <- function(bm, name) {
  # bench::mark returns bench_time objects; convert to nanoseconds
  median_ns <- as.numeric(bm$median, units = "secs") * 1e9
  min_ns <- as.numeric(bm$min, units = "secs") * 1e9
  # max may not be available with iterations = 1
  max_ns <- if ("max" %in% names(bm)) {
    as.numeric(bm$max, units = "secs") * 1e9
  } else {
    median_ns
  }

  entry <- list(
    latency = list(
      value = median_ns,
      lower_value = min_ns,
      upper_value = max_ns
    )
  )

  # Add memory allocation if available
  mem_bytes <- as.numeric(bm$mem_alloc)
  if (!is.na(mem_bytes)) {
    entry[["memory"]] <- list(value = mem_bytes)
  }

  stats::setNames(list(entry), name)
}

#' Write accumulated BMF results to a JSON file
#'
#' @param results Named list of BMF entries (as returned by bench_to_bmf)
#' @param path Output file path
write_bmf_json <- function(results, path) {
  # After accumulation via c(), results is already a flat named list:
  # list("bench_name1" = list(latency = ..., memory = ...),
  #      "bench_name2" = list(latency = ..., memory = ...), ...)
  # Write it directly as BMF JSON.
  jsonlite::write_json(results, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Safe benchmark runner
# ---------------------------------------------------------------------------

#' Run a single benchmark safely, returning a BMF entry or NULL on failure
#'
#' @param name Benchmark name
#' @param expr Expression to benchmark (quoted)
#' @param setup Expression to run before benchmarking (quoted)
#' @param iterations Number of iterations
#' @param env Environment for evaluation
#' @return A BMF entry (list) or NULL if the benchmark failed
run_one_benchmark <- function(
  name,
  expr,
  setup = NULL,
  iterations = 3L,
  env = parent.frame()
) {
  tryCatch(
    {
      if (!is.null(setup)) {
        eval(setup, envir = env)
      }
      bm <- bench::mark(
        eval(expr, envir = env),
        iterations = iterations,
        check = FALSE,
        filter_gc = FALSE
      )
      bench_to_bmf(bm, name)
    },
    error = function(e) {
      message("  [SKIP] ", name, ": ", conditionMessage(e))
      NULL
    }
  )
}

# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

#' Parse command-line arguments for the benchmark runner
parse_bench_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    n_obs = 2000L,
    n_vars = 1000L,
    iterations = 3L,
    suite = "all",
    output = "benchmarks/results.json"
  )

  i <- 1L
  while (i <= length(args)) {
    switch(
      args[i],
      "--n-obs" = {
        i <- i + 1L
        opts$n_obs <- as.integer(args[i])
      },
      "--n-vars" = {
        i <- i + 1L
        opts$n_vars <- as.integer(args[i])
      },
      "--iterations" = ,
      "-n" = {
        i <- i + 1L
        opts$iterations <- as.integer(args[i])
      },
      "--suite" = ,
      "-s" = {
        i <- i + 1L
        opts$suite <- args[i]
      },
      "--output" = ,
      "-o" = {
        i <- i + 1L
        opts$output <- args[i]
      },
      "--help" = ,
      "-h" = {
        cat(
          "Usage: Rscript benchmarks/run_benchmarks.R [OPTIONS]\n",
          "\n",
          "Options:\n",
          "  --n-obs N        Number of observations (default: 2000)\n",
          "  --n-vars N       Number of variables (default: 1000)\n",
          "  --iterations N   Iterations per benchmark (default: 3)\n",
          "  --suite NAME     Suite to run: all, read, write, get, set,\n",
          "                   convert, subset (default: all)\n",
          "  --output FILE    Output JSON path (default: benchmarks/results.json)\n",
          "  --help           Show this help\n",
          sep = ""
        )
        quit(status = 0)
      },
      {
        warning("Unknown argument: ", args[i])
      }
    )
    i <- i + 1L
  }
  opts
}
