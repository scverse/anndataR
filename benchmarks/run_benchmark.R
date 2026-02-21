#!/usr/bin/env Rscript

# =============================================================================
# anndataR Performance Benchmark
# =============================================================================
#
# Benchmarks read_h5ad() performance across different X matrix types using
# dummy-anndata to generate test datasets. Compares dense, CSC sparse, and
# CSR sparse matrices.
#
# Usage:
#   Rscript benchmarks/run_benchmark.R
#   Rscript benchmarks/run_benchmark.R --n_obs 100000 --n_vars 30000
#   Rscript benchmarks/run_benchmark.R --source installed

library(optparse)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--n_obs",
    type = "integer",
    default = 10000L,
    help = "Number of observations [default %default]"
  ),
  make_option(
    "--n_vars",
    type = "integer",
    default = 5000L,
    help = "Number of variables [default %default]"
  ),
  make_option(
    "--times",
    type = "integer",
    default = 3L,
    help = "Repetitions per benchmark [default %default]"
  ),
  make_option(
    "--source",
    type = "character",
    default = "local",
    help = "'local' (devtools::load_all) or 'installed' (library) [default %default]"
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opts <- parse_args(parser)

n_obs <- opts$n_obs
n_vars <- opts$n_vars
times <- opts$times

if (opts$source == "local") {
  cat("Loading anndataR from local source (devtools::load_all)...\n")
  devtools::load_all()
} else if (opts$source == "installed") {
  cat("Loading anndataR from installed package (library)...\n")
  library(anndataR)
} else {
  stop("--source must be 'local' or 'installed', got: ", opts$source)
}

cat(sprintf(
  "\n=== anndataR performance benchmark ===\nn_obs = %d, n_vars = %d, times = %d\n\n",
  n_obs,
  n_vars,
  times
))

# ---------------------------------------------------------------------------
# Generate test H5AD files using dummy-anndata
# ---------------------------------------------------------------------------
if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("reticulate is required")
}

da <- reticulate::import("dummy_anndata", convert = FALSE)

# X types to benchmark: dense, CSC sparse, CSR sparse (float & integer)
x_types <- c(
  "float_matrix",
  "float_csparse",
  "float_rsparse",
  "integer_matrix",
  "integer_csparse",
  "integer_rsparse"
)

cache_dir <- file.path(Sys.getenv("HOME"), ".cache", "anndataR_bench")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

h5ad_paths <- list()
for (x_type in x_types) {
  cache_key <- sprintf("bench_%d_%d_%s.h5ad", n_obs, n_vars, x_type)
  path <- file.path(cache_dir, cache_key)

  if (file.exists(path)) {
    cat(sprintf(
      "  [cached] %s: %s (%.1f MB)\n",
      x_type,
      path,
      file.size(path) / 1e6
    ))
  } else {
    cat(sprintf("  Generating %s dataset... ", x_type))
    gen_time <- system.time({
      adata <- da$generate_dataset(
        n_obs = as.integer(n_obs),
        n_vars = as.integer(n_vars),
        x_type = x_type,
        layer_types = list(),
        obs_types = list("categorical", "dense_array"),
        var_types = list("string_array", "boolean_array"),
        obsm_types = list(),
        varm_types = list(),
        obsp_types = list(),
        varp_types = list(),
        uns_types = list(),
        nested_uns_types = list()
      )
      adata$write_h5ad(path)
    })
    cat(sprintf(
      "done (%.1fs, %.1f MB)\n",
      gen_time[["elapsed"]],
      file.size(path) / 1e6
    ))
  }
  h5ad_paths[[x_type]] <- path
}
cat("\n")

# ---------------------------------------------------------------------------
# Benchmark helper
# ---------------------------------------------------------------------------
bench_one <- function(label, expr, times) {
  timings <- numeric(times)
  for (i in seq_len(times)) {
    gc(FALSE)
    timings[i] <- system.time(eval(expr))["elapsed"]
  }
  cat(sprintf(
    "  %-55s  median=%.3fs  [%s]\n",
    label,
    median(timings),
    paste(sprintf("%.3f", timings), collapse = ", ")
  ))
  invisible(median(timings))
}

# ---------------------------------------------------------------------------
# 1. read_h5ad as InMemoryAnnData by X type
# ---------------------------------------------------------------------------
cat("--- read_h5ad(as='InMemoryAnnData') by X type -------------------------\n")
for (x_type in x_types) {
  path <- h5ad_paths[[x_type]]
  bench_one(
    sprintf("X=%s", x_type),
    bquote(read_h5ad(.(path), as = "InMemoryAnnData")),
    times
  )
}
cat("\n")

# ---------------------------------------------------------------------------
# 2. read_h5ad as HDF5AnnData (lazy open) by X type
# ---------------------------------------------------------------------------
cat("--- read_h5ad(as='HDF5AnnData') by X type -----------------------------\n")
for (x_type in x_types) {
  path <- h5ad_paths[[x_type]]
  bench_one(
    sprintf("X=%s", x_type),
    bquote({
      ad <- read_h5ad(.(path), as = "HDF5AnnData")
      ad$close()
    }),
    times
  )
}
cat("\n")

# ---------------------------------------------------------------------------
# 3. Per-slot access on HDF5AnnData for each X type
# ---------------------------------------------------------------------------
cat("--- HDF5AnnData$X read by X type --------------------------------------\n")
for (x_type in x_types) {
  path <- h5ad_paths[[x_type]]
  hdf5_ad <- read_h5ad(path, as = "HDF5AnnData")
  bench_one(
    sprintf("hdf5_ad$X  (X=%s)", x_type),
    quote(hdf5_ad$X),
    times
  )
  hdf5_ad$close()
}
cat("\n")

# ---------------------------------------------------------------------------
# 4. Per-slot access for metadata (using float_rsparse file)
# ---------------------------------------------------------------------------
cat("--- Per-slot access (HDF5AnnData, X=float_rsparse) --------------------\n")
hdf5_ad <- read_h5ad(h5ad_paths[["float_rsparse"]], as = "HDF5AnnData")

bench_one("hdf5_ad$obs_names", quote(hdf5_ad$obs_names), times)
bench_one("hdf5_ad$var_names", quote(hdf5_ad$var_names), times)
bench_one("hdf5_ad$n_obs()", quote(hdf5_ad$n_obs()), times)
bench_one("hdf5_ad$n_vars()", quote(hdf5_ad$n_vars()), times)
bench_one("hdf5_ad$obs", quote(hdf5_ad$obs), times)
bench_one("hdf5_ad$var", quote(hdf5_ad$var), times)
bench_one("hdf5_ad$X", quote(hdf5_ad$X), times)

hdf5_ad$close()
cat("\n")

# ---------------------------------------------------------------------------
# 5. Rprof the hot path (InMemory read with float_rsparse)
# ---------------------------------------------------------------------------
cat("--- Rprof of read_h5ad(as='InMemoryAnnData', X=float_rsparse) ---------\n")
prof_path <- h5ad_paths[["float_rsparse"]]
prof_file <- tempfile(fileext = ".prof")
Rprof(prof_file, interval = 0.005)
for (i in seq_len(times)) {
  ad <- read_h5ad(prof_path, as = "InMemoryAnnData")
}
Rprof(NULL)

prof <- summaryRprof(prof_file)
cat("\nTop 15 functions by self.time:\n")
print(head(prof$by.self, 15))
cat("\nTop 15 functions by total.time:\n")
print(head(prof$by.total, 15))

# ---------------------------------------------------------------------------
# 6. Repeated name access
# ---------------------------------------------------------------------------
cat(
  "\n--- Repeated name access (HDF5AnnData, X=float_rsparse) ---------------\n"
)
hdf5_ad2 <- read_h5ad(h5ad_paths[["float_rsparse"]], as = "HDF5AnnData")

bench_one(
  "obs_names x 10 calls",
  quote(
    for (j in 1:10) {
      hdf5_ad2$obs_names
    }
  ),
  times
)
bench_one(
  "var_names x 10 calls",
  quote(
    for (j in 1:10) {
      hdf5_ad2$var_names
    }
  ),
  times
)

hdf5_ad2$close()

cat(sprintf("\nCache dir: %s\n", cache_dir))
cat("\n=== Benchmark complete ===\n")
