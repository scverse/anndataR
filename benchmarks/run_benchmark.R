#!/usr/bin/env Rscript

# =============================================================================
# anndataR Performance Benchmark
# =============================================================================
#
# This script creates a simple sparse H5AD file using Python anndata (via
# reticulate) and then benchmarks the key code-paths that matter for the paper
# benchmark:
#
#   1. read_h5ad(…, as = "InMemoryAnnData")   – the main user-facing path
#   2. read_h5ad(…, as = "HDF5AnnData")       – just opening the file
#   3. HDF5AnnData → InMemoryAnnData          – the conversion step
#   4. Individual slot access on HDF5AnnData   – profiling where time goes
#
# It also profiles the hot-path with Rprof so we can see exactly which
# internal functions dominate.
#
# Usage:
#   Rscript benchmarks/run_benchmark.R
#   Rscript benchmarks/run_benchmark.R --n_obs 50000
#   Rscript benchmarks/run_benchmark.R --source installed
#   Rscript benchmarks/run_benchmark.R --source local

library(optparse)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--n_obs",
    type = "integer",
    default = 50000L,
    help = "Number of observations (cells) [default %default]"
  ),
  make_option(
    "--n_vars",
    type = "integer",
    default = 20000L,
    help = "Number of variables (genes) [default %default]"
  ),
  make_option(
    "--density",
    type = "double",
    default = 0.05,
    help = "Sparse matrix density [default %default]"
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random seed for dataset generation [default %default]"
  ),
  make_option(
    "--source",
    type = "character",
    default = "local",
    help = paste(
      "Package source: 'local' (devtools::load_all) or",
      "'installed' (library) [default %default]"
    )
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opts <- parse_args(parser)

n_obs <- opts$n_obs
n_vars <- opts$n_vars
density <- opts$density
seed <- opts$seed

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
  "\n=== anndataR performance benchmark ===\nn_obs = %d, n_vars = %d, density = %.2f, seed = %d\n\n",
  n_obs,
  n_vars,
  density,
  seed
))

# ---------------------------------------------------------------------------
# Generate test H5AD file using scipy/anndata
# ---------------------------------------------------------------------------
generate_h5ad <- function(path, n_obs, n_vars, density, seed = 42L) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to generate test data")
  }
  sp <- reticulate::import("scipy.sparse")
  ad <- reticulate::import("anndata")
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")

  np$random$seed(as.integer(seed))

  X <- sp$random(
    as.integer(n_obs),
    as.integer(n_vars),
    density = density,
    format = "csr",
    dtype = "float32",
    random_state = as.integer(seed)
  )

  obs_names <- paste0("Cell", sprintf("%06d", seq_len(n_obs) - 1L))
  var_names <- paste0("Gene", sprintf("%05d", seq_len(n_vars) - 1L))

  obs <- pd$DataFrame(
    list(
      cell_type = np$random$choice(
        c("T", "B", "NK", "Mono"),
        size = as.integer(n_obs)
      ),
      n_counts = np$random$rand(as.integer(n_obs))
    ),
    index = obs_names
  )

  var <- pd$DataFrame(
    list(
      gene_name = var_names,
      highly_variable = np$random$choice(
        c(TRUE, FALSE),
        size = as.integer(n_vars)
      )
    ),
    index = var_names
  )

  adata <- ad$AnnData(X = X, obs = obs, var = var)
  adata$write_h5ad(path)
  invisible(path)
}

# Use a deterministic cache path based on parameters so the file is reused
cache_dir <- file.path(Sys.getenv("HOME"), ".cache", "anndataR_bench")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache_key <- sprintf("bench_%d_%d_%.4f_%d.h5ad", n_obs, n_vars, density, seed)
h5ad_path <- file.path(cache_dir, cache_key)

if (file.exists(h5ad_path)) {
  cat(sprintf("Using cached dataset: %s\n", h5ad_path))
  cat(sprintf("File size: %.1f MB\n\n", file.size(h5ad_path) / 1e6))
} else {
  cat("Generating test dataset... ")
  gen_time <- system.time(
    generate_h5ad(h5ad_path, n_obs, n_vars, density, seed)
  )
  cat(sprintf("done (%.1fs)\n", gen_time[["elapsed"]]))
  cat(sprintf("File size: %.1f MB\n\n", file.size(h5ad_path) / 1e6))
}

# ---------------------------------------------------------------------------
# Benchmark helper
# ---------------------------------------------------------------------------
bench_one <- function(label, expr, times = 3L) {
  timings <- numeric(times)
  for (i in seq_len(times)) {
    gc(FALSE)
    timings[i] <- system.time(eval(expr))["elapsed"]
  }
  cat(sprintf(
    "  %-50s  median=%.3fs  [%s]\n",
    label,
    median(timings),
    paste(sprintf("%.3f", timings), collapse = ", ")
  ))
  invisible(median(timings))
}

# ---------------------------------------------------------------------------
# 1. Top-level read_h5ad benchmarks
# ---------------------------------------------------------------------------
cat(
  "--- Top-level read_h5ad -------------------------------------------------\n"
)
bench_one(
  "read_h5ad(as='HDF5AnnData')",
  quote({
    ad <- read_h5ad(h5ad_path, as = "HDF5AnnData")
    ad$close()
  })
)

bench_one(
  "read_h5ad(as='InMemoryAnnData')",
  quote({
    ad <- read_h5ad(h5ad_path, as = "InMemoryAnnData")
  })
)

bench_one(
  "read_h5ad + sum(X@x)  [force materialisation]",
  quote({
    ad <- read_h5ad(h5ad_path, as = "InMemoryAnnData")
    sum(ad$X@x)
  })
)

cat("\n")

# ---------------------------------------------------------------------------
# 2. Per-slot access on HDF5AnnData
# ---------------------------------------------------------------------------
cat(
  "--- Per-slot access (HDF5AnnData) --------------------------------------\n"
)
hdf5_ad <- read_h5ad(h5ad_path, as = "HDF5AnnData")

bench_one("hdf5_ad$obs_names", quote(hdf5_ad$obs_names))
bench_one("hdf5_ad$var_names", quote(hdf5_ad$var_names))
bench_one("hdf5_ad$n_obs()", quote(hdf5_ad$n_obs()))
bench_one("hdf5_ad$n_vars()", quote(hdf5_ad$n_vars()))
bench_one("hdf5_ad$obs", quote(hdf5_ad$obs))
bench_one("hdf5_ad$var", quote(hdf5_ad$var))
bench_one("hdf5_ad$X", quote(hdf5_ad$X))
bench_one("hdf5_ad$uns", quote(hdf5_ad$uns))

cat("\n")

# ---------------------------------------------------------------------------
# 3. Conversion HDF5 -> InMemory
# ---------------------------------------------------------------------------
cat(
  "--- HDF5AnnData -> InMemoryAnnData --------------------------------------\n"
)

bench_one(
  "as_InMemoryAnnData (full)",
  quote({
    hdf5_tmp <- read_h5ad(h5ad_path, as = "HDF5AnnData")
    hdf5_tmp$as_InMemoryAnnData()
    hdf5_tmp$close()
  })
)

cat("\n")

# ---------------------------------------------------------------------------
# 4. Profiling the hot path
# ---------------------------------------------------------------------------
cat(
  "--- Rprof of read_h5ad(as='InMemoryAnnData') ---------------------------\n"
)
prof_file <- tempfile(fileext = ".prof")
Rprof(prof_file, interval = 0.005)
for (i in 1:3) {
  ad <- read_h5ad(h5ad_path, as = "InMemoryAnnData")
}
Rprof(NULL)

prof <- summaryRprof(prof_file)
cat("\nTop 20 functions by self.time:\n")
top20 <- head(prof$by.self, 20)
print(top20)

cat("\nTop 20 functions by total.time:\n")
top20_total <- head(prof$by.total, 20)
print(top20_total)

# ---------------------------------------------------------------------------
# 5. Detailed HDF5 internal timings
# ---------------------------------------------------------------------------
cat(
  "\n--- Detailed internal timings -------------------------------------------\n"
)
# Time just reading obs_names repeatedly (simulates what happens internally)
hdf5_ad2 <- read_h5ad(h5ad_path, as = "HDF5AnnData")

bench_one(
  "obs_names x 10 calls",
  quote({
    for (j in 1:10) {
      hdf5_ad2$obs_names
    }
  })
)

bench_one(
  "var_names x 10 calls",
  quote({
    for (j in 1:10) {
      hdf5_ad2$var_names
    }
  })
)

# Time raw rhdf5 read for comparison
bench_one(
  "raw rhdf5::h5read(X) [reference]",
  quote({
    file <- rhdf5::H5Fopen(h5ad_path, flags = "H5F_ACC_RDONLY")
    grp <- rhdf5::H5Gopen(file, "X")
    data <- as.vector(grp$data)
    indices <- as.vector(grp$indices)
    indptr <- as.vector(grp$indptr)
    rhdf5::H5Gclose(grp)
    rhdf5::H5Fclose(file)
  })
)

bench_one(
  "raw sparseMatrix construction [reference]",
  quote({
    file <- rhdf5::H5Fopen(h5ad_path, flags = "H5F_ACC_RDONLY")
    grp <- rhdf5::H5Gopen(file, "X")
    data <- as.vector(grp$data)
    indices <- as.vector(grp$indices)
    indptr <- as.vector(grp$indptr)
    attrs <- rhdf5::h5readAttributes(file, "X")
    shape <- as.vector(attrs[["shape"]])
    rhdf5::H5Gclose(grp)
    rhdf5::H5Fclose(file)
    mtx <- Matrix::sparseMatrix(
      j = indices,
      p = indptr,
      x = data,
      dims = shape,
      repr = "R",
      index1 = FALSE
    )
  })
)

hdf5_ad2$close()
hdf5_ad$close()

cat(sprintf("\nCached dataset: %s\n", h5ad_path))
cat("\n=== Benchmark complete ===\n")
