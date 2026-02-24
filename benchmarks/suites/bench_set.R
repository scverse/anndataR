# =============================================================================
# Benchmark suite: slot setters
# =============================================================================
# Benchmarks setting every AnnData slot on both InMemory and HDF5 backends.
# =============================================================================

bench_set <- function(h5ad_paths, iterations) {
  results <- list()
  path <- h5ad_paths[["float_csparse"]]

  slots <- c(
    "X",
    "obs",
    "var",
    "obs_names",
    "var_names",
    "layers",
    "obsm",
    "varm",
    "obsp",
    "varp",
    "uns"
  )

  for (backend in c("InMemoryAnnData", "HDF5AnnData")) {
    short <- if (backend == "InMemoryAnnData") "InMemory" else "HDF5"

    for (slot in slots) {
      # For HDF5, we need a fresh writable copy for each slot
      if (backend == "HDF5AnnData") {
        tmp <- tempfile(fileext = ".h5ad")
        file.copy(path, tmp)
        ad <- read_h5ad(tmp, as = "HDF5AnnData", mode = "r+")
      } else {
        ad <- read_h5ad(path, as = "InMemoryAnnData")
      }

      # Read the current value to use for setting
      val <- ad[[slot]]

      env <- new.env(parent = globalenv())
      env$.ad <- ad
      env$.val <- val
      env$.slot <- slot

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("set_", short, "_", slot),
          expr = quote(.ad[[.slot]] <- .val),
          iterations = iterations,
          env = env
        )
      )

      if (backend == "HDF5AnnData") {
        ad$close()
        unlink(tmp)
      }
    }
  }

  results
}
