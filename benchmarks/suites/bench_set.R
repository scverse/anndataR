# =============================================================================
# Benchmark suite: slot setters
# =============================================================================
# Benchmarks setting every AnnData slot on both InMemory and HDF5 backends.
# =============================================================================

bench_set <- function(h5ad_paths, iterations, zarr_paths) {
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

  for (backend in c("InMemoryAnnData", "HDF5AnnData", "ZarrAnnData")) {
    short <- switch(
      backend,
      InMemoryAnnData = "InMemory",
      HDF5AnnData = "HDF5",
      ZarrAnnData = "Zarr"
    )

    for (slot in slots) {
      # Each backend needs a fresh writable instance per slot
      if (backend == "ZarrAnnData") {
        # Copy Zarr store directory so each slot gets a fresh writable copy
        zarr_path <- zarr_paths[["float_csparse"]]
        tmp_parent <- tempfile()
        dir.create(tmp_parent, recursive = TRUE)
        file.copy(zarr_path, tmp_parent, recursive = TRUE)
        tmp <- file.path(tmp_parent, basename(zarr_path))
        ad <- read_zarr(tmp, as = "ZarrAnnData", mode = "r+")
      } else if (backend == "HDF5AnnData") {
        tmp <- tempfile(fileext = ".h5ad")
        file.copy(path, tmp)
        ad <- suppressWarnings(
          read_h5ad(tmp, as = "HDF5AnnData", mode = "r+")
        )
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

      if (backend == "ZarrAnnData") {
        unlink(tmp, recursive = TRUE)
      } else if (backend == "HDF5AnnData") {
        ad$close()
        unlink(tmp)
      }
    }
  }

  results
}
