# =============================================================================
# Benchmark suite: read_h5ad
# =============================================================================
# Benchmarks reading H5AD files into InMemoryAnnData and HDF5AnnData
# across different X matrix types.
# =============================================================================

bench_read <- function(h5ad_paths, iterations, x_types, zarr_paths) {
  results <- list()

  for (xt in x_types) {
    path <- h5ad_paths[[xt]]

    # Read → InMemoryAnnData
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("read_InMemory_", xt),
        expr = quote(read_h5ad(.path, as = "InMemoryAnnData")),
        setup = bquote(.path <- .(path)),
        iterations = iterations
      )
    )

    # Read → HDF5AnnData
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("read_HDF5_", xt),
        expr = quote({
          ad <- read_h5ad(.path, as = "HDF5AnnData")
          ad$close()
        }),
        setup = bquote(.path <- .(path)),
        iterations = iterations
      )
    )
  }

  # Read from Zarr store
  for (xt in x_types) {
    path <- zarr_paths[[xt]]

    # Read Zarr → InMemoryAnnData
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("read_zarr_InMemory_", xt),
        expr = quote(read_zarr(.path, as = "InMemoryAnnData")),
        setup = bquote(.path <- .(path)),
        iterations = iterations
      )
    )

    # Open Zarr lazily → ZarrAnnData
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("read_zarr_Zarr_", xt),
        expr = quote(read_zarr(.path, as = "ZarrAnnData")),
        setup = bquote(.path <- .(path)),
        iterations = iterations
      )
    )
  }

  results
}
