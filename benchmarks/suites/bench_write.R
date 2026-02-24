# =============================================================================
# Benchmark suite: write_h5ad
# =============================================================================
# Benchmarks writing AnnData objects to H5AD files from both backends,
# with different compression settings and X matrix types.
# =============================================================================

bench_write <- function(h5ad_paths, iterations, x_types) {
  results <- list()

  compressions <- c("none", "gzip")

  for (xt in x_types) {
    path <- h5ad_paths[[xt]]

    for (compression in compressions) {
      # Write from InMemoryAnnData
      env <- new.env(parent = globalenv())
      env$.ad <- read_h5ad(path, as = "InMemoryAnnData")
      env$.compression <- compression

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("write_InMemory_", xt, "_", compression),
          expr = quote({
            .tmp <- tempfile(fileext = ".h5ad")
            .ad$write_h5ad(.tmp, compression = .compression)
            unlink(.tmp)
          }),
          iterations = iterations,
          env = env
        )
      )

      # Write from HDF5AnnData
      env2 <- new.env(parent = globalenv())
      env2$.ad <- read_h5ad(path, as = "HDF5AnnData")
      env2$.compression <- compression

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("write_HDF5_", xt, "_", compression),
          expr = quote({
            .tmp <- tempfile(fileext = ".h5ad")
            .ad$write_h5ad(.tmp, compression = .compression)
            unlink(.tmp)
          }),
          iterations = iterations,
          env = env2
        )
      )

      # Close HDF5 handle
      env2$.ad$close()
    }
  }

  results
}
