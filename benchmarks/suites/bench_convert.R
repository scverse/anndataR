# =============================================================================
# Benchmark suite: conversions
# =============================================================================
# Benchmarks backend-to-backend conversions (HDF5↔InMemory) and
# format conversions (InMemory↔SCE, InMemory↔Seurat).
# =============================================================================

bench_convert <- function(h5ad_paths, iterations, x_types) {
  results <- list()

  # --- Backend conversions (per X type) ---
  for (xt in x_types) {
    path <- h5ad_paths[[xt]]

    # HDF5 → InMemory
    env <- new.env(parent = globalenv())
    env$.ad <- read_h5ad(path, as = "HDF5AnnData")

    results <- c(
      results,
      run_one_benchmark(
        name = paste0("convert_HDF5_to_InMemory_", xt),
        expr = quote(.ad$as_InMemoryAnnData()),
        iterations = iterations,
        env = env
      )
    )
    env$.ad$close()

    # InMemory → HDF5
    env2 <- new.env(parent = globalenv())
    env2$.ad <- read_h5ad(path, as = "InMemoryAnnData")

    results <- c(
      results,
      run_one_benchmark(
        name = paste0("convert_InMemory_to_HDF5_", xt),
        expr = quote({
          .tmp <- tempfile(fileext = ".h5ad")
          .result <- .ad$as_HDF5AnnData(.tmp)
          .result$close()
          unlink(.tmp)
        }),
        iterations = iterations,
        env = env2
      )
    )
  }

  # --- Format conversions (using float_csparse as representative) ---
  path <- h5ad_paths[["float_csparse"]]
  ad <- read_h5ad(path, as = "InMemoryAnnData")

  # InMemory → SingleCellExperiment
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    env_sce <- new.env(parent = globalenv())
    env_sce$.ad <- ad

    results <- c(
      results,
      run_one_benchmark(
        name = "convert_InMemory_to_SCE",
        expr = quote(as_SingleCellExperiment(.ad)),
        iterations = iterations,
        env = env_sce
      )
    )

    # SCE → InMemory
    sce <- as_SingleCellExperiment(ad)
    env_sce2 <- new.env(parent = globalenv())
    env_sce2$.sce <- sce

    results <- c(
      results,
      run_one_benchmark(
        name = "convert_SCE_to_InMemory",
        expr = quote(as_AnnData(.sce)),
        iterations = iterations,
        env = env_sce2
      )
    )
  } else {
    message("  [SKIP] SCE conversion: SingleCellExperiment not installed")
  }

  # InMemory → Seurat
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    env_seu <- new.env(parent = globalenv())
    env_seu$.ad <- ad

    results <- c(
      results,
      run_one_benchmark(
        name = "convert_InMemory_to_Seurat",
        expr = quote(suppressWarnings(as_Seurat(.ad))),
        iterations = iterations,
        env = env_seu
      )
    )

    # Seurat → InMemory
    seurat <- suppressWarnings(as_Seurat(ad))
    env_seu2 <- new.env(parent = globalenv())
    env_seu2$.seurat <- seurat

    results <- c(
      results,
      run_one_benchmark(
        name = "convert_Seurat_to_InMemory",
        expr = quote(as_AnnData(.seurat)),
        iterations = iterations,
        env = env_seu2
      )
    )
  } else {
    message("  [SKIP] Seurat conversion: SeuratObject not installed")
  }

  results
}
