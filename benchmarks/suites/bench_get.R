# =============================================================================
# Benchmark suite: slot getters
# =============================================================================
# Benchmarks accessing every AnnData slot on both InMemory and HDF5 backends.
# Also benchmarks all *_keys() methods, n_obs(), n_vars(), shape(), and
# S3 methods (dim, nrow, ncol, dimnames, rownames, colnames).
# =============================================================================

.bench_slots <- c(
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

.bench_key_methods <- c(
  "obs_keys",
  "var_keys",
  "layers_keys",
  "obsm_keys",
  "varm_keys",
  "obsp_keys",
  "varp_keys",
  "uns_keys"
)

.bench_dim_methods <- c("n_obs", "n_vars", "shape")

.bench_s3_methods <- list(
  dim = quote(dim(.ad)),
  nrow = quote(nrow(.ad)),
  ncol = quote(ncol(.ad)),
  dimnames = quote(dimnames(.ad)),
  rownames = quote(rownames(.ad)),
  colnames = quote(colnames(.ad))
)

bench_get <- function(h5ad_paths, iterations) {
  results <- list()
  path <- h5ad_paths[["float_csparse"]]

  for (backend in c("InMemoryAnnData", "HDF5AnnData")) {
    short <- if (backend == "InMemoryAnnData") "InMemory" else "HDF5"

    # Open the AnnData
    ad <- read_h5ad(path, as = backend)

    # --- Slot getters ---
    for (slot in .bench_slots) {
      env <- new.env(parent = globalenv())
      env$.ad <- ad
      env$.slot <- slot

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("get_", short, "_", slot),
          expr = quote(.ad[[.slot]]),
          iterations = iterations,
          env = env
        )
      )
    }

    # --- *_keys() methods ---
    for (method in .bench_key_methods) {
      env <- new.env(parent = globalenv())
      env$.ad <- ad

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("method_", short, "_", method),
          expr = bquote(.ad[[.(method)]]()),
          iterations = iterations,
          env = env
        )
      )
    }

    # --- n_obs, n_vars, shape ---
    for (method in .bench_dim_methods) {
      env <- new.env(parent = globalenv())
      env$.ad <- ad

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("method_", short, "_", method),
          expr = bquote(.ad[[.(method)]]()),
          iterations = iterations,
          env = env
        )
      )
    }

    # --- S3 methods ---
    for (s3name in names(.bench_s3_methods)) {
      env <- new.env(parent = globalenv())
      env$.ad <- ad

      results <- c(
        results,
        run_one_benchmark(
          name = paste0("s3_", short, "_", s3name),
          expr = .bench_s3_methods[[s3name]],
          iterations = iterations,
          env = env
        )
      )
    }

    # Clean up HDF5 handle
    if (backend == "HDF5AnnData") {
      ad$close()
    }
  }

  results
}
