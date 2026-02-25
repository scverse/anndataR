# =============================================================================
# Benchmark suite: subsetting & view materialization
# =============================================================================
# Benchmarks AnnDataView creation with different index types and
# materialization back to concrete implementations.
# =============================================================================

bench_subset <- function(h5ad_paths, iterations) {
  results <- list()
  path <- h5ad_paths[["float_csparse"]]

  for (backend in c("InMemoryAnnData", "HDF5AnnData")) {
    short <- if (backend == "InMemoryAnnData") "InMemory" else "HDF5"
    ad <- read_h5ad(path, as = backend)

    n_obs <- ad$n_obs()
    n_vars <- ad$n_vars()

    # Pre-compute indices
    obs_int <- sort(sample.int(n_obs, as.integer(n_obs * 0.5)))
    var_int <- sort(sample.int(n_vars, as.integer(n_vars * 0.5)))
    obs_char <- ad$obs_names[obs_int]
    var_char <- ad$var_names[var_int]
    obs_logical <- seq_len(n_obs) %in% obs_int

    # --- View creation: integer indices ---
    env <- new.env(parent = globalenv())
    env$.ad <- ad
    env$.obs_int <- obs_int
    env$.var_int <- var_int

    # Subset obs only
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("subset_obs_int_", short),
        expr = quote(.ad[.obs_int, ]),
        iterations = iterations,
        env = env
      )
    )

    # Subset both
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("subset_both_int_", short),
        expr = quote(.ad[.obs_int, .var_int]),
        iterations = iterations,
        env = env
      )
    )

    # --- View creation: character indices ---
    env2 <- new.env(parent = globalenv())
    env2$.ad <- ad
    env2$.obs_char <- obs_char
    env2$.var_char <- var_char

    results <- c(
      results,
      run_one_benchmark(
        name = paste0("subset_obs_char_", short),
        expr = quote(.ad[.obs_char, ]),
        iterations = iterations,
        env = env2
      )
    )

    results <- c(
      results,
      run_one_benchmark(
        name = paste0("subset_both_char_", short),
        expr = quote(.ad[.obs_char, .var_char]),
        iterations = iterations,
        env = env2
      )
    )

    # --- View creation: logical indices ---
    env3 <- new.env(parent = globalenv())
    env3$.ad <- ad
    env3$.obs_logical <- obs_logical

    results <- c(
      results,
      run_one_benchmark(
        name = paste0("subset_obs_logical_", short),
        expr = quote(.ad[.obs_logical, ]),
        iterations = iterations,
        env = env3
      )
    )

    # --- Materialize view → InMemory ---
    view <- ad[obs_int, var_int]
    env4 <- new.env(parent = globalenv())
    env4$.view <- view

    results <- c(
      results,
      run_one_benchmark(
        name = paste0("materialize_to_InMemory_", short),
        expr = quote(.view$as_InMemoryAnnData()),
        iterations = iterations,
        env = env4
      )
    )

    # --- Materialize view → HDF5 ---
    results <- c(
      results,
      run_one_benchmark(
        name = paste0("materialize_to_HDF5_", short),
        expr = quote({
          .tmp <- tempfile(fileext = ".h5ad")
          .result <- .view$as_HDF5AnnData(.tmp)
          .result$close()
          unlink(.tmp)
        }),
        iterations = iterations,
        env = env4
      )
    )

    # Clean up
    if (backend == "HDF5AnnData") {
      ad$close()
    }
  }

  results
}
