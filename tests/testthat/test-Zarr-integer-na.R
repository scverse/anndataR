skip_if_not_installed("Rarr")
skip_if_not_installed("rhdf5")

# Integer matrices with NA are written as for H5AD: as double, with NaN for NA
# (#542)
slots <- c("X", "layers", "obsm", "varm", "obsp", "varp")

for (slot in slots) {
  test_that(
    paste0(
      "an integer matrix with NA in ",
      slot,
      " reads back from Zarr as from H5AD"
    ),
    {
      types <- list(
        x_type = if (slot == "X") {
          "integer_matrix_with_nas"
        } else {
          "numeric_matrix"
        },
        layer_types = character(0),
        obs_types = character(0),
        var_types = character(0),
        obsm_types = character(0),
        varm_types = character(0),
        obsp_types = character(0),
        varp_types = character(0),
        uns_types = character(0)
      )
      if (slot != "X") {
        types[[paste0(
          sub("s$", "", slot),
          "_types"
        )]] <- "integer_matrix_with_nas"
      }
      adata <- do.call(
        generate_dataset,
        c(list(n_obs = 10L, n_vars = 20L, format = "AnnData"), types)
      )

      zarr <- withr::local_tempfile(fileext = ".zarr")
      h5ad <- withr::local_tempfile(fileext = ".h5ad")
      write_zarr(adata, zarr)
      write_h5ad(adata, h5ad)

      from_zarr <- read_zarr(zarr)
      from_h5ad <- read_h5ad(h5ad)
      expect_equal(from_zarr[[slot]], from_h5ad[[slot]])
    }
  )
}

test_that("Python anndata reads an integer matrix with NA written to Zarr", {
  skip_if_no_anndata_py()
  adata <- generate_dataset(
    n_obs = 10L,
    n_vars = 20L,
    x_type = "integer_matrix_with_nas",
    layer_types = character(0),
    obs_types = character(0),
    var_types = character(0),
    obsm_types = character(0),
    varm_types = character(0),
    obsp_types = character(0),
    varp_types = character(0),
    uns_types = character(0),
    format = "AnnData"
  )
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_zarr(adata, zarr)

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  x <- reticulate::py_to_r(np$asarray(ad$read_zarr(zarr)$X))
  expect_equal(is.na(x), is.na(adata$X), ignore_attr = TRUE)
  expect_equal(x[!is.na(x)], as.double(adata$X[!is.na(adata$X)]))
})
