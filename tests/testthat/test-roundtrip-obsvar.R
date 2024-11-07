skip_if_no_anndata()
skip_if_not_installed("reticulate")
requireNamespace("reticulate")
testthat::skip_if_not(
  reticulate::py_module_available("dummy_anndata"),
  message = "Python dummy_anndata module not available for testing"
)

ad <- reticulate::import("anndata")
da <- reticulate::import("dummy_anndata")

test_names <- names(da$vector_generators)

data <- da$generate_dataset(
  layer_types = list(),
  obsm_types = list(),
  varm_types = list(),
  obsp_types = list(),
  varp_types = list(),
  uns_types = list()
)

for (name in test_names) {
  test_that(paste0("roundtrip with obs and var '", name, "'"), {
    # subset data to make sure only one column is present
    adata_py <-
      ad$AnnData(
        X = data$X, # todo: remove X
        # shape = data$shape, # nolint
        obs = data$obs[, name, drop = FALSE],
        var = data$var[, name, drop = FALSE]
      )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    adata_py$write_h5ad(filename)

    # read from file
    adata_r <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      adata_r$obs[[name]],
      adata_py$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
    expect_equal(
      adata_r$var[[name]],
      adata_py$var[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )

    # write back to file
    filename2 <- withr::local_file(tempfile(fileext = ".h5ad"))
    write_h5ad(adata_r, filename2)

    # read from file
    adata_py2 <- ad$read_h5ad(filename2)

    # expect slots are unchanged
    expect_equal(
      adata_py2$obs[[name]],
      adata_py$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )

    expect_equal(
      adata_py2$var[[name]],
      adata_py$var[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
