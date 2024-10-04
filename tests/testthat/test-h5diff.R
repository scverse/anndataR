skip_if_not_installed("hdf5r")

requireNamespace("reticulate")
testthat::skip_if_not(
  reticulate::py_module_available("dummy_anndata"),
  message = "Python dummy_anndata module not available for testing"
)

test_that("h5diff", {
  requireNamespace("processx")
  
  data_R <- generate_dataset(10L, 20L, format = "AnnData")

  da <- reticulate::import("dummy_anndata")
  data_python <- da$generate_dataset(10L, 20L)

  h5ad_file1 <- tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
  h5ad_file2 <- tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")

  write_h5ad(data_R, h5ad_file1)
  data_python$write_h5ad(h5ad_file2)
  res <- processx::run("h5diff", c("-v", h5ad_file1, h5ad_file2), error_on_status = FALSE)

  expect_equal(res$status, 0, info = res$stdout)

})

