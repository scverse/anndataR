# helper function to skip tests if we don't have the Python 'anndata' module
skip_if_no_anndata_py <- function() {
  testthat::skip_if_not_installed("reticulate")
  requireNamespace("reticulate")
  reticulate::py_require("anndata")
  testthat::skip_if_not(
    reticulate::py_module_available("anndata"),
    message = "Python anndata module not available for testing"
  )
}
