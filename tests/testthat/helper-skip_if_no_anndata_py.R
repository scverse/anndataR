# helper function to skip tests if we don't have the Python 'anndata' module
skip_if_no_anndata_py <- function() {
  testthat::skip_if_not_installed("reticulate")
  requireNamespace("reticulate")
  # TEMP(anndata<0.13.0): pin added in PR #481, revert once anndataR supports
  # Python anndata >= 0.13.0
  reticulate::py_require("anndata<0.13.0")
  testthat::skip_if_not(
    reticulate::py_module_available("anndata"),
    message = "Python anndata module not available for testing"
  )
}
