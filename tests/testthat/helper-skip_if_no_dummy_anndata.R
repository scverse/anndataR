# helper function to skip tests if we don't have the Python 'anndata' module
# or the R {anndata} package
skip_if_no_dummy_anndata <- function() {
  testthat::skip_if_not_installed("reticulate")
  requireNamespace("reticulate")
  reticulate::py_require("dummy_anndata")
  testthat::skip_if_not(
    reticulate::py_module_available("dummy_anndata"),
    message = "Python dummy_anndata module not available for testing"
  )
}
