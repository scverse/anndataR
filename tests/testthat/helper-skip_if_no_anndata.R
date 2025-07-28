# helper function to skip tests if we don't have the Python 'anndata' module
# or the R {anndata} package
skip_if_no_anndata <- function() {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("anndata")
  requireNamespace("reticulate")
  reticulate::py_require("anndata<=0.11.4")
  testthat::skip_if_not(
    reticulate::py_module_available("anndata"),
    message = "Python anndata module not available for testing"
  )
}
