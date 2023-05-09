# helper function to skip tests if we don't have the 'anndata' module
skip_if_no_anndata <- function() {
  requireNamespace("reticulate")
  testthat::skip_if_not(
    reticulate::py_module_available("anndata"),
    message = "anndata not available for testing"
  )
}