# helper function to skip tests if a Python module is not available
#
# `requirement` is passed to `reticulate::py_require()`, set it to `NULL` to
# skip declaring a requirement
skip_if_no_py_module <- function(module, requirement = module) {
  testthat::skip_if_not_installed("reticulate")
  if (!is.null(requirement)) {
    reticulate::py_require(requirement)
  }
  testthat::skip_if_not(
    reticulate::py_module_available(module),
    message = paste("Python", module, "module not available for testing")
  )
}
