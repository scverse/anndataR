# helper function to skip tests if a Python module is not available
#
# `requirement` is passed to `reticulate::py_require()`, set it to `NULL` to
# skip declaring a requirement
#
# Set the `ANNDATAR_TEST_REQUIRE_PYTHON` environment variable to `true` to
# error instead of skipping, e.g. when testing against a specific version of
# Python anndata in CI
skip_if_no_py_module <- function(module, requirement = module) {
  if (!py_tests_required()) {
    testthat::skip_if_not_installed("reticulate")
  }
  if (!is.null(requirement)) {
    reticulate::py_require(requirement)
  }
  if (py_tests_required()) {
    # import directly, py_module_available() would hide the error
    reticulate::import(module)
    return(invisible(TRUE))
  }
  testthat::skip_if_not(
    reticulate::py_module_available(module),
    message = paste("Python", module, "module not available for testing")
  )
}

py_tests_required <- function() {
  isTRUE(as.logical(Sys.getenv("ANNDATAR_TEST_REQUIRE_PYTHON", "false")))
}
