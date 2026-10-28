# helper function to skip tests if we don't have the Python 'anndata' module
# or the R {anndata} package
skip_if_no_dummy_anndata <- function() {
  skip_if_no_py_module("dummy_anndata") # nolint: object_usage_linter
}
