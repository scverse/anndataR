# helper function to skip tests if we don't have the Python 'anndata' module
skip_if_no_anndata_py <- function() {
  skip_if_no_py_module("anndata") # nolint: object_usage_linter
}
