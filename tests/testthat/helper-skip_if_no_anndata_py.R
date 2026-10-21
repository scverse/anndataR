# helper function to skip tests if we don't have the Python 'anndata' module
skip_if_no_anndata_py <- function() {
  # TEMP(anndata<0.13.0): pin added in PR #481, revert once anndataR supports
  # Python anndata >= 0.13.0
  skip_if_no_py_module("anndata", requirement = "anndata<0.13.0")
}
