# helper function to skip tests if we don't have the Python 'zarr' module
# or the R {anndata} package
skip_if_no_zarr <- function() {
  testthat::skip_if_not_installed("reticulate")
  reticulate::py_require("zarr")
  testthat::skip_if_not(
    reticulate::py_module_available("zarr"),
    message = "Python zarr module not available for testing"
  )

  # TODO: Remove when this warning is removed from anndata
  wn <- reticulate::import("warnings")
  wn$filterwarnings(
    "ignore",
    message = "Writing zarr v2 data will no longer be the default"
  )
}
