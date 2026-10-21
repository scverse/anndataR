# helper function to skip tests if we don't have the Python 'zarr' module
# or the R {anndata} package
skip_if_no_zarr <- function() {
  skip_if_no_py_module("zarr")

  # TODO: Remove when this warning is removed from anndata
  wn <- reticulate::import("warnings")
  wn$filterwarnings(
    "ignore",
    message = "Writing zarr v2 data will no longer be the default"
  )

  # Warnings emitted by Python zarr when it writes the v3 stores used by the
  # roundtrip tests. Filter these by message rather than ignoring everything,
  # so that warnings we do care about still show up.
  wn$filterwarnings(
    "ignore",
    message = "zarr v3 autosharding will be the default"
  )
  wn$filterwarnings(
    "ignore",
    message = "Consolidated metadata is currently not part in the Zarr format 3"
  )
}
