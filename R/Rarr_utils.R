# Zarr metadata files used to identify valid Zarr nodes (arrays or groups)
ZARR_METADATA_FILES <- c(".zarray", ".zattrs", ".zgroup", "zarr.json")

#' is_zarr_empty
#'
#' Check if a Zarr store is empty
#'
#' @param store The location of the Zarr store
#'
#' @return Returns `TRUE` if the Zarr store is empty
#'
#' @noRd
is_zarr_empty <- function(store) {
  files <- list.files(store, recursive = FALSE, full.names = FALSE)
  all(files %in% ZARR_METADATA_FILES)
}

#' Zarr path exists
#'
#' Check that a path in Zarr exists
#'
#' @return Whether the `target_path` exists in `store`
#' @noRd
#'
#' @param store Path to a Zarr store
#' @param target_path The path within the store to test for
zarr_path_exists <- function(store, target_path) {
  zarr <- file.path(store, target_path)
  if (!dir.exists(zarr)) {
    FALSE
  } else {
    list_files <- list.files(
      path = zarr,
      full.names = FALSE,
      recursive = FALSE,
      all.files = TRUE
    )
    if (any(ZARR_METADATA_FILES %in% list_files)) {
      TRUE
    } else {
      FALSE
    }
  }
}
