# Zarr metadata files used to identify valid Zarr nodes (arrays or groups)
ZARR_METADATA_FILES <- c(".zarray", ".zattrs", ".zgroup", "zarr.json")

#' check_zarr_format
#'
#' Check that a Zarr format is one of the supported versions
#'
#' @param format Zarr format
#'
#' @return `format` as an integer
#'
#' @noRd
check_zarr_format <- function(format, call = rlang::caller_env()) {
  if (
    length(format) != 1L ||
    !is.numeric(format) ||
    is.na(format) ||
    !format %in% c(2L, 3L)
  ) {
    cli_abort(
      "{.arg zarr_format} must be either {.val {2L}} or {.val {3L}}, \\
      not {.val {format}}.",
      call = call
    )
  }
  
  as.integer(format)
}

#' get_zarr_format
#'
#' Determine the Zarr format of an existing store from its root metadata
#'
#' @param store The location of the Zarr store
#'
#' @return The Zarr format of `store` as an integer
#'
#' @noRd
get_zarr_format <- function(store, call = rlang::caller_env()) {
  files <- list.files(store, all.files = TRUE, recursive = FALSE)

  is_v2 <- any(c(".zgroup", ".zarray") %in% files)
  is_v3 <- "zarr.json" %in% files

  if (is_v2 && is_v3) {
    cli_abort(
      c(
        "Store {.file {store}} contains both Zarr v2 and v3 root metadata.",
        "i" = "Found both {.file .zgroup}/{.file .zarray} and {.file zarr.json}."
      ),
      call = call
    )
  }

  if (!is_v2 && !is_v3) {
    cli_abort(
      "Could not determine the Zarr format of {.file {store}}.",
      call = call
    )
  }

  if (is_v2) 2L else 3L
}

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
