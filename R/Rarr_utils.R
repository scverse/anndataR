# Zarr metadata files used to identify valid Zarr nodes (arrays or groups)
ZARR_METADATA_FILES <- c(".zarray", ".zattrs", ".zgroup", "zarr.json")

#' create_zarr_group
#'
#' Create a Zarr group
#'
#' @param store The location of the Zarr store
#' @param name Name of the group
#' @param version Zarr version
#'
#' @return `NULL`
#'
#' @noRd
create_zarr_group <- function(store, name, version = "v2") {
  # Split "a/b/c" into c("a", "b", "c")
  split_name <- strsplit(name, split = "/", fixed = TRUE)[[1]]
  if (length(split_name) > 1) {
    # Build cumulative paths: c("a", "a/b", "a/b/c")
    split_name <- vapply(
      seq_along(split_name),
      function(x) paste(split_name[seq_len(x)], collapse = "/"),
      FUN.VALUE = character(1)
    )
    # Keep only the target and its immediate parent:
    # split_name[1] = "a/b/c" (target), split_name[2] = "a/b" (parent)
    split_name <- rev(tail(split_name, 2))
    # Recursively ensure the parent group exists before creating the target
    if (!dir.exists(file.path(store, split_name[2]))) {
      create_zarr_group(store = store, name = split_name[2])
    }
  }
  dir.create(file.path(store, split_name[1]), showWarnings = FALSE)
  switch(
    version,
    v2 = {
      write(
        "{\"zarr_format\":2}",
        file = file.path(store, split_name[1], ".zgroup")
      )
    },
    v3 = {
      cli_abort("Currently only zarr v2 is supported!")
    },
    cli_abort("Only zarr v2 is supported. Use version = 'v2'")
  )
}

#' create_zarr
#'
#' Create Zarr store
#'
#' @param store The location of the Zarr store
#' @param version Zarr version
#'
#' @return `NULL`
#'
#' @noRd
create_zarr <- function(store, version = "v2") {
  prefix <- basename(store)
  dir <- gsub(paste0(prefix, "$"), "", store)
  create_zarr_group(store = dir, name = prefix, version = version)
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
