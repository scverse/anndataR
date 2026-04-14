#' create_zarr_group
#'
#' Create a zarr group
#'
#' @param store the location of (zarr) store
#' @param name name of the group
#' @param version zarr version
#'
#' @return `NULL`
#'
#' @noRd
#'
#' @examples
#' store <- tempfile(fileext = ".zarr")
#' create_zarr(store)
#' create_zarr_group(store, "gp")
create_zarr_group <- function(store, name, version = "v2") {
  split_name <- strsplit(name, split = "/", fixed = TRUE)[[1]]
  if (length(split_name) > 1) {
    split_name <- vapply(
      seq_along(split_name),
      function(x) paste(split_name[seq_len(x)], collapse = "/"),
      FUN.VALUE = character(1)
    )
    split_name <- rev(tail(split_name, 2))
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
#' Create zarr store
#'
#' @param store the location of zarr store
#' @param version zarr version
#'
#' @return `NULL`
#'
#' @noRd
#'
#' @examples
#' store <- tempfile(fileext = ".zarr")
#' create_zarr(store)
create_zarr <- function(store, version = "v2") {
  prefix <- basename(store)
  dir <- gsub(paste0(prefix, "$"), "", store)
  create_zarr_group(store = dir, name = prefix, version = version)
}

#' is_zarr_empty
#'
#' check if a zarr store is empty or not.
#'
#' @param store the location of zarr store
#'
#' @return returns TRUE if zarr store is empty
#'
#' @noRd
#'
#' @examples
#' store <- tempfile(fileext = ".zarr")
#' create_zarr(store)
#' is_zarr_empty(store)
is_zarr_empty <- function(store) {
  files <- list.files(store, recursive = FALSE, full.names = FALSE)
  all(files %in% c(".zarray", ".zattrs", ".zgroup", "zarr.json"))
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
    if (any(c(".zarray", ".zattrs", ".zgroup", "zarr.json") %in% list_files)) {
      TRUE
    } else {
      FALSE
    }
  }
}
