zarr_node_format <- function(store, name) {
  files <- list.files(
    file.path(store, name),
    all.files = TRUE
  )
  if ("zarr.json" %in% files) {
    3L
  } else if (any(c(".zgroup", ".zarray") %in% files)) {
    2L
  } else {
    stop("zarr format cannot be determined!")
  }
}
