#' Determine the Zarr format of a single node in a store
#'
#' Fails if the node carries both v2 and v3 metadata, as such a node is not a
#' valid Zarr node and would otherwise silently be reported as v3.
zarr_node_format <- function(store, name) {
  files <- list.files(
    file.path(store, name),
    all.files = TRUE
  )

  is_v2 <- any(c(".zgroup", ".zarray") %in% files)
  is_v3 <- "zarr.json" %in% files

  if (is_v2 && is_v3) {
    stop("node '", name, "' contains both v2 and v3 metadata!")
  }
  if (!is_v2 && !is_v3) {
    stop("zarr format cannot be determined!")
  }

  if (is_v2) 2L else 3L
}

#' Determine the Zarr format of every node in a store
#'
#' Walks the whole store rather than just the top level, so that a v3 array
#' nested inside a v2 group is detected.
zarr_store_formats <- function(store) {
  dirs <- list.dirs(store, recursive = TRUE, full.names = TRUE)
  nodes <- dirs[
    vapply(
      dirs,
      function(dir) {
        files <- list.files(dir, all.files = TRUE)
        any(ZARR_METADATA_FILES %in% files)
      },
      logical(1)
    )
  ]

  vapply(
    nodes,
    function(node) zarr_node_format(node, ""),
    integer(1),
    USE.NAMES = FALSE
  )
}

#' Expect every node in a store to be written in `format`
expect_zarr_store_format <- function(store, format) {
  formats <- zarr_store_formats(store)
  testthat::expect_gt(length(formats), 0L)
  testthat::expect_equal(unique(formats), as.integer(format))
}
