#' Build an HDF5 file with an H5AD written at a non-root group
#'
#' Copies the bundled `example.h5ad` fixture into a group at `root` inside a
#' new file, optionally alongside unrelated sibling content, to test reading
#' and writing an `AnnData` from/to a sub-path of an HDF5 file.
#'
#' @param root The group path to copy the example H5AD into
#' @param with_sibling Whether to add unrelated sibling content to the file
#'
#' @return Path to the new HDF5 file
make_nested_h5ad_fixture <- function(root = "mod/rna", with_sibling = TRUE) {
  src <- system.file("extdata", "example.h5ad", package = "anndataR")
  dest <- tempfile(fileext = ".h5ad")

  rhdf5::h5createFile(dest)
  segments <- strsplit(root, "/", fixed = TRUE)[[1]]
  for (i in seq_along(segments)) {
    rhdf5::h5createGroup(dest, paste(segments[seq_len(i)], collapse = "/"))
  }

  if (with_sibling) {
    rhdf5::h5createGroup(dest, "mod/other_stuff")
    rhdf5::h5write("hello", dest, "mod/other_stuff/marker")
  }

  h5src <- rhdf5::H5Fopen(src, flags = "H5F_ACC_RDONLY")
  h5dest <- rhdf5::H5Fopen(dest, flags = "H5F_ACC_RDWR")
  withr::defer(rhdf5::H5Fclose(h5src))
  withr::defer(rhdf5::H5Fclose(h5dest))

  items <- rhdf5::h5ls(h5src, recursive = FALSE)$name
  for (item in items) {
    rhdf5::H5Ocopy(h5src, item, h5dest, paste0(root, "/", item))
  }

  attrs <- rhdf5::h5readAttributes(h5src, "/")
  h5obj <- rhdf5::H5Oopen(h5dest, root)
  withr::defer(rhdf5::H5Oclose(h5obj))
  for (nm in names(attrs)) {
    rhdf5::h5writeAttribute(attrs[[nm]], h5obj, name = nm, asScalar = TRUE)
  }

  dest
}
