#' Create a dataset in a HDF5 file
#'
#' Write HDF5 dataset with chosen compression (can be none)
#'
#' @noRd
#'
#' @param file Path to a HDF5 file
#' @param name Name of the element within the H5AD file containing the data
#' frame
#' @param value Value to write. Must be a vector to the same length as the data
#' frame.
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param scalar Whether to write the value as a scalar or not
#' @param dtype The data type of the attribute to write. If `NULL` then the data
#'  type is guessed using `hdf5r::guess_dtype()`.
#' @param space The space of the attribute to write. If `NULL` then the space
#'   is guessed using `hdf5r::guess_space()`.
#'
#' @return Whether the `path` exists in `file`
hdf5_create_dataset <- function(
    file,
    name,
    value,
    compression = c("none", "gzip", "lzf"),
    scalar = FALSE,
    dtype = NULL,
    space = NULL) {
  compression <- match.arg(compression)

  if (!is.null(dim(value))) {
    dims <- rev(dim(value))
  } else {
    dims <- length(value)
  }

  if (is.null(dtype)) {
    dtype <- hdf5r::guess_dtype(value, scalar = scalar, string_len = Inf)
  }

  if (is.null(space)) {
    space <- hdf5r::guess_space(value, dtype = dtype, chunked = FALSE)
  }

  # TODO: lzf compression is not supported in hdf5r
  # TODO: allow the user to choose compression level
  gzip_level <- if (compression == "none") 0 else 9

  out <- file$create_dataset(
    name = name,
    dims = dims,
    gzip_level = gzip_level,
    robj = value,
    chunk_dims = NULL,
    space = space,
    dtype = dtype
  )
  # todo: create attr?

  out
}


#' HDF5 path exists
#'
#' Check that a path in HDF5 exists
#'
#' @noRd
#'
#' @param file Path to a HDF5 file
#' @param target_path The path within the file to test for
#'
#' @return Whether the `path` exists in `file`
hdf5_path_exists <- function(file, target_path) {
  tryCatch(
    {
      file$exists(target_path)
    },
    error = function(e) {
      FALSE
    }
  )
}

#' Create a HDF5 attribute
#'
#' Create a HDF5 attribute in a HDF5 file
#'
#' @noRd
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param attr_name Name of the attribute to write
#' @param attr_value Value of the attribute to write
#' @param is_scalar Whether to write attributes as scalar values. Can be `TRUE`
#' to write all attributes as scalars, `FALSE` to write no attributes as
#' scalars, or a vector of the names of `attributes` that should be written.
#' @param dtype The data type of the attribute to write. If `NULL` then the data
#'   type is guessed using `hdf5r::guess_dtype()`.
#' @param space The space of the attribute to write. If `NULL` and `is_scalar` then
#'   the space is set to a scalar space. If `NULL` and `!is_scalar` then the space
#'   is guessed using `hdf5r::guess_space()`.
hdf5_create_attribute <- function(
    file,
    name,
    attr_name,
    attr_value,
    is_scalar = TRUE,
    dtype = NULL,
    space = NULL) {
  if (!inherits(file, "H5File")) {
    stop("file must be an open H5AD handle")
  }

  if (is.null(dtype)) {
    dtype <- hdf5r::guess_dtype(attr_value, scalar = is_scalar, string_len = Inf)
  }
  if (is.null(space)) {
    space <-
      if (is_scalar) {
        hdf5r::H5S$new(type = "scalar")
      } else {
        hdf5r::guess_space(attr_value, dtype = dtype, chunked = FALSE)
      }
  }
  file$create_attr_by_name(
    attr_name = attr_name,
    obj_name = name,
    robj = attr_value,
    dtype = dtype,
    space = space
  )
}
