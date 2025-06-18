#' Write a dataset to a HDF5 file
#'
#' Write a HDF5 dataset with chosen compression (can be none)
#'
#' @param file A H5AD file handle
#' @param name Name of the element within the H5AD file
#' @param value Value to write
#' @param compression The compression to use when writing the element. Can be
#'   one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param ... Other arguments passed to [rhdf5::h5write()]
#'
#' @noRd
rhdf5_hdf5_write_dataset <- function(
  file,
  name,
  value,
  compression = c("none", "gzip", "lzf"),
  ...
) {
  compression <- match.arg(compression)

  if (!is.null(dim(value))) {
    dims <- dim(value)
  } else {
    dims <- length(value)
  }

  rhdf5::h5createDataset(
    file,
    name,
    dims,
    storage.mode = storage.mode(value),
    filter = toupper(compression),
    native = TRUE
  )

  rhdf5::h5write(value, file, name, native = TRUE, ...)
}


#' HDF5 path exists
#'
#' Check that a path in HDF5 exists
#'
#' @param file Path to a HDF5 file
#' @param target_path The path within the file to test for
#'
#' @return Whether `target_path` exists in `file`
#' @noRd
rhdf5_hdf5_path_exists <- function(file, target_path) {
  tryCatch(
    {
      rhdf5::H5Lexists(file, target_path)
    },
    error = function(e) {
      FALSE
    }
  )
}

#' Write a HDF5 attribute
#'
#' Write a HDF5 attribute to a HDF5 file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param attr_name Name of the attribute to write
#' @param attr_value Value of the attribute to write
#' @param is_scalar Whether to write the attribute as a scalar value
#'
#' @noRd
rhdf5_hdf5_write_attribute <- function(
  file,
  name,
  attr_name,
  attr_value,
  is_scalar = TRUE
) {
  if (!inherits(file, "H5IdComponent")) {
    cli_abort("{.arg file} must be an open H5AD handle")
  }

  if (name != "/") {
    h5obj <- file&name
    on.exit(rhdf5::H5Oclose(h5obj), add = TRUE)
  } else {
    h5obj <- file
  }

  rhdf5::h5writeAttribute(
    attr = attr_value,
    h5obj = h5obj,
    name = attr_name,
    asScalar = is_scalar
  )
}
