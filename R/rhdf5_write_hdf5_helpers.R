#' Write a dataset to a HDF5 file
#'
#' Write a HDF5 dataset with chosen compression (can be none)
#'
#' @param file A H5AD file handle
#' @param name Name of the element within the H5AD file
#' @param value Value to write
#' @param H5type Datatype to write, see [rhdf5::h5createDataset()]
#' @param compression The compression to use when writing the element. Can be
#'   one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param ... Other arguments passed to [rhdf5::h5write()]
#'
#' @noRd
rhdf5_hdf5_write_dataset <- function(
  file,
  name,
  value,
  H5type = NULL,
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
    H5type = H5type,
    filter = toupper(compression),
    native = FALSE
  )

  rhdf5::h5write(value, file, name, native = FALSE, ...)
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
    h5obj <- rhdf5::H5Oopen(file, name)
    on.exit(rhdf5::H5Oclose(h5obj), add = TRUE)
  } else {
    h5obj <- file
  }

  rhdf5::h5writeAttribute(
    attr = attr_value,
    h5obj = h5obj,
    name = attr_name,
    asScalar = is_scalar,
    encoding = "UTF-8",
    variableLengthString = TRUE
  )
}

#' Clear HDF5 file
#'
#' Remove all contents from a HDF5 file
#'
#' @param h5file A HDF5 file handle
#'
#' @noRd
rhdf5_hdf5_clear_file <- function(h5file) {
  if (!inherits(h5file, "H5IdComponent")) {
    cli_abort("{.arg h5file} must be an open H5AD handle")
  }

  contents <- rhdf5::h5ls(h5file, recursive = FALSE)

  for (item in contents$name) {
    rhdf5::h5delete(h5file, item)
  }
}

#' Clear rhdf5 attributes
#'
#' Remove attributes added by rhdf5 from an element in a HDF5 file
#'
#' @param h5file A HDF5 file handle or a path to a HDF5 file
#' @name name Name of the element within the HDF5 file
#'
#' @details
#' rhdf5 adds a `rhdf5-NA.OK` attribute to elements which doesn't affect the
#' validity of the HDF5 but shows up when comparing to a file created by Python.
#' This function removes those attributes from an element so files are
#' comparable.
#'
#' @noRd
rhdf5_hdf5_clear_rhdf5_attributes <- function(h5file, name) {

  if (!inherits(h5file, "H5IdComponent")) {
    h5file <- rhdf5::H5Fopen(h5file)
    on.exit(rhdf5::H5Fclose(h5file), add = TRUE)
  }

  h5obj <- rhdf5::H5Oopen(h5file, name)
  h5type <- rhdf5::H5Iget_type(h5obj)
  rhdf5_na_ok_exists <- rhdf5::H5Aexists(h5obj, "rhdf5-NA.OK")
  rhdf5::H5Oclose(h5obj)

  if (rhdf5_na_ok_exists) {
    rhdf5::h5deleteAttribute(h5file, name, "rhdf5-NA.OK")
  }

  if (h5type == "H5I_GROUP") {
    h5group <- rhdf5::H5Gopen(h5file, name)
    contents <- rhdf5::h5ls(h5group, recursive = FALSE)
    rhdf5::H5Gclose(h5group)

    for (item in contents$name) {
      rhdf5_hdf5_clear_rhdf5_attributes(h5file, paste0(name, "/", item))
    }
  }
}
