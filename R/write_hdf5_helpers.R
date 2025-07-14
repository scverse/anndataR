#' Write a dataset to a HDF5 file
#'
#' Write a HDF5 dataset with chosen compression (can be none)
#'
#' @param file A HDF5 file handle
#' @param name Name of the element within the H5AD file
#' @param value Value to write
#' @param H5type Datatype to write, see [rhdf5::h5createDataset()]
#' @param compression The compression to use when writing the element. Can be
#'   one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param ... Other arguments passed to [rhdf5::h5write()]
#'
#' @noRd
hdf5_write_dataset <- function(
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

#' Write a scalar to a HDF5 file
#'
#' Write a HDF5 scalar
#'
#' @param file A HDF5 file handle
#' @param name Name of the element within the H5AD file
#' @param value Value to write
#'
#' @noRd
hdf5_write_scalar <- function(
  file,
  name,
  value
) {
  # Missing values need to be stored as floats
  if (is.integer(value) && is.na(value)) {
    value <- as.numeric(value)
  }

  h5space <- rhdf5::H5Screate("H5S_SCALAR", native = TRUE)
  on.exit(rhdf5::H5Sclose(h5space), add = TRUE)

  # Adjusted based on rhdf5:::.setDataType
  tid <- switch(
    storage.mode(value)[1],
    integer = "H5T_STD_I64LE",
    double = "H5T_IEEE_F64LE",
    character = {
      tid <- rhdf5::H5Tcopy("H5T_C_S1")
      rhdf5::H5Tset_strpad(tid, strpad = "NULLTERM")
      rhdf5::H5Tset_size(tid, NULL)
      rhdf5::H5Tset_cset(tid, "UTF-8")
      tid
    },
    {
      cli_abort(c(
        "{.arg value} must be of type {.val integer}, {.val double} or {.val character}",
        "i" = "Got {.val {storage.mode(value)[1]}}"
      ))
    }
  )

  dcpl <- rhdf5::H5Pcreate("H5P_DATASET_CREATE")
  on.exit(rhdf5::H5Pclose(dcpl), add = TRUE)
  rhdf5::H5Pset_fill_time(dcpl, "H5D_FILL_TIME_ALLOC")
  rhdf5::H5Pset_obj_track_times(dcpl, FALSE)

  # Create the dataset with this new datatype
  h5dataset <- rhdf5::H5Dcreate(
    h5loc = file,
    name = name,
    dtype_id = tid,
    h5space = h5space,
    dcpl = dcpl
  )
  on.exit(rhdf5::H5Dclose(h5dataset), add = TRUE)

  rhdf5::H5Dwrite(h5dataset, value)
}

#' Write a boolean dataset to a HDF5 array
#'
#' Write a boolean dataset with chosen compression (can be none)
#'
#' @param file A HDF5 file handle
#' @param name Name of the element within the H5AD file
#' @param value Value to write
#' @param is_scalar Whether to use a scalar data space
#' @param compression The compression to use when writing the element. Can be
#'   one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`. CURRENTLY NOT
#'   IMPLEMENTED!
#'
#' @noRd
hdf5_write_boolean_dataset <- function(
  file,
  name,
  value,
  is_scalar = FALSE,
  compression = c("none", "gzip", "lzf")
) {
  # Based on https://stackoverflow.com/a/74653515/4384120

  # TODO: Add compression
  # nolint start
  # dcpl <- rhdf5::H5Pcreate("H5P_DATASET_CREATE")
  # on.exit(rhdf5::H5Pclose(dcpl), add = TRUE)
  # if (compression == "gzip") {
  #   rhdf5::H5Pset_deflate(dcpl, 6L) # 6 is a good default compression level
  # } else if (compression == "lzf") {
  #   rhdf5::H5Pset_lzf(dcpl, h5tid = tid)
  # }
  # nolint end

  value <- as.integer(value)

  # Create the dataspace for the data to write
  if (is_scalar) {
    h5space <- rhdf5::H5Screate("H5S_SCALAR", native = TRUE)
  } else {
    if (!is.null(dim(value))) {
      dims <- dim(value)
    } else {
      dims <- length(value)
    }

    h5space <- rhdf5::H5Screate_simple(dims = dims, NULL, native = TRUE)
  }
  on.exit(rhdf5::H5Sclose(h5space), add = TRUE)

  # Create the Boolean ENUM datatype (TRUE = 0 FALSE = 1)
  tid <- rhdf5::H5Tenum_create(dtype_id = "H5T_NATIVE_SCHAR")
  rhdf5::H5Tenum_insert(tid, name = "FALSE", value = 0L)
  rhdf5::H5Tenum_insert(tid, name = "TRUE", value = 1L)

  dcpl <- rhdf5::H5Pcreate("H5P_DATASET_CREATE")
  on.exit(rhdf5::H5Pclose(dcpl), add = TRUE)
  rhdf5::H5Pset_fill_time(dcpl, "H5D_FILL_TIME_ALLOC")
  rhdf5::H5Pset_obj_track_times(dcpl, FALSE)

  # Create the dataset with this new datatype
  h5dataset <- rhdf5::H5Dcreate(
    h5loc = file,
    name = name,
    dtype_id = tid,
    h5space = h5space,
    dcpl = dcpl
  )
  on.exit(rhdf5::H5Dclose(h5dataset), add = TRUE)

  # Write the data. We have to use as.raw() because our base type is 8-bit and
  # R integers are 32-bit
  rhdf5::H5Dwrite(h5dataset, as.raw(value), h5type = tid)
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
hdf5_path_exists <- function(file, target_path) {
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
hdf5_write_attribute <- function(
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
hdf5_clear_file <- function(h5file) {
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
hdf5_clear_rhdf5_attributes <- function(h5file, name) {
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
      hdf5_clear_rhdf5_attributes(h5file, paste0(name, "/", item))
    }
  }
}
