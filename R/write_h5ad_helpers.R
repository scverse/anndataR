#' Write H5AD element
#'
#' Write an element to an H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' #' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param ... Additional arguments passed to writing functions
#'
#' @noRd
#'
#' @details
#' `write_h5ad_element()` should always be used instead of any of the specific
#' writing functions as it contains additional boilerplate to make sure
#' elements are written correctly.
# nolint start cyclocomp_linter
write_h5ad_element <- function(
    value,
    file,
    name,
    compression = c("none", "gzip", "lzf"),
    stop_on_error = FALSE,
    ...) {
  compression <- match.arg(compression)

  # Delete the path if it already exists
  if (hdf5_path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

  # Sparse matrices
  write_fun <-
    if (inherits(value, "sparseMatrix")) { # Sparse matrices
      write_h5ad_sparse_array
    } else if (is.factor(value)) { # Categoricals
      write_h5ad_categorical
    } else if (is.list(value)) { # Lists and data frames
      if (is.data.frame(value)) {
        write_h5ad_data_frame
      } else {
        write_h5ad_mapping
      }
    } else if (is.character(value)) { # Character values
      if (length(value) == 1 && !is.matrix(value)) {
        write_h5ad_string_scalar
      } else {
        write_h5ad_string_array
      }
    } else if (is.numeric(value) || inherits(value, "denseMatrix")) { # Numeric values
      if (length(value) == 1 && !is.matrix(value)) {
        write_h5ad_numeric_scalar
      } else if (is.integer(value) && any(is.na(value))) {
        write_h5ad_nullable_integer
      } else {
        write_h5ad_dense_array
      }
    } else if (is.logical(value)) { # Logical values
      if (any(is.na(value))) {
        write_h5ad_nullable_boolean
      } else if (length(value) == 1) {
        # Single Booleans should be written as numeric scalars
        write_h5ad_numeric_scalar
      } else {
        write_h5ad_dense_array
      }
    } else { # Fail if unknown
      stop("Writing '", class(value), "' objects to H5AD files is not supported")
    }

  tryCatch(
    {
      write_fun(
        value = value,
        file = file,
        name = name,
        compression = compression,
        ...
      )
    },
    error = function(e) {
      message <- paste0(
        "Could not write element '", name, "' of type '", class(value), "':\n",
        conditionMessage(e)
      )
      if (stop_on_error) {
        stop(message)
      } else {
        warning(message)
        return(NULL)
      }
    }
  )
}
# nolint end cyclocomp_linter

#' Write H5AD attributes
#'
#' Write H5AD attributes to an element in an H5AD file
#'
#' @noRd
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param attributes Named list of attributes to write
#' @param is_scalar Whether to write attributes as scalar values. Can be `TRUE`
#' to write all attributes as scalars, `FALSE` to write no attributes as
#' scalars, or a vector of the names of `attributes` that should be written.
write_h5ad_attributes <- function(file, name, attributes, is_scalar = TRUE) {
  if (is.character(file) && length(file) == 1 && !is.na(file)) {
    h5file <- rhdf5::H5Fopen(file)
    on.exit(rhdf5::H5Fclose(h5file))
  } else if (inherits(file, "H5IdComponent")) {
    h5file <- file
  } else {
    stop("file must be a path to an H5AD file or an open H5AD handle")
  }

  oid <- rhdf5::H5Oopen(h5file, name)
  type <- rhdf5::H5Iget_type(oid)
  rhdf5::H5Oclose(oid)

  if (isTRUE(is_scalar)) {
    is_scalar <- names(attributes)
  } else if (isFALSE(is_scalar)) {
    is_scalar <- c()
  } else if (is.character(is_scalar)) {
    is_scalar <- is_scalar
  } else {
    stop("is_scalar must be TRUE, FALSE or a character vector")
  }

  if (type == "H5I_GROUP") {
    h5obj <- rhdf5::H5Gopen(h5file, name)
    on.exit(rhdf5::H5Gclose(h5obj), add = TRUE)
  } else {
    h5obj <- rhdf5::H5Dopen(h5file, name)
    on.exit(rhdf5::H5Dclose(h5obj), add = TRUE)
  }

  for (attr_name in names(attributes)) {
    attr_value <- attributes[[attr_name]]
    if (isTRUE(attr_value) || isFALSE(attr_value)) {
      write_h5ad_boolean_attribute(attr_value, h5obj, attr_name)
    } else {
      scalar_value <- attr_name %in% is_scalar
      rhdf5::h5writeAttribute(
        attr_value, h5obj, attr_name, asScalar = scalar_value
      ) # nolint
    }
  }
}

#' Write a H5AD Boolean attribute
#'
#' Write a Boolean attribute to a HDF5 element
#'
#' @noRd
#'
#' @param attr_value Boolean value to write
#' @param h5obj Object representing the HDF5 element to write the attribute to
#' @param attr_name Name of the attribute to write
write_h5ad_boolean_attribute <- function(attr_value, h5obj, attr_name) {
  # Based on https://github.com/grimbough/rhdf5/issues/136#issuecomment-1940860832

  # Create an ENUM datatype and insert our two key:value pairs
  type_id <- rhdf5::H5Tenum_create(dtype_id = "H5T_NATIVE_UCHAR")
  rhdf5::H5Tenum_insert(type_id, name = "TRUE", value = 1L)
  rhdf5::H5Tenum_insert(type_id, name = "FALSE", value = 0L)

  # Create the attribute using our new data type
  space_id <- rhdf5::H5Screate()
  on.exit(rhdf5::H5Sclose(space_id), add = TRUE)
  attribute_id <- rhdf5::H5Acreate(
    h5obj = h5obj, name = attr_name, dtype_id = type_id, h5space = space_id
  )
  on.exit(rhdf5::H5Aclose(attribute_id), add = TRUE)

  rhdf5::H5Awrite(attribute_id, buf = attr_value)
}

#' Write H5AD encoding
#'
#' Write H5AD encoding attributes to an element in an H5AD file
#'
#' @noRd
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param encoding The encoding type to set
#' @param version The encoding version to set
write_h5ad_encoding <- function(file, name, encoding, version) {
  write_h5ad_attributes(
    file = file,
    name = name,
    attributes = list(
      "encoding-type" = encoding,
      "encoding-version" = version
    )
  )
}

#' Write H5AD dense array
#'
#' Write a dense array to an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
#'
#' @noRd
write_h5ad_dense_array <- function(value, file, name, compression, version = "0.2.0") {
  version <- match.arg(version)

  if (!is.vector(value)) {
    # Transpose the value because writing with native=TRUE does not
    # seem to work as expected
    value <- t(value)
  }

  hdf5_write_compressed(file, name, value, compression)

  # Write encoding
  write_h5ad_encoding(file, name, "array", version)
}

#' Write H5AD sparse array
#'
#' Write a sparse array to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_sparse_array <- function(value, file, name, compression, version = "0.1.0") {
  version <- match.arg(version)

  # check types
  stopifnot(inherits(value, "sparseMatrix"))

  if (inherits(value, "RsparseMatrix")) {
    type <- "csr_matrix"
    indices_attr <- "j"
  } else if (inherits(value, "CsparseMatrix")) {
    type <- "csc_matrix"
    indices_attr <- "i"
  } else {
    stop(
      "Unsupported matrix format in ", name, ".",
      "Supported formats are RsparseMatrix and CsparseMatrix",
      "(and objects that inherit from those)."
    )
  }

  # Write sparse matrix
  rhdf5::h5createGroup(file, name)
  hdf5_write_compressed(file, paste0(name, "/indices"), attr(value, indices_attr), compression)
  hdf5_write_compressed(file, paste0(name, "/indptr"), value@p, compression)
  hdf5_write_compressed(file, paste0(name, "/data"), value@x, compression)

  # Add encoding
  write_h5ad_encoding(file, name, type, version)

  # Write shape attribute
  write_h5ad_attributes(
    file, name, list("shape" = dim(value)),
    is_scalar = FALSE
  )
}

#' Write H5AD nullable boolean
#'
#' Write a nullable boolean to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_nullable_boolean <- function(value, file, name, compression, version = "0.1.0") {
  # Write mask and values
  rhdf5::h5createGroup(file, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE

  # hdf5_write_compressed(file, paste0(name, "/values"), value_no_na, compression)
  # hdf5_write_compressed(file, paste0(name, "/mask"), is.na(value), compression)

  write_h5ad_boolean_array(value_no_na, file, paste0(name, "/values"))
  write_h5ad_boolean_array(is.na(value), file, paste0(name, "/mask"))

  # Write encoding
  write_h5ad_encoding(file, name, "nullable-boolean", version)
}

#' Write H5AD boolean array
#'
#' Write a boolean array to an H5AD file as a custom ENUM dataset
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
write_h5ad_boolean_array <- function(value, file, name) {
  # Based on https://stackoverflow.com/a/74653515/4384120

  h5file <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(h5file))

  value <- as.integer(value)
  if (!is.null(dim(value))) {
    dims <- dim(value)
  } else {
    dims <- length(value)
  }

  # Create the dataspace for the data to write
  h5space <- rhdf5::H5Screate_simple(dims = dims, NULL, native = TRUE)
  on.exit(rhdf5::H5Sclose(h5space), add = TRUE)

  # Create the Boolean ENUM datatype (TRUE = 1 FALSE = 0)
  tid <- rhdf5::H5Tenum_create(dtype_id = "H5T_NATIVE_UCHAR")
  rhdf5::H5Tenum_insert(tid, name = "TRUE", value = 1L)
  rhdf5::H5Tenum_insert(tid, name = "FALSE", value = 0L)

  # Create the dataset with this new datatype
  h5dataset <- rhdf5::H5Dcreate(h5file, name, tid, h5space)
  on.exit(rhdf5::H5Dclose(h5dataset), add = TRUE)

  # Write the data. We have to use as.raw() because our base type is 8-bit and
  # R integers are 32-bit
  rhdf5::H5Dwrite(h5dataset, as.raw(value), h5type = tid)
}


#' Write H5AD nullable integer
#'
#' Write a nullable integer to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_nullable_integer <- function(value, file, name, compression, version = "0.1.0") {
  # write mask and values
  rhdf5::h5createGroup(file, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- -1L

  hdf5_write_compressed(file, paste0(name, "/values"), value_no_na, compression)
  hdf5_write_compressed(file, paste0(name, "/mask"), is.na(value), compression)

  # Write encoding
  write_h5ad_encoding(file, name, "nullable-integer", version)
}

#' Write H5AD string array
#'
#' Write a string array to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_string_array <- function(value, file, name, compression, version = "0.2.0") {
  rhdf5::h5write(
    value,
    file,
    name,
    variableLengthString = TRUE,
    encoding = "UTF-8"
  )

  write_h5ad_encoding(file, name, "string-array", version)
}

#' Write H5AD categorical
#'
#' Write a categorical to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_categorical <- function(value, file, name, compression, version = "0.2.0") {
  rhdf5::h5createGroup(file, name)

  categories <- levels(value)

  # Use zero-indexed values
  codes <- as.integer(value) - 1

  # Set missing values to -1
  codes[is.na(codes)] <- -1

  # write values to file
  hdf5_write_compressed(file, paste0(name, "/categories"), categories, compression)
  hdf5_write_compressed(file, paste0(name, "/codes"), codes, compression)

  # Write encoding
  write_h5ad_encoding(file = file, name = name, encoding = "categorical", version = version)

  # Write ordered attribute
  write_h5ad_attributes(file, name, list("ordered" = is.ordered(value)))
}

#' Write H5AD string scalar
#'
#' Write a string scalar to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_string_scalar <- function(value, file, name, compression, version = "0.2.0") {
  # Write scalar
  rhdf5::h5write(
    value,
    file,
    name,
    variableLengthString = TRUE,
    encoding = "UTF-8"
  )

  # Write encoding
  write_h5ad_encoding(file, name, "string", version)
}

#' Write H5AD numeric scalar
#'
#' Write a numeric scalar to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_numeric_scalar <- function(value, file, name, compression, version = "0.2.0") {
  # Write scalar
  hdf5_write_compressed(file, name, value, compression)

  # Write encoding
  write_h5ad_encoding(file, name, "numeric-scalar", version)
}

#' Write H5AD mapping
#'
#' Write a mapping to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_h5ad_mapping <- function(value, file, name, compression, version = "0.1.0") {
  rhdf5::h5createGroup(file, name)

  # Write mapping elements
  for (key in names(value)) {
    write_h5ad_element(value[[key]], file, paste0(name, "/", key), compression)
  }

  write_h5ad_encoding(file, name, "dict", version)
}

#' Write H5AD data frame
#'
#' Write a data frame to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param index The index to write. Can either be a vector of length equal to
#' the number of rows in `values` or a single character string giving the name
#' of a column in `values`. If `NULL` then `rownames(value)` is used.
#' @param version Encoding version of the element to write
write_h5ad_data_frame <- function(value, file, name, compression, index = NULL,
                                  version = "0.2.0") {
  rhdf5::h5createGroup(file, name)
  write_h5ad_encoding(file, name, "dataframe", version)

  if (is.null(index)) {
    index_name <- "_index"
    index_value <- rownames(value)
  } else if (length(index) == nrow(value)) {
    index_name <- "_index"
    index_value <- index
  } else if (length(index) == 1 && index %in% colnames(value)) {
    index_name <- index
    index_value <- value[[index_name]]
    value[[index_name]] <- NULL
  } else {
    stop(
      "index must be a vector with length `nrow(value)` or a single character",
      "string giving the name of a column in `value`"
    )
  }

  # Write index
  write_h5ad_data_frame_index(index_value, file, name, compression, index_name)

  # Write data frame columns
  for (col in colnames(value)) {
    write_h5ad_element(value[[col]], file, paste0(name, "/", col), compression)
  }

  # Write additional data frame attributes
  col_order <- colnames(value)
  col_order <- col_order[col_order != index_name]
  # If there are no columns other than the index we set column order to an
  # empty numeric vector
  if (length(col_order) == 0) {
    col_order <- numeric()
  }

  write_h5ad_attributes(
    file, name, list("column-order" = col_order),
    is_scalar = FALSE
  )
}

#' Write H5AD data frame index
#'
#' Write an for a data frame to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write. Must be a vector to the same length as the data
#' frame.
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file containing the data
#' frame
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param index_name Name of the data frame column storing the index
write_h5ad_data_frame_index <- function(value, file, name, compression, index_name) {
  if (!hdf5_path_exists(file, name)) {
    stop("The data frame '", name, "' does not exist in '", file, "'")
  }

  encoding <- read_h5ad_encoding(file, name)
  if (encoding$type != "dataframe") {
    stop("'", name, "' in '", file, "' is not a data frame")
  }

  # Write index columns
  write_h5ad_element(value, file, paste0(name, "/", index_name))

  # Write data frame index attribute
  write_h5ad_attributes(file, name, list("_index" = index_name))
}

#' Write empty H5AD
#'
#' Write a new empty H5AD file
#'
#' @noRd
#'
#' @param file Path to the H5AD file to write
#' @param obs_names Vector containing observation names
#' @param var_names Vector containing variable names
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version The H5AD version to write
write_empty_h5ad <- function(file, obs_names, var_names, compression, version = "0.1.0") {
  h5file <- rhdf5::H5Fcreate(file)
  rhdf5::H5Fclose(h5file)

  write_h5ad_encoding(file, "/", "anndata", "0.1.0")

  write_h5ad_element(data.frame(row.names = obs_names), file, "/obs", compression)
  write_h5ad_element(data.frame(row.names = var_names), file, "/var", compression)

  rhdf5::h5createGroup(file, "layers")
  write_h5ad_encoding(file, "/layers", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "obsm")
  write_h5ad_encoding(file, "/obsm", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "obsp")
  write_h5ad_encoding(file, "/obsp", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "uns")
  write_h5ad_encoding(file, "/uns", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "varm")
  write_h5ad_encoding(file, "/varm", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "varp")
  write_h5ad_encoding(file, "/varp", "dict", "0.1.0")
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
  if (substr(target_path, 1, 1) != "/") {
    target_path <- paste0("/", target_path)
  }

  content <- rhdf5::h5ls(file)

  paths <- file.path(content$group, content$name)
  paths <- gsub("//", "/", paths) # Remove double slash for root paths

  target_path %in% paths
}

#' HDF5 write compressed
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
#'
#' @return Whether the `path` exists in `file`
hdf5_write_compressed <- function(file, name, value, compression = c("none", "gzip", "lzf")) {

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
    filter = toupper(compression)
  )

  rhdf5::h5write(value, file, name)
}
