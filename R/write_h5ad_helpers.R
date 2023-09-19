#' Write H5AD element
#'
#' Write an element to an H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param ... Additional arguments passed to writing functions
#'
#' @noRd
#'
#' @details
#' `write_h5ad_element()` should always be used instead of any of the specific
#' writing functions as it contains additional boilerplate to make sure
#' elements are written correctly.
write_h5ad_element <- function(value, file, name, compression = "NONE",...) { # nolint

  # Delete the path if it already exists
  if (hdf5_path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

  # Sparse matrices
  if (inherits(value, "sparseMatrix")) {
    write_fun <- write_h5ad_sparse_array
    # Categoricals
  } else if (is.factor(value)) {
    write_fun <- write_h5ad_categorical
    # Lists and data frames
  } else if (is.list(value)) {
    if (is.data.frame(value)) {
      write_fun <- write_h5ad_data_frame
    } else {
      write_fun <- write_h5ad_mapping
    }
    # Character values
  } else if (is.character(value)) {
    if (length(value) == 1) {
      write_fun <- write_h5ad_string_scalar
    } else {
      write_fun <- write_h5ad_string_array
    }
    # Numeric values
  } else if (is.numeric(value)) {
    if (length(value) == 1) {
      write_fun <- write_h5ad_numeric_scalar
    } else if (is.integer(value) && any(is.na(value))) {
      write_fun <- write_h5ad_nullable_integer
    } else {
      write_fun <- write_h5ad_dense_array
    }
    # Logical values
  } else if (is.logical(value)) {
    if (any(is.na(value))) {
      write_fun <- write_h5ad_nullable_boolean
    } else {
      write_fun <- write_h5ad_dense_array
    }
    # Fail if unknown
  } else {
    stop("Writing '", class(value), "' objects to H5AD files is not supported")
  }

  write_fun(value = value, file = file, name = name, compression = compression, ...)
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
  h5file <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(h5file))

  oid <- rhdf5::H5Oopen(h5file, name)
  type <- rhdf5::H5Iget_type(oid)
  rhdf5::H5Oclose(oid)

  if (type == "H5I_GROUP") {
    h5obj <- rhdf5::H5Gopen(h5file, name)
    on.exit(rhdf5::H5Gclose(h5obj), add = TRUE)
  } else {
    h5obj <- rhdf5::H5Dopen(h5file, name)
    on.exit(rhdf5::H5Dclose(h5obj), add = TRUE)
  }

  rhdf5::h5writeAttribute(encoding, h5obj, "encoding-type", asScalar = TRUE) # nolint
  rhdf5::h5writeAttribute(version, h5obj, "encoding-version", asScalar = TRUE) # nolint
}

#' Write H5AD dense array
#'
#' Write a dense array to an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_dense_array <- function(value, file, name, compression, version = "0.2.0") {
  version <- match.arg(version)

  if (!is.vector(value)) {
    # Transpose the value because writing with native=TRUE does not
    # seem to work as expected
    value <- t(value)
  } 

  hdf5_write_compressed(file, name, value, compression)

  # Write attributes
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
  h5file <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(h5file))

  h5obj <- rhdf5::H5Gopen(h5file, name)
  on.exit(rhdf5::H5Gclose(h5obj), add = TRUE)

  rhdf5::h5writeAttribute(dim(value), h5obj, "shape")
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
#' @param version Encoding version of the element to write
write_h5ad_nullable_boolean <- function(value, file, name, compression, version = "0.1.0") {
  # write mask and values
  rhdf5::h5createGroup(file, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE
  
  hdf5_write_compressed(file, paste0(name, "/values"), value_no_na, compression)
  hdf5_write_compressed(file, paste0(name, "/mask"), is.na(value), compression)

  # Write attributes
  write_h5ad_encoding(file, name, "nullable-boolean", version)
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
#' @param version Encoding version of the element to write
write_h5ad_nullable_integer <- function(value, file, name, compression, version = "0.1.0") {
  # write mask and values
  rhdf5::h5createGroup(file, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- -1L
  
  hdf5_write_compressed(file, paste0(name, "/values"), value_no_na, compression)
  hdf5_write_compressed(file, paste0(name, "/mask"), is.na(value), compression)

  # Write attributes
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
#' @param version Encoding version of the element to write
write_h5ad_string_array <- function(value, file, name, compression, version = "0.2.0") {
  rhdf5::h5write(
    value,
    file,
    name,
    variableLengthString = TRUE,
    encoding = "UTF-8")

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
#' @param version Encoding version of the element to write
write_h5ad_categorical <- function(value, file, name, compression, version = "0.2.0") {
  rhdf5::h5createGroup(file, name)
  hdf5_write_compressed(file, paste0(name, "/categories"), levels(value), compression)
  hdf5_write_compressed(file, paste0(name, "/codes"), as.integer(value), compression)
  hdf5_write_compressed(file, paste0(name, "/ordered"), is.ordered(value), compression)

  write_h5ad_encoding(file, name, "categorical", version)
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

  # Write attributes
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
#' @param version Encoding version of the element to write
write_h5ad_numeric_scalar <- function(value, file, name, compression, version = "0.2.0") {
  # Write scalar
  
  # TODO: figure out whether or not a scalar can be compressed? Are there errors?
  # rhdf5::h5write(value, file, name)
  hdf5_write_compressed(file, name, value, compression)
  

  # Write attributes
  write_h5ad_encoding(file, name, "numeric", version)
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
  h5file <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(h5file))

  h5obj <- rhdf5::H5Gopen(h5file, name)
  on.exit(rhdf5::H5Gclose(h5obj), add = TRUE)

  col_order <- colnames(value)
  col_order <- col_order[col_order != index_name]
  # If there are no columns other than the index we set column order to an
  # empty numeric vector
  if (length(col_order) == 0) {
    col_order <- numeric()
  }
  rhdf5::h5writeAttribute(col_order, h5obj, "column-order") # nolint
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
  h5file <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(h5file))

  h5obj <- rhdf5::H5Gopen(h5file, name)
  on.exit(rhdf5::H5Gclose(h5obj), add = TRUE)

  rhdf5::h5writeAttribute(index_name, h5obj, "_index", asScalar = TRUE)
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

hdf5_write_compressed <- function(file, name, value, compression){
  if (!is.null(dim(value))) {
    dims <- dim(value)
  } else {
    dims <- length(value)
  }
  rhdf5::h5createDataset(file, name, dims, filter=compression)
  rhdf5::h5write(value, file, name)
}

