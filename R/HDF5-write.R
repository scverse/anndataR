#' Write H5AD element
#'
#' Write an element to a H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param ... Additional arguments passed to writing functions
#' 
#' @details
#' `write_h5ad_element()` should always be used instead of any of the specific
#' writing functions as it contains additional boilerplate to make sure
#' elements are written correctly.
write_h5ad_element <- function(value, file, name, ...) { # nolint
  
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
    stop("Write '", class(value), "' objects to H5AD files is not supported")
  }
  
  write_fun(value = value, file = file, name = name, ...)
}

#' Write H5AD encoding
#'
#' Write H5AD encoding attributes to an element in an H5AD file
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
  
  rhdf5::h5writeAttribute(encoding, h5obj, "encoding-type") # nolint
  rhdf5::h5writeAttribute(version, h5obj, "encoding-version") # nolint
}

#' Write H5AD dense array
#'
#' Write a dense array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_dense_array <- function(value, file, name, version = "0.2.0") {
  requireNamespace("rhdf5")

  version <- match.arg(version)

  # Transpose the value because writing with native=TRUE does not
  # seem to work as expected
  rhdf5::h5write(t(value), file, name)

  # Write attributes
  write_h5ad_encoding(file, name, "array", version)
}

#' Write H5AD sparse array
#'
#' Write a sparse array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_sparse_array <- function(value, file, name, version = "0.1.0") {
  requireNamespace("rhdf5")

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
      "(and objects that inherit from those).")
  }

  # Write sparse matrix
  rhdf5::h5createGroup(file, name)
  rhdf5::h5write(attr(value, indices_attr), file, paste0(name, "/indices"))
  rhdf5::h5write(value@p, file, paste0(name, "/indptr"))
  rhdf5::h5write(value@x, file, paste0(name, "/data"))

  # Add encoding
  write_h5ad_encoding(file, name, type, version)
}

#' Write H5AD nullable boolean
#'
#' Write a nullable boolean from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_nullable_boolean <- function(value, file, name, version = "0.1.0") {
  requireNamespace("rhdf5")
  
  # write mask and values
  rhdf5::h5createGroup(file, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE
  rhdf5::h5write(value_no_na, file, paste0(name, "/values"))
  rhdf5::h5write(is.na(value), file, paste0(name, "/mask"))

  # Write attributes
  write_h5ad_encoding(file, name, "nullable-boolean", version)
}

#' Write H5AD nullable integer
#'
#' Write a nullable integer from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_nullable_integer <- function(value, file, name, version = "0.1.0") {
  requireNamespace("rhdf5")

  # write mask and values
  rhdf5::h5createGroup(file, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- -1L
  rhdf5::h5write(value_no_na, file, paste0(name, "/values"))
  rhdf5::h5write(is.na(value), file, paste0(name, "/mask"))

  # Write attributes
  write_h5ad_encoding(file, name, "nullable-integer", version)
}

#' Write H5AD string array
#'
#' Write a string array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_string_array <- function(value, file, name, version = "0.2.0") {
  requireNamespace("rhdf5")

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
#' Write a categorical from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_categorical <- function(value, file, name, version = "0.2.0") {
  rhdf5::h5createGroup(file, name)
  rhdf5::h5write(levels(value), file, paste0(name, "/categories"))
  rhdf5::h5write(as.integer(value) - 1, file, paste0(name, "/codes"))
  rhdf5::h5write(is.ordered(value), file, paste0(name, "/ordered"))

  write_h5ad_encoding(file, name, "categorical", version)
}

#' Write H5AD string scalar
#'
#' Write a string scalar from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_string_scalar <- function(value, file, name, version = "0.2.0") {
  requireNamespace("rhdf5")

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
#' Write a numeric scalar from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_numeric_scalar <- function(value, file, name, version = "0.2.0") {
  requireNamespace("rhdf5")

  # Write scalar
  rhdf5::h5write(value, file, name)

  # Write attributes
  write_h5ad_encoding(file, name, "numeric", version)
}

#' Write H5AD mapping
#'
#' Write a mapping from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_mapping <- function(value, file, name, version = "0.1.0") {
  rhdf5::h5createGroup(file, name)
  
  # Write mapping elements
  for (key in names(value)) {
    write_h5ad_element(value[[key]], file, paste0(name, "/", key))
  }

  write_h5ad_encoding(file, name, "dict", version)
}

#' Write H5AD data frame
#'
#' Write a data frame from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param index The index to write. Can either be a vector of length equal to
#' the number of rows in `values` or a single character string giving the name
#' of a column in `values`. If `NULL` then `rownames(value)` is used.
#' @param version Encoding version of the element to write
write_h5ad_data_frame <- function(value, file, name, index = NULL,
                                  version = "0.2.0") {
  
  rhdf5::h5createGroup(file, name)
  
  if (is.null(index)) {
    index_name <- "_index"
    value["_index"] <- rownames(value)
  } else if (length(index) == nrow(value)) {
    index_name <- "_index"
    value["_index"] <- index
  } else if (length(index) == 1 && index %in% colnames(value)) {
    index_name <- index
  } else {
    stop(
      "index must be a vector with length `nrow(value)` or a single character",
      "string giving the name of a column in `value`"
    )
  }
  
  # Write data frame columns
  for (col in colnames(value)) {
    write_h5ad_element(value[[col]], file, paste0(name, "/", col))
  }
  
  # Write encoding
  write_h5ad_encoding(file, name, "dataframe", version)
  
  # Write additional data frame attributes
  h5file <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(h5file))
  
  h5obj <- rhdf5::H5Gopen(h5file, name)
  on.exit(rhdf5::H5Gclose(h5obj), add = TRUE)
  
  rhdf5::h5writeAttribute(index_name, h5obj, "_index")
  col_order <- colnames(value)
  col_order <- col_order[col_order != index_name]
  rhdf5::h5writeAttribute(col_order, h5obj, "column-order") #nolint
}

#' HDF5 path exists
#'
#' Check that a path in HDF5 exists
#'
#' @param file Path to a HDF5 file
#' @param target_path The path within the file to test for
hdf5_path_exists <- function(file, target_path) {
  content <- rhdf5::h5ls(file)
  
  paths <- file.path(content$group, content$name)
  paths <- gsub("//", "/", paths) # Remove double slash for root paths
  
  target_path %in% paths
}
