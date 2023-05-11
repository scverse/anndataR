#' Set H5AD encoding
#'
#' Add the H5AD encoding to the attributes of an object
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param encoding The encoding type to set
#' @param version The encoding version to set
#'
#' @return The object with additional encoding attributes
write_h5ad_encoding <- function(file, name, encoding, version) {
  requireNamespace("rhdf5")

  dataset <- file & name
  rhdf5::h5writeAttribute(encoding, dataset, "encoding-type") # nolint
  rhdf5::h5writeAttribute(version, dataset, "encoding-version") # nolint
}

#' Write H5AD element
#'
#' Write an element to a H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
write_h5ad_element <- function(value, file, name) { # nolint
  write_fun <-
    if (is.matrix(value) || is.vector(value)) {
      write_h5ad_dense_array
    } else if (inherits(value, "sparseMatrix")) { # nolint
      write_h5ad_sparse_array
    } else if (is.factor(value)) {
      write_h5ad_categorical
    } else if (is.data.frame(value)) {
      write_h5ad_data_frame
    } else if (is.list(value)) {
      write_h5ad_mapping
    } else if (is.character(value) && length(value) == 1) {
      write_h5ad_string_scalar
    } else if (is.numeric(value) && length(value) == 1) {
      write_h5ad_numeric_scalar
    } else if (is.logical(value) && any(is.na(value))) {
      write_h5ad_nullable_boolean
    } else if (is.integer(value) && any(is.na(value))) {
      write_h5ad_nullable_integer
    } else {
      stop("Unsupported data type: ", class(value)) # nolint
    }
  write_fun(value = value, file = file, name = name)
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

path_exists <- function(file, target_path) {
  requireNamespace("rhdf5")

  content <- rhdf5::h5ls(file)

  any(paste0(content$group, content$name) == target_path)
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
  } else if (inherits(value, "RsparseMatrix")) {
    type <- "csc_matrix"
    indices_attr <- "i"
  } else {
    stop("Unsupported matrix format in ", name, ". Supported formats are RsparseMatrix and CsparseMatrix.")
  }

  if (path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

  # Write sparse matrix
  rhdf5::h5createGroup(file, name)
  rhdf5::h5write(attr(value, indices_attr), file, paste0(name, "/indices"))
  rhdf5::h5write(value@p, file, paste0(name, "/indptr"))
  rhdf5::h5write(value@x, file, paste0(name, "/data"))

  # add encoding
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

  # remove first if it already exists
  if (path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

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

  # remove first if it already exists
  if (path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

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

  # Write scalar
  rhdf5::h5write(value, file, name)

  # Write attributes
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
  requireNamespace("rhdf5")

  # remove first if it already exists
  if (path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

  # Write sparse matrix
  rhdf5::h5createGroup(file, name)
  rhdf5::h5write(levels(value), file, paste0(name, "/categories"))
  rhdf5::h5write(as.integer(value) - 1, file, paste0(name, "/codes"))
  rhdf5::h5write(inherits(value, "ordered"), file, paste0(name, "/ordered"))

  # add encoding
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
  rhdf5::h5write(value, file, name)

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
  requireNamespace("rhdf5")

  # delete name if it already exists?

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
#' @param version Encoding version of the element to write
write_h5ad_data_frame <- function(value, file, name, version = "0.2.0") {
  stop("Writing H5AD element not yet implemented")
}
