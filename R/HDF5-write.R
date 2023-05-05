#' Set H5AD encoding
#'
#' Add the H5AD encoding to the attributes of an object
#'
#' @param x Object to set encoding on
#' @param encoding The encoding type to set
#' @param version The encoding version to set
#'
#' @return The object with additional encoding attributes
set_h5ad_encoding <- function(x, encoding, version) {
  # nolint start
  attr(x, "encoding-type") <- encoding
  attr(x, "encoding-version") <- version
  # nolint end

  return(x)
}

#' Write H5AD element
#'
#' Write an element to a H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
write_h5ad_element <- function(value, file, name) {
  encoding <- switch(class(value)[1],
    "matrix" = "array",
    stop("Unknown encoding for class '", class(value)[1], "'")
  )

  switch(encoding,
    "array" = write_h5ad_dense_array(file, name, value = value),
    "csr_matrix" = write_h5ad_sparse_array(file, name,
      value = value,
      type = "csr"
    ),
    "csc_matrix" = write_h5ad_sparse_array(file, name,
      value = value,
      type = "csc"
    ),
    "dataframe" = write_h5ad_data_frame(file, name, value = value),
    "dict" = write_h5ad_mapping(file, name, value = value),
    "string" = write_h5ad_string_scalar(file, name, value = value),
    "numeric-scalar" = write_h5ad_numeric_scalar(file, name, value = value),
    "categorical" = write_h5ad_categorical(file, name, value = value),
    "string-array" = write_h5ad_string_array(file, name, value = value),
    "nullable-integer" = write_h5ad_nullable_integer(file, name, value = value),
    "nullable-integer" = write_h5ad_nullable_boolean(file, name, value = value),
    stop("No function for writing H5AD encoding: ", encoding)
  )
}

#' Write H5AD dense array
#'
#' Write a dense array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_dense_array <- function(value, file, name, version = c("0.2.0")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)

  # Transpose the value because writing with native=TRUE doesn't seem to work
  # as expected
  value <- t(value)

  value <- set_h5ad_encoding(value, encoding = "array", version = version)

  rhdf5::h5write(value, file, name, write.attributes = TRUE)
}

path_exists <- function(file, target_path) {
  requireNamespace("rhdf5")
  content <- rhdf5::h5ls(file)
  return(any(content$path == target_path))
}

#' Write H5AD sparse array
#'
#' Write a sparse array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
#'
#' @examples
#' value <- Matrix::rsparsematrix(10, 12, .1)
#' file_path <- system.file("extdata", "krumsiek11_augmented_sparse_v0-8.h5ad", package = "anndataR")
#' file <- rhdf5::H5Fopen(file_path)
#' name <- "/X"
#' write_h5ad_sparse_array(value, file, name)
write_h5ad_sparse_array <- function(value, file, name, version = c("0.1.0"),
                                    type = c("csr", "csc")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  type <- match.arg(type)

  stopifnot(inherits(value, "sparseMatrix"))

  if (type == "csr") {
    value <- as(value, "RsparseMatrix")
    indices_attr <- "j"
  } else if (type == "csc") {
    value <- as(value, "CsparseMatrix")
    indices_attr <- "i"
  }

  out <- list(
    data = value@x,
    indices = attr(value, indices_attr),
    indptr = value@p
  )
  out <- set_h5ad_encoding(out, encoding = type, version = version)

  if (path_exists(file, name)) {
    rhdf5::h5delete(file, name)
  }

  rhdf5::h5write(out, file, name, write.attributes = TRUE)
}

#' Write H5AD nullable boolean
#'
#' Write a nullable boolean from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_nullable_boolean <- function(value, file, name,
                                        version = c("0.1.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD nullable integer
#'
#' Write a nullable integer from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_nullable_integer <- function(value, file, name,
                                        version = c("0.1.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD string array
#'
#' Write a string array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_string_array <- function(value, file, name, version = c("0.2.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD categorical
#'
#' Write a categorical from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_categorical <- function(value, file, name, version = c("0.2.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD string scalar
#'
#' Write a string scalar from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_string_scalar <- function(value, file, name, version = c("0.2.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD numeric scalar
#'
#' Write a numeric scalar from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_numeric_scalar <- function(value, file, name, version = c("0.2.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD mapping
#'
#' Write a mapping from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_mapping <- function(value, file, name, version = c("0.1.0")) {
  stop("Writing H5AD element not yet implemented")
}

#' Write H5AD data frame
#'
#' Write a data frame from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_data_frame <- function(value, file, name, version = c("0.2.0")) {
  stop("Writing H5AD element not yet implemented")
}
