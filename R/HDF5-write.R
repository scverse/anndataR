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
  attributes(x) <- list(
    dim = attr(x, "dim"),
    "encoding-type" = encoding,
    "encoding-version" = version
  )

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
  version <- match.arg(version)

  value <- set_h5ad_encoding(value, encoding = "array", version = version)

  # Transpose the value because writing with native=TRUE doesn't seem to work
  # as expected
  rhdf5::h5write(t(value), file, name, write.attributes = TRUE)
}

#' Write H5AD sparse array
#'
#' Write a sparse array from an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to write
write_h5ad_sparse_array <- function(value, file, name, version = c("0.1.0"),
                                    type = c("csr", "csc")) {
  stop("Writing H5AD element not yet implemented")
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
