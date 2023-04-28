#' Write H5AD encoding
#'
#' Write the encoding and version of an element to a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param encoding The encoding to write
#' @param version The encoding version to write
write_h5ad_encoding <- function(file, name, encoding, version) {
  NULL
}

#' Write H5AD element
#'
#' Write an element to a H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param encoding The encoding of the element to read
#' @param version The encoding version of the element to read
write_h5ad_element <- function(value, file, name, encoding, version) {
  switch (encoding,
    "array" = write_h5ad_dense_array(file, name, value = value,
                                     version = version),
    "csr_matrix" = write_h5ad_sparse_array(file, name, value = value,
                                           version = version, type = "csr"),
    "csc_matrix" = write_h5ad_sparse_array(file, name, value = value,
                                           version = version, type = "csc"),
    "dataframe" = write_h5ad_data_frame(file, name, value = value, 
                                        version = version),
    "dict" = write_h5ad_mapping(file, name, value = value, 
                                version = version),
    "string" = write_h5ad_string_scalar(file, name, value = value, 
                                        version = version),
    "numeric-scalar" = write_h5ad_numeric_scalar(file, name, value = value,
                                                 version = version),
    "categorical" = write_h5ad_categorical(file, name, value = value,
                                           version = version),
    "string-array" = write_h5ad_string_array(file, name, value = value,
                                             version = version),
    "nullable-integer" = write_h5ad_nullable_integer(file, name, value = value,
                                                     version = version),
    "nullable-integer" = write_h5ad_nullable_boolean(file, name, value = value,
                                                     version = version),
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
  NULL
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
  NULL
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
  NULL
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
  NULL
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
  NULL
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
  NULL
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
  NULL
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
  NULL
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
  NULL
}
