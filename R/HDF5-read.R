#' Read H5AD encoding
#'
#' Read the encoding and version of an element in a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#'
#' @return A named list with the encoding and version
read_hdf5_encoding <- function(file, path) {
  tryCatch(
    error = function(cnd){
      print(paste0("An error occurred when reading ", name, " in file ", file, "."))
      conditionMessage(cnd)
    },
    h5readAttributes(test_file, "X")
  )
}

#' Read H5AD element
#'
#' Read an element from a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param encoding The encoding of the element to read
#' @param version The encoding version of the element to read
#'
#' @return Value depending on the encoding
read_h5ad_element <- function(file, name, encoding, version) {
  switch (encoding,
    "array" = read_h5ad_dense_array(file, name, version = version),
    "csr_matrix" = read_h5ad_sparse_array(file, name, version = version,
                                          type = "csr"),
    "csc_matrix" = read_h5ad_sparse_array(file, name, version = version,
                                          type = "csc"),
    "dataframe" = read_h5ad_data_frame(file, name, version = version),
    "dict" = read_h5ad_mapping(file, name, version = version),
    "string" = read_h5ad_string_scalar(file, name, version = version),
    "numeric-scalar" = read_h5ad_numeric_scalar(file, name, version = version),
    "categorical" = read_h5ad_categorical(file, name, version = version),
    "string-array" = read_h5ad_string_array(file, name, version = version),
    "nullable-integer" = read_h5ad_nullable_integer(file, name,
                                                    version = version),
    "nullable-boolean" = read_h5ad_nullable_boolean(file, name,
                                                    version = version),
    stop("No function for reading H5AD encoding: ", encoding)
  )
}

#' Read H5AD dense array
#'
#' Read a dense array from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a Matrix/sparse matrix/DelayedArray???
read_hdf5_dense_array <- function(file, name, version = c("0.2.0")) {
  tryCatch(
    error = function(cnd){
      print(paste0("An error occurred when reading ", name, " in file."))
      conditionMessage(cnd)
    },
    {
      if(version == c("0.2.0")){
        h5read(file, name)
      }
    }
  )
}

#' Read H5AD sparse array
#'
#' Read a sparse array from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#' @param type Type of the sparse matrix, either "csr" or "csc"
#'
#' @return a sparse matrix/DelayedArray???
read_h5ad_sparse_array <- function(file, name, version = c("0.1.0"),
                                   type = c("csr", "csc")) {
  NULL
}

#' Read H5AD nullable boolean
#'
#' Read a nullable boolean from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a boolean vector
read_h5ad_nullable_boolean <- function(file, name, version = c("0.1.0")) {
  
  version <- match.arg(version)
  
  element <- rhdf5::h5read(file, name)
  
  # Get mask and convert to Boolean
  mask <- as.logical(element[["mask"]])
  
  # Get values and set missing
  element <- as.integer(element[["values"]])
  element[mask] <- NA
  
  return(element)
}

#' Read H5AD nullable integer
#'
#' Read a nullable integer from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return an integer vector
read_h5ad_nullable_integer <- function(file, name, version = c("0.1.0")) {
  
  version <- match.arg(version)
  
  element <- rhdf5::h5read(file, name)
  
  # Get mask and convert to Boolean
  mask <- as.logical(element[["mask"]])
  
  # Get values and set missing
  element <- as.logical(element[["values"]])
  element[mask] <- NA
  
  return(element)
}

#' Read H5AD string array
#'
#' Read a string array from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a character vector/matrix???
read_h5ad_string_array <- function(file, name, version = c("0.2.0")) {
  NULL
}

#' Read H5AD categorical
#'
#' Read a categorical from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a factor
read_h5ad_categorical <- function(file, name, version = c("0.2.0")) {
  
  version <- match.arg(version)
  
  element <- rhdf5::h5read(file, name)
  
  # Get codes and convert to 1-based indexing
  codes <- element[["codes"]] + 1
  
  if (!is.vector(codes)) {
    stop("There is currently no support for multidimensional categorical arrays")
  }
  
  # Set missing values
  codes[codes == 0] <- NA
  
  levels <- element[["categories"]]
  
  ordered <- element[["ordered"]]
  if (is.na(ordered)) {
    # This version of {rhdf5} doesn't yet support ENUM type attributes so we
    # can't tell if the categorical should be ordered, 
    # see https://github.com/grimbough/rhdf5/issues/125
    warning(
      "Unable to determine if categorical '", name,
      "' is ordered, assuming it isn't"
    )
    ordered <- FALSE
  }
  
  factor(levels[codes], levels=levels, ordered=ordered)
}

#' Read H5AD string scalar
#'
#' Read a string scalar from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a character vector of length 1
read_h5ad_string_scalar <- function(file, name, version = c("0.2.0")) {
  NULL
}

#' Read H5AD numeric scalar
#'
#' Read a numeric scalar from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a numeric vector of length 1
read_h5ad_numeric_scalar <- function(file, name, version = c("0.2.0")) {
  NULL
}

#' Read H5AD mapping
#'
#' Read a mapping from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a numeric vector of length 1
read_h5ad_mapping <- function(file, name, version = c("0.1.0")) {
  NULL
}

#' Read H5AD data frame
#'
#' Read a data frame from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a data.frame
read_h5ad_data_frame <- function(file, name, version = c("0.2.0")) {
  NULL
}

