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
    "nullable-integer" = read_h5ad_nullable_boolean(file, name,
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
#' @return a Matrix/sparse matrix/DelayedArray???, or a vector if 1D
read_hdf5_dense_array <- function(file, name, version = c("0.2.0")) {
  tryCatch(
    error = function(cnd){
      print(paste0("An error occurred when reading ", name, " in file."))
      conditionMessage(cnd)
    },
    {
      if(version == c("0.2.0")){
        # TODO: ideally, native = TRUE should take care of the row order and column order, 
        # but it doesn't
        darr <- t(h5read(file, name))
        # If the dense array is a 1D matrix, convert to vector
        if(dim(darr)[2] == 1){
          darr <- as.vector(darr)
        }
        darr
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
#' @return a sparse matrix/DelayedArray???, or a vector if 1D
read_h5ad_sparse_array <- function(file, name, version = c("0.1.0"),
                                   type = c("csr", "csc")) {
  tryCatch(
    error = function(cnd){
      conditionMessage(cnd)
    },
    if(version == c("0.1.0")){
      data <- h5read(file, paste0(name, "/data"))
      indices <- h5read(file, paste0(name, "/indices"))
      indptr <- h5read(file, paste0(name, "/indptr"))
      shape <- h5readAttributes(file, name)$shape
      
      if(type == "csc"){
        mtx <- sparseMatrix(i = indices, p = indptr, x = data, dims = shape, repr = "C", index1 = F)
      } else {
        mtx <- sparseMatrix(j = indices, p = indptr, x = data, dims = shape, repr = "R", index1 = F)
      }
      if(dim(mtx)[2] == 1){
        mtx <- as.vector(mtx)
      }
      mtx
    }
  )
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
  NULL
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
  NULL
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
  NULL
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

