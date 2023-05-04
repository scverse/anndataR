#' Read H5AD encoding
#'
#' Read the encoding and version of an element in a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#'
#' @return A named list with the encoding and version
read_h5ad_encoding <- function(file, name) {
  requireNamespace("rhdf5")
  attrs <- rhdf5::h5readAttributes(file, name)

  if (!"encoding-type" %in% names(attrs) ||
    !"encoding-version" %in% names(attrs)) {
    stop(
      "Encoding not found for element '", name, "' in '", file, "'"
    )
  }

  list(
    encoding = attrs[["encoding-type"]],
    version = attrs[["encoding-version"]]
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
#' @details
#' If `encoding` is `NULL` the encoding and version are read from the element
#' using `read_h5ad_encoding()`
#'
#' @return Value depending on the encoding
read_h5ad_element <- function(file, name, encoding = NULL, version = NULL) {
  if (is.null(encoding)) {
    encoding_list <- read_h5ad_encoding(file, name)
    encoding <- encoding_list$encoding
    version <- encoding_list$version
  }

  switch(encoding,
    "array" = read_h5ad_dense_array(file, name, version = version),
    "csr_matrix" = read_h5ad_sparse_array(file, name,
      version = version,
      type = "csr"
    ),
    "csc_matrix" = read_h5ad_sparse_array(file, name,
      version = version,
      type = "csc"
    ),
    "dataframe" = read_h5ad_data_frame(file, name, version = version),
    "dict" = read_h5ad_mapping(file, name, version = version),
    "string" = read_h5ad_string_scalar(file, name, version = version),
    "numeric-scalar" = read_h5ad_numeric_scalar(file, name, version = version),
    "categorical" = read_h5ad_categorical(file, name, version = version),
    "string-array" = read_h5ad_string_array(file, name, version = version),
    "nullable-integer" = read_h5ad_nullable_integer(file, name,
      version = version
    ),
    "nullable-boolean" = read_h5ad_nullable_boolean(file, name,
      version = version
    ),
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
#' @return a matrix or a vector if 1D
read_h5ad_dense_array <- function(file, name, version = c("0.2.0")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  # TODO: ideally, native = TRUE should take care of the row order and column order,
  # but it doesn't
  darr <- t(rhdf5::h5read(file, name))
  # If the dense array is a 1D matrix, convert to vector
  if (any(dim(darr) == 1)) {
    darr <- as.vector(darr)
  }
  darr
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
#' @importFrom rhdf5 h5read h5readAttributes
#' @importFrom Matrix sparseMatrix
read_h5ad_sparse_array <- function(file, name, version = c("0.1.0"),
                                   type = c("csr", "csc")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  type <- match.arg(type)

  data <- rhdf5::h5read(file, paste0(name, "/data"))
  indices <- rhdf5::h5read(file, paste0(name, "/indices"))
  indptr <- rhdf5::h5read(file, paste0(name, "/indptr"))
  shape <- rhdf5::h5readAttributes(file, name)$shape

  if (type == "csc") {
    mtx <- sparseMatrix(i = indices, p = indptr, x = data, dims = shape, repr = "C", index1 = FALSE)
  } else {
    mtx <- sparseMatrix(j = indices, p = indptr, x = data, dims = shape, repr = "R", index1 = FALSE)
  }
  if (dim(mtx)[2] == 1) {
    mtx <- as.vector(mtx)
  }
  mtx
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
  requireNamespace("rhdf5")

  version <- match.arg(version)

  element <- rhdf5::h5read(file, name)

  # Get mask and convert to Boolean
  mask <- as.logical(element[["mask"]])

  # Get values and set missing
  element <- as.logical(element[["values"]])
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
  requireNamespace("rhdf5")

  version <- match.arg(version)

  element <- rhdf5::h5read(file, name)

  # Get mask and convert to Boolean
  mask <- as.logical(element[["mask"]])

  # Get values and set missing
  element <- as.integer(element[["values"]])
  element[mask] <- NA_integer_

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
  requireNamespace("rhdf5")

  version <- match.arg(version)

  rhdf5::h5read(file, name)
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
  requireNamespace("rhdf5")

  version <- match.arg(version)

  element <- rhdf5::h5read(file, name)

  # Get codes and convert to 1-based indexing
  codes <- element[["codes"]] + 1

  if (!length(dim(codes)) == 1) {
    stop("There is currently no support for multidimensional categorical arrays")
  }

  # Set missing values
  codes[codes == 0] <- NA

  levels <- element[["categories"]]

  ordered <- element[["ordered"]]
  if (is.null(ordered)) {
    # This version of {rhdf5} doesn't yet support ENUM type attributes so we
    # can't tell if the categorical should be ordered,
    # see https://github.com/grimbough/rhdf5/issues/125
    warning(
      "Unable to determine if categorical '", name,
      "' is ordered, assuming it isn't"
    )
    ordered <- FALSE
  }

  factor(levels[codes], levels = levels, ordered = ordered)
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
#' @importFrom rhdf5 h5read
read_h5ad_string_scalar <- function(file, name, version = c("0.2.0")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  rhdf5::h5read(file, name)
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
  requireNamespace("rhdf5")

  version <- match.arg(version)
  rhdf5::h5read(file, name)
}

#' Read H5AD mapping
#'
#' Read a mapping from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a named list
read_h5ad_mapping <- function(file, name, version = c("0.1.0")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  groupname <- paste0("/", name)
  columns <- subset(h5ls(file, recursive = T), group == groupname)$name
  
  read_h5ad_collection(file, name, columns)
}

# TODO: read index, return dataframe with columns

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
  requireNamespace("rhdf5")
  
  version <- match.arg(version)

  attributes <- rhdf5::h5readAttributes(file, name)
  index_name <- attributes$`_index`
  column_order <- attributes$`column-order`
  column_order <- append(column_order, index_name)
  
  columns <- read_h5ad_collection(file, name, column_order)
  
  index <- columns[[index_name]]
  columns[[index_name]] <- NULL
  
  if(length(columns) == 0){
    data.frame(row.names = index)
  } else {
    data.frame(columns, row.names = index)
  }
}

#' Read multiple H5AD datatypes
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#'
#' @return a named list
read_h5ad_collection <- function(file, name, column_order){
  columns <- list()
  for (col_name in column_order) {
    new_name <- file.path(name, col_name)
    encoding <- rhdf5::h5readAttributes(file, new_name)
    columns[[col_name]] <- read_h5ad_element(file, new_name, encoding$`encoding-type`, encoding$`encoding-version`)
  } 
  columns
}
