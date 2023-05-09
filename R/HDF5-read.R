#' Read H5AD encoding
#'
#' Read the encoding and version of an element in a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#'
#' @return A named list with names type and version
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
    type = attrs[["encoding-type"]],
    version = attrs[["encoding-version"]]
  )
}

#' Read H5AD element
#'
#' Read an element from a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param type The encoding type of the element to read
#' @param version The encoding version of the element to read
#'
#' @details
#' If `encoding` is `NULL` the encoding and version are read from the element
#' using `read_h5ad_encoding()`
#'
#' @return Value depending on the encoding
read_h5ad_element <- function(file, name, type = NULL, version = NULL) {
  if (is.null(type)) {
    encoding_list <- read_h5ad_encoding(file, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <-
    switch(type,
      "array" = read_h5ad_dense_array,
      "csr_matrix" = read_h5ad_csr_matrix,
      "csc_matrix" = read_h5ad_csc_matrix,
      "dataframe" = read_h5ad_data_frame,
      "dict" = read_h5ad_mapping,
      "string" = read_h5ad_string_scalar,
      "numeric-scalar" = read_h5ad_numeric_scalar,
      "categorical" = read_h5ad_categorical,
      "string-array" = read_h5ad_string_array,
      "nullable-integer" = read_h5ad_nullable_integer,
      "nullable-boolean" = read_h5ad_nullable_boolean,
      stop("No function for reading H5AD encoding: ", type)
    )
  read_fun(file = file, name = name, version = version)
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
read_h5ad_dense_array <- function(file, name, version = "0.2.0") {
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

read_h5ad_csr_matrix <- function(file, name, version) {
  read_h5ad_sparse_array(
    file = file,
    name = name,
    version = version,
    type = "csr_matrix"
  )
}
read_h5ad_csc_matrix <- function(file, name, version) {
  read_h5ad_sparse_array(
    file = file,
    name = name,
    version = version,
    type = "csc_matrix"
  )
}

#' Read H5AD sparse array
#'
#' Read a sparse array from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#' @param type Type of the sparse matrix, either "csr_matrix" or "csc_matrix"
#'
#' @return a sparse matrix/DelayedArray???, or a vector if 1D
#' @importFrom Matrix sparseMatrix
read_h5ad_sparse_array <- function(file, name, version = "0.1.0",
                                   type = c("csr_matrix", "csc_matrix")) {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  type <- match.arg(type)

  data <- as.vector(rhdf5::h5read(file, paste0(name, "/data")))
  indices <- as.vector(rhdf5::h5read(file, paste0(name, "/indices")))
  indptr <- as.vector(rhdf5::h5read(file, paste0(name, "/indptr")))
  shape <- as.vector(rhdf5::h5readAttributes(file, name)$shape)

  if (type == "csc_matrix") {
    mtx <- Matrix::sparseMatrix(
      i = indices,
      p = indptr,
      x = data,
      dims = shape,
      repr = "C",
      index1 = FALSE
    )
  } else if (type == "csr_matrix") {
    mtx <- Matrix::sparseMatrix(
      j = indices,
      p = indptr,
      x = data,
      dims = shape,
      repr = "R",
      index1 = FALSE
    )
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
read_h5ad_nullable_boolean <- function(file, name, version = "0.1.0") {
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
read_h5ad_nullable_integer <- function(file, name, version = "0.1.0") {
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
read_h5ad_string_array <- function(file, name, version = "0.2.0") {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  # reads in transposed
  string_array <- rhdf5::h5read(file, name)
  if (is.matrix(string_array)) {
    string_array <- t(string_array)
  }
  string_array
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
read_h5ad_categorical <- function(file, name, version = "0.2.0") {
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

  factor(codes, labels = levels, ordered = ordered)
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
read_h5ad_string_scalar <- function(file, name, version = "0.2.0") {
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
read_h5ad_numeric_scalar <- function(file, name, version = "0.2.0") {
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
read_h5ad_mapping <- function(file, name, version = "0.1.0") {
  requireNamespace("rhdf5")

  version <- match.arg(version)
  groupname <- paste0("/", name)

  file_structure <- rhdf5::h5ls(file, recursive = TRUE)
  columns <- file_structure[file_structure$group == groupname, "name"]

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
read_h5ad_data_frame <- function(file, name, version = "0.2.0") {
  requireNamespace("rhdf5")

  version <- match.arg(version)

  attributes <- rhdf5::h5readAttributes(file, name)
  index_name <- attributes$`_index`
  column_order <- attributes$`column-order`
  column_order <- append(column_order, index_name)

  columns <- read_h5ad_collection(file, name, column_order)

  index <- columns[[index_name]]
  columns[[index_name]] <- NULL

  # NOTE: anndata will soon support non-character indices.
  # therefore we shouldn't set the index as row names of the data frame?
  if (length(columns) == 0) {
    data.frame(row.names = index)
  } else {
    data.frame(columns, row.names = index)
  }
}

#' Read multiple H5AD datatypes
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param column_order Vector of item names (in order)
#'
#' @return a named list
read_h5ad_collection <- function(file, name, column_order) {
  columns <- list()
  for (col_name in column_order) {
    new_name <- paste0(name, "/", col_name)
    encoding <- rhdf5::h5readAttributes(file, new_name)
    columns[[col_name]] <- read_h5ad_element(
      file = file,
      name = new_name,
      type = encoding$`encoding-type`,
      version = encoding$`encoding-version`
    )
  }
  columns
}
