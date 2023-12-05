#' Read H5AD encoding
#'
#' Read the encoding and version of an element in a H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#'
#' @return A named list with names type and version
#'
#' @noRd
read_h5ad_encoding <- function(file, name) {
  attrs <- rhdf5::h5readAttributes(file, name)

  if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
    path <- if (is.character(file)) file else rhdf5::H5Fget_name(file)
    stop(
      "Encoding attributes not found for element '", name, "' ",
      "in '", path, "'"
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
#' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param ... Extra arguments passed to individual reading functions
#'
#' @details
#' Encoding is automatically determined from the element using
#' `read_h5ad_encoding` and used to select the appropriate reading function.
#'
#' @return Value depending on the encoding
#'
#' @noRd
read_h5ad_element <- function(file, name, type = NULL, version = NULL, stop_on_error = FALSE, ...) {
  if (is.null(type)) {
    encoding_list <- read_h5ad_encoding(file, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(type,
    "array" = read_h5ad_dense_array,
    "rec-array" = read_h5ad_rec_array,
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
    stop(
      "No function for reading H5AD encoding '", type,
      "' for element '", name, "'"
    )
  )

  tryCatch(
    {
      read_fun(file = file, name = name, version = version, ...)
    },
    error = function(e) {
      message <- paste0(
        "Error reading element '", name, "' of type '", type, "':\n",
        conditionMessage(e)
      )
      if (stop_on_error) {
        stop(message)
      } else {
        warning(message)
        return(NULL)
      }
    }
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
#'
#' @noRd
read_h5ad_dense_array <- function(file, name, version = "0.2.0") {
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
#'
#' @noRd
read_h5ad_sparse_array <- function(file, name, version = "0.1.0",
                                   type = c("csr_matrix", "csc_matrix")) {
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

#' Read H5AD recarray
#'
#' Read a recarray from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @details
#' A "record array" (recarray) is a Python NumPy array type that contains
#' "fields" that can be indexed using attributes (similar to columns in a
#' spreadsheet). See https://numpy.org/doc/stable/reference/generated/numpy.recarray.html
#' for details.
#'
#' They are used by **scanpy** to score marker gene testing results.
#'
#' @return a named list of 1D arrays
#'
#' @noRd
read_h5ad_rec_array <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  rhdf5::h5read(file, name, compoundAsDataFrame = FALSE)
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
#' 
#' @noRd
read_h5ad_nullable_boolean <- function(file, name, version = "0.1.0") {
  as.logical(read_h5ad_nullable(file, name, version))
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
#'
#' @noRd
read_h5ad_nullable_integer <- function(file, name, version = "0.1.0") {
  as.integer(read_h5ad_nullable(file, name, version))
}

#' Read H5AD nullable
#'
#' Read a nullable vector (boolean or integer) from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a nullable vector
#'
#' @noRd
read_h5ad_nullable <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)

  element <- rhdf5::h5read(file, name)

  # Some versions of rhdf5 automatically apply mask, in which case
  # there is no 'mask' element
  if (!is.null(names(element))) {
    # Get mask and convert to Boolean
    mask <- as.logical(element[["mask"]])
    # Get values and set missing
    element <- as.vector(element[["values"]])
    element[mask] <- NA
  }

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
#' @return a character vector/matrix
#'
#' @noRd
read_h5ad_string_array <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)
  # reads in transposed
  string_array <- rhdf5::h5read(file, name)
  if (is.matrix(string_array)) {
    string_array <- t(string_array)
  }

  # If the array is 1D, convert to vector
  if (length(dim(string_array)) == 1) {
    string_array <- as.vector(string_array)
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
#'
#' @noRd
read_h5ad_categorical <- function(file, name, version = "0.2.0") {
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

  attributes <- rhdf5::h5readAttributes(file, name)
  ordered <- attributes[["ordered"]]
  if (is.null(ordered) || is.na(ordered)) {
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
#'
#' @noRd
read_h5ad_string_scalar <- function(file, name, version = "0.2.0") {
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
#'
#' @noRd
read_h5ad_numeric_scalar <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)
  scalar <- rhdf5::h5read(file, name)

  # If the numeric vector is Boolean it gets read as a factor by {rhdf5}
  if (is.factor(scalar)) {
    scalar <- as.logical(scalar)
  }

  return(scalar)
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
#'
#' @noRd
read_h5ad_mapping <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)
  groupname <- paste0("/", name)

  file_structure <- rhdf5::h5ls(file, recursive = TRUE)
  columns <- file_structure[file_structure$group == groupname, "name"]

  read_h5ad_collection(file, name, columns)
}

#' Read H5AD data frame
#'
#' Read a data frame from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#' @param include_index Whether or not to include the index as a column
#'
#' @details
#' If `include_index == TRUE` the index stored in the HDF5 file is added as a
#' column to output `data.frame` using the defined index name as the column
#' name and this is set as an attribute. If `include_index == FALSE` the index
#' is not provided in the output. In either case row names are not set.
#'
#' @return a data.frame
#'
#' @noRd
read_h5ad_data_frame <- function(file, name, include_index = FALSE,
                                 version = "0.2.0") {
  version <- match.arg(version)

  attributes <- rhdf5::h5readAttributes(file, name)
  index_name <- attributes$`_index`
  column_order <- attributes$`column-order`

  columns <- read_h5ad_collection(file, name, column_order)

  if (length(columns) == 0) {
    index <- read_h5ad_data_frame_index(file, name)
    df <- data.frame(row.names = seq_along(index))
  } else {
    df <- data.frame(columns)
  }

  if (isTRUE(include_index)) {
    index <- read_h5ad_data_frame_index(file, name)
    df <- cbind(index, df)

    # The default index name is not allowed as a column name so adjust it
    if (index_name == "_index") {
      index_name <- ".index"
      colnames(df)[1] <- index_name
    }

    attr(df, "_index") <- index_name # nolint
  }

  df
}

#' Read H5AD data frame index
#'
#' Read the index of a data frame from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return an object containing the index
#'
#' @noRd
read_h5ad_data_frame_index <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  attributes <- rhdf5::h5readAttributes(file, name)
  index_name <- attributes$`_index`

  read_h5ad_element(file, file.path(name, index_name))
}

#' Read multiple H5AD datatypes
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param column_order Vector of item names (in order)
#'
#' @return a named list
#'
#' @noRd
read_h5ad_collection <- function(file, name, column_order) {
  columns <- list()
  for (col_name in column_order) {
    new_name <- paste0(name, "/", col_name)
    encoding <- read_h5ad_encoding(file, new_name)
    columns[[col_name]] <- read_h5ad_element(
      file = file,
      name = new_name,
      type = encoding$type,
      version = encoding$version
    )
  }
  columns
}
