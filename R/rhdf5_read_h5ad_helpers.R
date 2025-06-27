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
rhdf5_read_h5ad_encoding <- function(file, name) {
  tryCatch(
    {
      attrs <- rhdf5::h5readAttributes(file, name)
      list(
        type = attrs[["encoding-type"]],
        version = attrs[["encoding-version"]]
      )
    },
    error = function(e) {
      path <- if (is.character(file)) file else rhdf5::H5Fget_name(file) # nolint object_usage_linter
      cli_abort(
        "Encoding attributes not found for element {.val {name}} in {.path {path}}"
      )
    }
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
#' `rhdf5_read_h5ad_encoding` and used to select the appropriate reading function.
#'
#' @return Value depending on the encoding
#'
#' @noRd
rhdf5_read_h5ad_element <- function(
  file,
  name,
  type = NULL,
  version = NULL,
  stop_on_error = FALSE,
  ...
) {
  if (!rhdf5_hdf5_path_exists(file, name)) {
    return(NULL)
  }

  if (is.null(type)) {
    encoding_list <- rhdf5_read_h5ad_encoding(file, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(
    type,
    "array" = rhdf5_read_h5ad_dense_array,
    "rec-array" = rhdf5_read_h5ad_rec_array,
    "csr_matrix" = rhdf5_read_h5ad_csr_matrix,
    "csc_matrix" = rhdf5_read_h5ad_csc_matrix,
    "dataframe" = rhdf5_read_h5ad_data_frame,
    "dict" = rhdf5_read_h5ad_mapping,
    "string" = rhdf5_read_h5ad_string_scalar,
    "numeric-scalar" = rhdf5_read_h5ad_numeric_scalar,
    "categorical" = rhdf5_read_h5ad_categorical,
    "string-array" = rhdf5_read_h5ad_string_array,
    "nullable-integer" = rhdf5_read_h5ad_nullable_integer,
    "nullable-boolean" = rhdf5_read_h5ad_nullable_boolean,
    cli_abort(
      "No function for reading H5AD encoding {.cls {type}} for element {.val {name}}"
    )
  )

  tryCatch(
    {
      read_fun(file = file, name = name, version = version, ...)
    },
    error = function(e) {
      msg <- cli::cli_fmt(cli::cli_bullets(c(
        paste0("Error reading element {.field {name}} of type {.cls {type}}"),
        "i" = conditionMessage(e)
      )))
      if (stop_on_error) {
        cli_abort(msg)
      } else {
        cli_warn(msg)
        NULL
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
rhdf5_read_h5ad_dense_array <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  data <- rhdf5::h5read(file, name, native = FALSE)

  # If the array is 1D, explicitly add a dimension
  if (is.null(dim(data))) {
    data <- as.vector(data)
    dim(data) <- length(data)
  }

  # transpose the matrix if need be
  if (is.matrix(data)) {
    data <- t(data)
  } else if (is.array(data) && length(dim(data)) > 1) {
    data <- aperm(data)
  }

  # Reverse {rhdf5} coercion to factors
  if (is.factor(data) && all(levels(data) %in% c("TRUE", "FALSE"))) {
    data <- as.logical(data)
    dim(data) <- length(data)
  }

  data
}

rhdf5_read_h5ad_csr_matrix <- function(file, name, version) {
  rhdf5_read_h5ad_sparse_array(
    file = file,
    name = name,
    version = version,
    type = "csr_matrix"
  )
}

rhdf5_read_h5ad_csc_matrix <- function(file, name, version) {
  rhdf5_read_h5ad_sparse_array(
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
rhdf5_read_h5ad_sparse_array <- function(
  file,
  name,
  version = "0.1.0",
  type = c("csr_matrix", "csc_matrix")
) {
  version <- match.arg(version)
  type <- match.arg(type)

  h5group <- rhdf5::H5Gopen(file, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)
  attrs <- rhdf5::h5readAttributes(file, name, native = FALSE)

  data <- as.vector(h5group$data)
  indices <- as.vector(h5group$indices)
  indptr <- as.vector(h5group$indptr)
  shape <- as.vector(attrs[["shape"]])

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
rhdf5_read_h5ad_rec_array <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  rhdf5::h5read(file, name, native = FALSE, compoundAsDataFrame = FALSE) |>
    lapply(as.vector)
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
# nolint start: object_length_linter
rhdf5_read_h5ad_nullable_boolean <- function(file, name, version = "0.1.0") {
  # nolint end: object_length_linter
  as.logical(rhdf5_read_h5ad_nullable(file, name, version))
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
# nolint start: object_length_linter
rhdf5_read_h5ad_nullable_integer <- function(file, name, version = "0.1.0") {
  # nolint end: object_length_linter
  as.integer(rhdf5_read_h5ad_nullable(file, name, version))
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
rhdf5_read_h5ad_nullable <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)

  h5group <- rhdf5::H5Gopen(file, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)

  data <- as.vector(h5group$values)

  mask <- as.logical(h5group$mask)

  data[mask] <- NA

  data
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
rhdf5_read_h5ad_string_array <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  data <- rhdf5::h5read(file, name, native = FALSE)

  if (is.null(dim(data)) || length(dim(data)) == 1) {
    data <- as.vector(data)
    dim(data) <- length(data)
  }

  # transpose the matrix if need be
  if (is.matrix(data)) {
    data <- t(data)
  } else if (is.array(data) && length(dim(data)) > 1) {
    data <- aperm(data)
  }

  data
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
rhdf5_read_h5ad_categorical <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  h5group <- rhdf5::H5Gopen(file, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)

  # Get codes and convert to 1-based indexing
  codes <- h5group$codes + 1L

  # Set missing values
  codes[codes == 0L] <- NA_integer_

  levels <- h5group$categories

  attrs <- rhdf5::h5readAttributes(file, name, native = FALSE)
  ordered <- attrs[["ordered"]]

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
#'
#' @noRd
rhdf5_read_h5ad_string_scalar <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  rhdf5::h5read(file, name, native = FALSE)
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
rhdf5_read_h5ad_numeric_scalar <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  value <- rhdf5::h5read(file, name, native = FALSE)

  if (is.factor(value) && all(levels(value) %in% c("TRUE", "FALSE"))) {
    value <- as.logical(value)
  }

  value
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
rhdf5_read_h5ad_mapping <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)

  h5group <- rhdf5::H5Gopen(file, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)
  items <- rhdf5::h5ls(h5group, recursive = FALSE)$name

  rhdf5_read_h5ad_collection(file, name, items)
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
#'
#' @noRd
rhdf5_read_h5ad_data_frame <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  attrs <- rhdf5::h5readAttributes(file, name, native = FALSE)
  index_name <- attrs[["_index"]]
  column_order <- attrs[["column-order"]]

  index <- rhdf5_read_h5ad_element(file, file.path(name, index_name))
  data <- rhdf5_read_h5ad_collection(file, name, column_order)

  as.data.frame(
    row.names = index,
    data,
    check.names = FALSE,
    fix.empty.names = FALSE
  )
}

#' Read multiple H5AD datatypes
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param item_names Vector of item names (in order)
#'
#' @return a named list
#'
#' @noRd
rhdf5_read_h5ad_collection <- function(file, name, item_names) {
  items <- list()
  for (item_name in item_names) {
    new_name <- paste0(name, "/", item_name)
    encoding <- rhdf5_read_h5ad_encoding(file, new_name)
    items[[item_name]] <- rhdf5_read_h5ad_element(
      file = file,
      name = new_name,
      type = encoding$type,
      version = encoding$version
    )
  }
  items
}
