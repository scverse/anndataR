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
#' `read_h5ad_encoding` and used to select the appropriate reading function.
#'
#' @return Value depending on the encoding
#'
#' @noRd
read_h5ad_element <- function(
  file,
  name,
  type = NULL,
  version = NULL,
  stop_on_error = FALSE,
  ...
) {
  if (!hdf5_path_exists(file, name)) {
    return(NULL)
  }

  if (is.null(type)) {
    encoding_list <- read_h5ad_encoding(file, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(
    type,
    "null" = read_h5ad_null,
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

#' Read H5AD element keys
#'
#' Read the keys of an element from a H5AD file
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
read_h5ad_element_keys <- function(
  file,
  name,
  type = NULL,
  version = NULL,
  stop_on_error = FALSE,
  ...
) {
  if (!hdf5_path_exists(file, name)) {
    return(NULL)
  }

  if (is.null(type)) {
    encoding_list <- read_h5ad_encoding(file, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(
    type,
    "dataframe" = read_h5ad_data_frame_keys,
    "dict" = read_h5ad_mapping_keys,
    cli_abort(
      "No function for reading keys for H5AD encoding {.cls {type}} for element {.val {name}}"
    )
  )

  tryCatch(
    {
      read_fun(file = file, name = name, version = version, ...)
    },
    error = function(e) {
      msg <- cli::cli_fmt(cli::cli_bullets(c(
        paste0(
          "Error reading element keys for {.field {name}} of type {.cls {type}}"
        ),
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

#' Read H5AD null
#'
#' Read a null value from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return `NULL`
#' @noRd
read_h5ad_null <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)

  NULL
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

  data <- rhdf5::h5read(file, name, native = FALSE)

  # If the array is 1D, explicitly add a dimension
  if (is.null(dim(data))) {
    data <- as.vector(data)
    dim(data) <- length(data)
  }

  # Transpose the matrix if need be
  if (is.matrix(data)) {
    data <- t(data)
  } else if (is.array(data) && length(dim(data)) > 1) {
    data <- aperm(data)
  }

  # Reverse {rhdf5} coercion to factors
  if (is.factor(data) && all(levels(data) %in% c("TRUE", "FALSE"))) {
    dims <- dim(data)
    data <- as.logical(data)
    dim(data) <- dims
  }

  data
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
read_h5ad_sparse_array <- function(
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

  x_data <- as.vector(h5group$data)
  # dgCMatrix/dgRMatrix x slot must be double
  if (!is.double(x_data)) {
    x_data <- as.double(x_data)
  }
  indices <- as.integer(as.vector(h5group$indices))
  indptr <- as.integer(as.vector(h5group$indptr))
  shape <- as.integer(as.vector(attrs[["shape"]]))

  # The Matrix package validity checks require that indices are sorted within
  # each major axis group (row indices within columns for CSC, column indices
  # within rows for CSR). For sparse matrices in Python order isn't guaranteed,
  # so we sort if needed.
  if (length(indices) > 1L) {
    row_lengths <- diff(indptr)
    group_ids <- rep.int(seq_along(row_lengths), row_lengths)
    ord <- order(group_ids, indices)
    if (is.unsorted(ord)) {
      indices <- indices[ord]
      x_data <- x_data[ord]
    }
  }

  if (type == "csc_matrix") {
    # Directly construct dgCMatrix (CSC format) to avoid overhead of constructing
    # a general sparseMatrix and then coercing to dgCMatrix
    # Slots: i = row indices (0-based), p = col pointers, x = values, Dim
    mtx <- new("dgCMatrix", i = indices, p = indptr, x = x_data, Dim = shape)
  } else if (type == "csr_matrix") {
    # Directly construct dgRMatrix (CSR format)
    # Slots: j = column indices (0-based), p = row pointers, x = values, Dim
    mtx <- new("dgRMatrix", j = indices, p = indptr, x = x_data, Dim = shape)
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
read_h5ad_string_array <- function(file, name, version = "0.2.0") {
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
read_h5ad_categorical <- function(file, name, version = "0.2.0") {
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
read_h5ad_string_scalar <- function(file, name, version = "0.2.0") {
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
read_h5ad_numeric_scalar <- function(file, name, version = "0.2.0") {
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
read_h5ad_mapping <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)

  items <- read_h5ad_mapping_keys(file, name, version)

  read_h5ad_collection(file, name, items)
}

#' Read H5AD mapping keys
#'
#' Read keys for a mapping from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#'
#' @return a character vector
#'
#' @noRd
read_h5ad_mapping_keys <- function(file, name, version = "0.1.0") {
  version <- match.arg(version)

  h5group <- rhdf5::H5Gopen(file, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)

  rhdf5::h5ls(h5group, recursive = FALSE)$name
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
read_h5ad_data_frame <- function(file, name, version = "0.2.0") {
  version <- match.arg(version)

  dim_keys <- read_h5ad_data_frame_keys(file, name, version)
  data <- read_h5ad_collection(file, name, dim_keys$cols)

  as.data.frame(
    row.names = dim_keys$rows,
    data,
    check.names = FALSE,
    fix.empty.names = FALSE
  )
}

#' Read H5AD data frame keys
#'
#' Read keys for a data frame from an H5AD file
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#' @param dim Dimension to read keys for, either "both, "rows" or "cols"
#'
#' @return a character vector if dim is "rows" or "cols", or a list with
#' elements "rows" and "cols" if dim is "both"
#'
#' @noRd
read_h5ad_data_frame_keys <- function(
  file,
  name,
  version = "0.2.0",
  dim = c("both", "rows", "cols")
) {
  version <- match.arg(version)
  dim <- match.arg(dim)

  attrs <- rhdf5::h5readAttributes(file, name, native = FALSE)
  index_name <- attrs[["_index"]]
  column_order <- attrs[["column-order"]]

  if (dim == "both") {
    list(
      rows = as.character(read_h5ad_element(file, file.path(name, index_name))),
      cols = as.character(column_order)
    )
  } else if (dim == "rows") {
    as.character(read_h5ad_element(file, file.path(name, index_name)))
  } else if (dim == "cols") {
    as.character(column_order)
  }
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
read_h5ad_collection <- function(file, name, item_names) {
  items <- lapply(
    item_names,
    function(item_name) {
      new_name <- paste0(name, "/", item_name)
      encoding <- read_h5ad_encoding(file, new_name)
      read_h5ad_element(
        file = file,
        name = new_name,
        type = encoding$type,
        version = encoding$version
      )
    }
  )
  names(items) <- item_names

  items
}
