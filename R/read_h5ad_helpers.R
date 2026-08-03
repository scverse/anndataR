#' Read H5AD element
#'
#' Read an element from a H5AD file
#'
#' @param hdf5_file An `HDF5File` object
#' @param name Name of the element within the H5AD file
#' @param type The encoding type of the element to read
#' @param version The encoding version of the element to read
#' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param backed Whether to return a DelayedArray for array/matrix types
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
  hdf5_file,
  name,
  type = NULL,
  version = NULL,
  stop_on_error = FALSE,
  backed = FALSE,
  ...
) {
  if (!hdf5_path_exists(hdf5_file, name)) {
    return(NULL)
  }

  if (is.null(type)) {
    encoding_list <- read_h5ad_encoding(hdf5_file, name)
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
    "nullable-string-array" = read_h5ad_nullable_string,
    cli_abort(
      "No function for reading H5AD encoding {.cls {type}} for element {.val {name}}"
    )
  )

  tryCatch(
    {
      read_fun(
        hdf5_file = hdf5_file,
        name = name,
        version = version,
        backed = backed,
        ...
      )
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
#' @inheritParams read_h5ad_element
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
  file$open_and_defer_close(readonly = TRUE)

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
      read_fun(hdf5_file = file, name = name, version = version, ...)
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

#' Read H5AD encoding
#'
#' Read the encoding and version of an element in a H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return A named list with names type and version
#'
#' @noRd
read_h5ad_encoding <- function(hdf5_file, name) {
  hdf5_file$open_and_defer_close(readonly = TRUE)

  tryCatch(
    {
      attrs <- rhdf5::h5readAttributes(hdf5_file$handle, name)
      list(
        type = attrs[["encoding-type"]],
        version = attrs[["encoding-version"]]
      )
    },
    error = function(e) {
      cli_abort(
        "Encoding attributes not found for element {.val {name}} in {.path {hdf5_file$path}}"
      )
    }
  )
}

#' Read H5AD null
#'
#' Read a null value from an H5AD file
#'
#' @param hdf5_file An `HDF5File` object
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#' @param ... Extra arguments for other reading functions (not used)
#'
#' @return `NULL`
#' @noRd
read_h5ad_null <- function(hdf5_file, name, version = "0.1.0", ...) {
  version <- match.arg(version)

  NULL
}

#' Read H5AD dense array
#'
#' Read a dense array from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a matrix/array or a DelayedArray if `backed = TRUE`
#'
#' @noRd
read_h5ad_dense_array <- function(
  hdf5_file,
  name,
  backed = FALSE,
  version = "0.2.0",
  ...
) {
  version <- match.arg(version)

  # Eager reads use base {rhdf5} directly
  if (!isTRUE(backed)) {
    return(read_h5ad_dense_array_base(hdf5_file, name, version))
  }

  dataset_type <- hdf5_get_dataset_type(hdf5_file, name)
  if (dataset_type == "H5T_INTEGER") {
    dataset_type <- "integer"
  } else if (dataset_type == "H5T_ENUM") {
    # HDF5Array doesn't support this type so fall back to base rhdf5
    return(read_h5ad_dense_array_base(hdf5_file, name, version))
  } else {
    dataset_type <- NA
  }

  # File handle must be closed for use by HDF5Array
  hdf5_file$close_and_defer_open()
  data <- HDF5Array::HDF5Array(hdf5_file$path, name, type = dataset_type)

  if (length(dim(data)) == 2) {
    data <- t(data)
  } else if (length(dim(data)) > 2) {
    data <- aperm(data)
  }

  data
}

#' Read H5AD dense array (base)
#'
#' Read a dense array from an H5AD file using base {rhdf5}
#'
#' @inheritParams read_h5ad_element
#'
#' @return a matrix/array
#'
#' @noRd
read_h5ad_dense_array_base <- function(
  hdf5_file,
  name,
  version = "0.2.0",
  ...
) {
  hdf5_file$open_and_defer_close(readonly = TRUE)
  data <- rhdf5::h5read(hdf5_file$handle, name, native = FALSE)

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

read_h5ad_csr_matrix <- function(hdf5_file, name, version, backed, ...) {
  read_h5ad_sparse_array(
    hdf5_file = hdf5_file,
    name = name,
    backed = backed,
    version = version,
    type = "csr_matrix"
  )
}

read_h5ad_csc_matrix <- function(hdf5_file, name, version, backed, ...) {
  read_h5ad_sparse_array(
    hdf5_file = hdf5_file,
    name = name,
    backed = backed,
    version = version,
    type = "csc_matrix"
  )
}

#' Read H5AD sparse array
#'
#' Read a sparse array from an H5AD file
#'
#' @inheritParams read_h5ad_element
#' @param type Type of the sparse matrix, either "csr_matrix" or "csc_matrix"
#'
#' @return a sparse matrix/array or a DelayedArray if `backed = TRUE`
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
read_h5ad_sparse_array <- function(
  hdf5_file,
  name,
  backed = FALSE,
  version = "0.1.0",
  type = c("csr_matrix", "csc_matrix"),
  ...
) {
  version <- match.arg(version)
  type <- match.arg(type)

  # Eager reads use base {rhdf5} directly
  if (!isTRUE(backed)) {
    return(read_h5ad_sparse_array_base(hdf5_file, name, version, type))
  }

  hdf5_file$close_and_defer_open()
  t(HDF5Array::H5SparseMatrix(hdf5_file$path, name))
}

#' Read H5AD sparse array (base)
#'
#' Read a sparse array from an H5AD file using base {rhdf5}
#'
#' @inheritParams read_h5ad_element
#' @inheritParams read_h5ad_sparse_array type
#'
#' @return a sparse matrix
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
read_h5ad_sparse_array_base <- function(
  hdf5_file,
  name,
  version = "0.1.0",
  type = c("csr_matrix", "csc_matrix"),
  ...
) {
  version <- match.arg(version)
  type <- match.arg(type)

  hdf5_file$open_and_defer_close(readonly = TRUE)
  attrs <- rhdf5::h5readAttributes(hdf5_file$handle, name, native = FALSE)

  h5group <- rhdf5::H5Gopen(hdf5_file$handle, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)

  construct_sparse_matrix(
    data = as.vector(h5group$data),
    indices = as.vector(h5group$indices),
    indptr = as.vector(h5group$indptr),
    shape = as.vector(attrs[["shape"]]),
    type = type
  )
}

#' Read H5AD recarray
#'
#' Read a recarray from an H5AD file
#'
#' @inheritParams read_h5ad_element
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
read_h5ad_rec_array <- function(hdf5_file, name, version = "0.2.0", ...) {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  rhdf5::h5read(
    hdf5_file$handle,
    name,
    native = FALSE,
    compoundAsDataFrame = FALSE
  ) |>
    lapply(as.vector)
}

#' Read H5AD nullable boolean
#'
#' Read a nullable boolean from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a boolean vector
#'
#' @noRd
read_h5ad_nullable_boolean <- function(
  hdf5_file,
  name,
  version = "0.1.0",
  ...
) {
  as.logical(read_h5ad_nullable(hdf5_file, name, version))
}

#' Read H5AD nullable integer
#'
#' Read a nullable integer from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return an integer vector
#'
#' @noRd
read_h5ad_nullable_integer <- function(
  hdf5_file,
  name,
  version = "0.1.0",
  ...
) {
  as.integer(read_h5ad_nullable(hdf5_file, name, version))
}

#' Read H5AD nullable string
#'
#' Read a nullable string from an H5AD file
#'
#' @param hdf5_file An `HDF5File` object
#' @param name Name of the element within the H5AD file
#' @param version Encoding version of the element to read
#' @param ... Extra arguments for other reading functions (not used)
#'
#' @details
#' **NOTE:** This implementation ignores the `"na-value"` attribute, see
#' https://anndata.scverse.org/en/stable/fileformat-prose.html#nullable-string-specifications-v0-1-0
#'
#' @return a character vector
#'
#' @noRd
read_h5ad_nullable_string <- function(
  hdf5_file,
  name,
  version = "0.1.0",
  ...
) {
  as.character(read_h5ad_nullable(hdf5_file, name, version))
}

#' Read H5AD nullable
#'
#' Read a nullable vector (boolean or integer) from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a nullable vector
#'
#' @noRd
read_h5ad_nullable <- function(hdf5_file, name, version = "0.1.0", ...) {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)
  h5group <- rhdf5::H5Gopen(hdf5_file$handle, name)
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
#' @inheritParams read_h5ad_element
#'
#' @return a character vector/matrix
#'
#' @noRd
read_h5ad_string_array <- function(hdf5_file, name, version = "0.2.0", ...) {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  data <- rhdf5::h5read(hdf5_file$handle, name, native = FALSE)

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
#' @inheritParams read_h5ad_element
#'
#' @return a factor
#'
#' @noRd
read_h5ad_categorical <- function(hdf5_file, name, version = "0.2.0", ...) {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  attrs <- rhdf5::h5readAttributes(hdf5_file$handle, name, native = FALSE)

  h5group <- rhdf5::H5Gopen(hdf5_file$handle, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)

  # Get codes and convert to 1-based indexing
  codes <- h5group$codes + 1L

  # Set missing values
  codes[codes == 0L] <- NA_integer_

  levels <- h5group$categories

  ordered <- attrs[["ordered"]]

  factor(levels[codes], levels = levels, ordered = ordered)
}

#' Read H5AD string scalar
#'
#' Read a string scalar from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a character vector of length 1
#'
#' @noRd
read_h5ad_string_scalar <- function(hdf5_file, name, version = "0.2.0", ...) {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  rhdf5::h5read(hdf5_file$handle, name, native = FALSE)
}

#' Read H5AD numeric scalar
#'
#' Read a numeric scalar from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a numeric vector of length 1
#'
#' @noRd
read_h5ad_numeric_scalar <- function(hdf5_file, name, version = "0.2.0", ...) {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  value <- rhdf5::h5read(hdf5_file$handle, name, native = FALSE)

  if (is.factor(value) && all(levels(value) %in% c("TRUE", "FALSE"))) {
    value <- as.logical(value)
  }

  value
}

#' Read H5AD mapping
#'
#' Read a mapping from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a named list
#'
#' @noRd
read_h5ad_mapping <- function(
  hdf5_file,
  name,
  version = "0.1.0",
  backed = FALSE,
  ...
) {
  version <- match.arg(version)

  items <- read_h5ad_mapping_keys(hdf5_file, name, version)

  read_h5ad_collection(hdf5_file, name, items, backed = backed, ...)
}

#' Read H5AD mapping keys
#'
#' Read keys for a mapping from an H5AD file
#'
#' @inheritParams read_h5ad_element
#'
#' @return a character vector
#'
#' @noRd
read_h5ad_mapping_keys <- function(hdf5_file, name, version = "0.1.0") {
  version <- match.arg(version)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  h5group <- rhdf5::H5Gopen(hdf5_file$handle, name)
  on.exit(rhdf5::H5Gclose(h5group), add = TRUE)

  rhdf5::h5ls(h5group, recursive = FALSE)$name
}

#' Read H5AD data frame
#'
#' Read a data frame from an H5AD file
#'
#' @inheritParams read_h5ad_element
#' @param backed Whether to return DelayedArrays for array/matrix types (always
#' `FALSE`)
#'
#' @details
#' The `backed` argument is included to avoid `backed` being passed via `...`
#' and is always overwritten to `FALSE` since data frame columns should not be
#' DelayedArrays.
#'
#' @return a data.frame
#'
#' @noRd
read_h5ad_data_frame <- function(
  hdf5_file,
  name,
  version = "0.2.0",
  backed = FALSE,
  ...
) {
  version <- match.arg(version)

  dim_keys <- read_h5ad_data_frame_keys(hdf5_file, name, version)
  data <- read_h5ad_collection(hdf5_file, name, dim_keys$cols)

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
#' @inheritParams read_h5ad_element
#' @param dim Dimension to read keys for, either "both, "rows" or "cols"
#'
#' @return a character vector if dim is "rows" or "cols", or a list with
#' elements "rows" and "cols" if dim is "both"
#'
#' @noRd
read_h5ad_data_frame_keys <- function(
  hdf5_file,
  name,
  version = "0.2.0",
  dim = c("both", "rows", "cols")
) {
  version <- match.arg(version)
  dim <- match.arg(dim)

  hdf5_file$open_and_defer_close(readonly = TRUE)

  attrs <- rhdf5::h5readAttributes(hdf5_file$handle, name, native = FALSE)
  index_name <- attrs[["_index"]]
  column_order <- attrs[["column-order"]]

  if (dim == "both") {
    list(
      rows = as.character(read_h5ad_element(
        hdf5_file,
        file.path(name, index_name)
      )),
      cols = as.character(column_order)
    )
  } else if (dim == "rows") {
    as.character(read_h5ad_element(hdf5_file, file.path(name, index_name)))
  } else if (dim == "cols") {
    as.character(column_order)
  }
}

#' Read multiple H5AD datatypes
#'
#' @inheritParams read_h5ad_element
#' @param item_names Vector of item names (in order)
#'
#' @return a named list
#'
#' @noRd
read_h5ad_collection <- function(
  hdf5_file,
  name,
  item_names,
  backed = FALSE,
  ...
) {
  hdf5_file$open_and_defer_close(readonly = TRUE)

  lapply(
    item_names,
    function(item_name) {
      new_name <- paste0(name, "/", item_name)
      encoding <- read_h5ad_encoding(hdf5_file, new_name)
      read_h5ad_element(
        hdf5_file = hdf5_file,
        name = new_name,
        type = encoding$type,
        version = encoding$version,
        backed = backed,
        ...
      )
    }
  ) |>
    setNames(item_names)
}
