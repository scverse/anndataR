#' Read Zarr encoding
#'
#' Read the encoding and version of an element in a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#'
#' @return A named list with names type and version
#'
#' @noRd
read_zarr_encoding <- function(store, name, stop_on_error = TRUE) {
  # Path can be to array or group
  g <- pizzarr::zarr_open(store, path = name)
  attrs <- g$get_attrs()$to_list()

  if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
    if (stop_on_error) {
      stop(
        "Encoding attributes not found for element '", name, "' "
      )
    } else {
      return(NULL)
    }
  }

  list(
    type = attrs[["encoding-type"]],
    version = attrs[["encoding-version"]]
  )
}

#' Read Zarr element
#'
#' Read an element from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param type The encoding type of the element to read
#' @param version The encoding version of the element to read
#' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param ... Extra arguments passed to individual reading functions
#'
#' @details
#' Encoding is automatically determined from the element using
#' `read_zarr_encoding` and used to select the appropriate reading function.
#'
#' @return Value depending on the encoding
#'
#' @noRd
read_zarr_element <- function(store, name, type = NULL, version = NULL, stop_on_error = FALSE, ...) {
  if (is.null(type)) {
    encoding_list <- read_zarr_encoding(store, name, stop_on_error = stop_on_error)
    if (is.null(encoding_list)) {
      if (stop_on_error) {
        stop("No encoding info found for element '", name, "'")
      } else {
        warning("No encoding found for element '", name, "'")
        return(NULL)
      }
    }
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(type,
    "array" = read_zarr_dense_array,
    "rec-array" = read_zarr_rec_array,
    "csr_matrix" = read_zarr_csr_matrix,
    "csc_matrix" = read_zarr_csc_matrix,
    "dataframe" = read_zarr_data_frame,
    "dict" = read_zarr_mapping,
    "string" = read_zarr_string_scalar,
    "numeric-scalar" = read_zarr_numeric_scalar,
    "categorical" = read_zarr_categorical,
    "string-array" = read_zarr_string_array,
    "nullable-integer" = read_zarr_nullable_integer,
    "nullable-boolean" = read_zarr_nullable_boolean,
    stop(
      "No function for reading H5AD encoding '", type,
      "' for element '", name, "'"
    )
  )

  tryCatch(
    {
      read_fun(store = store, name = name, version = version, ...)
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

read_zarr_array <- function(store, name) {
  zarr_arr <- pizzarr::zarr_open_array(store, path = name)
  nested_arr <- zarr_arr$get_item("...")
  return(nested_arr$data)
}

#' Read Zarr dense array
#'
#' Read a dense array from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a matrix or a vector if 1D
#'
#' @noRd
read_zarr_dense_array <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  # Extract the NestedArray contents as a base R array.
  darr <- read_zarr_array(store, name)


  # TODO: ideally, native = TRUE should take care of the row order and column order,
  # but it doesn't
  # If the dense array is a 1D matrix, convert to vector
  if (length(dim(darr)) == 1) {
    darr <- as.vector(darr)
  }
  darr
}

read_zarr_csr_matrix <- function(store, name, version) {
  read_zarr_sparse_array(
    store = store,
    name = name,
    version = version,
    type = "csr_matrix"
  )
}

read_zarr_csc_matrix <- function(store, name, version) {
  read_zarr_sparse_array(
    store = store,
    name = name,
    version = version,
    type = "csc_matrix"
  )
}

#' Read Zarr sparse array
#'
#' Read a sparse array from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#' @param type Type of the sparse matrix, either "csr_matrix" or "csc_matrix"
#'
#' @return a sparse matrix/DelayedArray???, or a vector if 1D
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
read_zarr_sparse_array <- function(store, name, version = "0.1.0",
                                   type = c("csr_matrix", "csc_matrix")) {
  version <- match.arg(version)
  type <- match.arg(type)

  g <- pizzarr::zarr_open_group(store, path = name)

  data <- as.vector(read_zarr_array(store, paste0(name, "/data")))
  indices <- as.vector(read_zarr_array(store, paste0(name, "/indices")))
  indptr <- as.vector(read_zarr_array(store, paste0(name, "/indptr")))
  shape <- as.vector(unlist(g$get_attrs()$to_list()$shape, use.names = FALSE))

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

#' Read Zarr recarray
#'
#' Read a recarray from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
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
read_zarr_rec_array <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  stop("Reading recarrays is not yet implemented")
}

#' Read Zarr nullable boolean
#'
#' Read a nullable boolean from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a boolean vector
#'
#' @noRd
read_zarr_nullable_boolean <- function(store, name, version = "0.1.0") {
  as.logical(read_zarr_nullable(store, name, version))
}

#' Read Zarr nullable integer
#'
#' Read a nullable integer from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return an integer vector
#'
#' @noRd
read_zarr_nullable_integer <- function(store, name, version = "0.1.0") {
  as.integer(read_zarr_nullable(store, name, version))
}

#' Read Zarr nullable
#'
#' Read a nullable vector (boolean or integer) from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a nullable vector
#'
#' @noRd
read_zarr_nullable <- function(store, name, version = "0.1.0") {
  version <- match.arg(version)

  mask <- read_zarr_array(store, paste0(name, "/mask"))
  values <- read_zarr_array(store, paste0(name, "/values"))

  # Get values and set missing
  element <- values
  element[mask] <- NA

  return(element)
}

#' Read Zarr string array
#'
#' Read a string array from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a character vector/matrix
#'
#' @noRd
read_zarr_string_array <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)
  # reads in transposed
  string_array <- read_zarr_array(store, name)

  # If the array is 1D, convert to vector
  if (length(dim(string_array)) == 1) {
    string_array <- as.vector(string_array)
  }

  string_array
}

#' Read Zarr categorical
#'
#' Read a categorical from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a factor
#'
#' @noRd
read_zarr_categorical <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  codes <- read_zarr_array(store, paste0(name, "/codes"))
  categories <- read_zarr_array(store, paste0(name, "/categories"))

  # Get codes and convert to 1-based indexing
  codes <- codes + 1

  if (!length(dim(codes)) == 1) {
    stop("There is currently no support for multidimensional categorical arrays")
  }

  # Set missing values
  codes[codes == 0] <- NA

  levels <- categories

  g <- pizzarr::zarr_open_group(store, path = name)

  attributes <- g$get_attrs()$to_list()
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

#' Read Zarr string scalar
#'
#' Read a string scalar from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a character vector of length 1
#'
#' @noRd
read_zarr_string_scalar <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)
  scalar <- as.character(read_zarr_array(store, name))
  return(scalar)
}

#' Read Zarr numeric scalar
#'
#' Read a numeric scalar from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a numeric vector of length 1
#'
#' @noRd
read_zarr_numeric_scalar <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)
  scalar <- as.numeric(read_zarr_array(store, name))
  return(scalar)
}

#' Read Zarr mapping
#'
#' Read a mapping from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return a named list
#'
#' @noRd
read_zarr_mapping <- function(store, name, version = "0.1.0") {
  version <- match.arg(version)

  g <- pizzarr::zarr_open(store)
  columns <- g$get_store()$listdir(name)

  # Omit Zarr metadata files from the list of columns.
  columns <- columns[!columns %in% c(".zgroup", ".zattrs", ".zarray")]

  read_zarr_collection(store, name, columns)
}

#' Read Zarr data frame
#'
#' Read a data frame from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#' @param include_index Whether or not to include the index as a column
#'
#' @details
#' If `include_index == TRUE` the index stored in the Zarr store is added as a
#' column to output `data.frame` using the defined index name as the column
#' name and this is set as an attribute. If `include_index == FALSE` the index
#' is not provided in the output. In either case row names are not set.
#'
#' @return a data.frame
#'
#' @noRd
read_zarr_data_frame <- function(store, name, include_index = TRUE,
                                 version = "0.2.0") {
  version <- match.arg(version)

  g <- pizzarr::zarr_open_group(store, path = name)

  attributes <- g$get_attrs()$to_list()
  index_name <- attributes$`_index`
  column_order <- attributes$`column-order`

  columns <- read_zarr_collection(store, name, column_order)

  if (length(columns) == 0) {
    index <- read_zarr_data_frame_index(store, name)
    df <- data.frame(row.names = seq_along(index))
  } else {
    df <- data.frame(columns)
  }

  if (isTRUE(include_index)) {
    index <- read_zarr_data_frame_index(store, name)

    # The default index name is not allowed as a column name so adjust it
    if (index_name == "_index") {
      rownames(df) <- index
    }

  }

  df
}

#' Read Zarr data frame index
#'
#' Read the index of a data frame from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return an object containing the index
#'
#' @noRd
read_zarr_data_frame_index <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  g <- pizzarr::zarr_open_group(store, path = name)

  attributes <- g$get_attrs()$to_list()
  index_name <- attributes$`_index`

  read_zarr_element(store, file.path(name, index_name))
}

#' Read multiple Zarr datatypes
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param column_order Vector of item names (in order)
#'
#' @return a named list
#'
#' @noRd
read_zarr_collection <- function(store, name, column_order) {
  columns <- list()
  for (col_name in column_order) {
    new_name <- paste0(name, "/", col_name)
    tryCatch({
      encoding <- read_zarr_encoding(store, new_name)
      columns[[col_name]] <- read_zarr_element(
        store = store,
        name = new_name,
        type = encoding$type,
        version = encoding$version
      )
    }, error = function(cond) {
      warning("Not reading file '", new_name, "' in collection")
    })
  }
  columns
}
