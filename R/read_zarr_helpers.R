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
read_zarr_encoding <- function(store, name) {
  tryCatch(
    {
      attrs <- Rarr::read_zarr_attributes(file.path(store, name))
      list(
        type = attrs[["encoding-type"]],
        version = attrs[["encoding-version"]]
      )
    },
    error = function(e) {
      cli_abort(
        "Encoding attributes not found for element {.val {name}} in {.path {store}}"
      )
    }
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
read_zarr_element <- function(
  store,
  name,
  type = NULL,
  version = NULL,
  stop_on_error = FALSE,
  ...
) {
  if (!zarr_path_exists(store, name)) {
    return(NULL)
  }

  if (is.null(type)) {
    encoding_list <- read_zarr_encoding(store, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(
    type,
    "null" = read_zarr_null,
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
    cli_abort(
      "No function for reading Zarr encoding {.cls {type}} for element {.val {name}}"
    )
  )

  tryCatch(
    {
      read_fun(store = store, name = name, version = version, ...)
    },
    error = function(e) {
      msg <- cli::cli_fmt(cli::cli_bullets(c(
        paste0(
          "Error reading element {.field {name}} of type {.cls {type}}"
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

#' Read Zarr null
#'
#' Read a null value from an Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return `NULL`
#' @noRd
read_zarr_null <- function(store, name, version = "0.1.0") {
  version <- match.arg(version)

  NULL
}

#' Read Zarr dense array
#'
#' Read a dense array from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A matrix or a vector if 1D
#'
#' @noRd
read_zarr_dense_array <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  data <- Rarr::read_zarr_array(file.path(store, name))

  data
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
#' @return A sparse matrix/DelayedArray???, or a vector if 1D
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
read_zarr_sparse_array <- function(
  store,
  name,
  version = "0.1.0",
  type = c("csr_matrix", "csc_matrix")
) {
  version <- match.arg(version)
  type <- match.arg(type)

  attrs <- Rarr::read_zarr_attributes(file.path(store, name))

  construct_sparse_matrix(
    data = as.vector(Rarr::read_zarr_array(file.path(store, name, "data"))),
    indices = as.vector(Rarr::read_zarr_array(file.path(
      store,
      name,
      "indices"
    ))),
    indptr = as.vector(Rarr::read_zarr_array(file.path(store, name, "indptr"))),
    shape = as.vector(unlist(attrs$shape, use.names = FALSE)),
    type = type
  )
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
#' @return A named list of 1D arrays
#'
#' @noRd
read_zarr_rec_array <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)
  Rarr::read_zarr_array(file.path(store, name)) |>
    lapply(as.vector)
}

#' Read Zarr nullable boolean
#'
#' Read a nullable boolean from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A boolean vector
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
#' @return An integer vector
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
#' @return A nullable vector
#'
#' @noRd
read_zarr_nullable <- function(store, name, version = "0.1.0") {
  version <- match.arg(version)

  mask <- Rarr::read_zarr_array(file.path(store, paste0(name, "/mask")))
  values <- Rarr::read_zarr_array(file.path(store, paste0(name, "/values")))

  # Get values and set missing
  element <- values
  element[mask] <- NA

  element
}

#' Read Zarr string array
#'
#' Read a string array from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A character vector/matrix
#'
#' @noRd
read_zarr_string_array <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  data <- Rarr::read_zarr_array(file.path(store, name))

  # convert "NA" to NA (as in rhdf5:::.h5postProcessDataset)
  data[data == "NA"] <- NA

  data
}

#' Read Zarr categorical
#'
#' Read a categorical from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A factor
#'
#' @noRd
read_zarr_categorical <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  codes <- Rarr::read_zarr_array(file.path(store, paste0(name, "/codes")))
  categories <- Rarr::read_zarr_array(file.path(
    store,
    paste0(name, "/categories")
  ))

  # Get codes and convert to 1-based indexing
  codes <- codes + 1L

  # Set missing values
  codes[codes == 0L] <- NA_integer_

  levels <- categories

  attributes <- Rarr::read_zarr_attributes(file.path(store, name))
  ordered <- attributes[["ordered"]]

  factor(levels[codes], levels = levels, ordered = ordered)
}

#' Read Zarr string scalar
#'
#' Read a string scalar from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A character vector of length 1
#'
#' @noRd
read_zarr_string_scalar <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)
  as.character(Rarr::read_zarr_array(file.path(store, name)))
}

#' Read Zarr numeric scalar
#'
#' Read a numeric scalar from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A numeric vector of length 1
#'
#' @noRd
read_zarr_numeric_scalar <- function(store, name, version = "0.2.0") {
  version <- match.arg(version)

  value <- Rarr::read_zarr_array(file.path(store, name))

  # convert array to vector
  value <- as.vector(value)

  value
}

#' Read Zarr mapping
#'
#' Read a mapping from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A named list
#'
#' @noRd
read_zarr_mapping <- function(store, name, version = "0.1.0") {
  version <- match.arg(version)
  items <- read_zarr_mapping_keys(store, name, version)
  read_zarr_collection(store, name, items)
}

#' Read Zarr data frame
#'
#' Read a data frame from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A data.frame
#'
#' @noRd
read_zarr_data_frame <- function(
  store,
  name,
  version = "0.2.0"
) {
  version <- match.arg(version)

  dim_keys <- read_zarr_data_frame_keys(store, name, version)
  data <- read_zarr_collection(store, name, dim_keys$cols)

  as.data.frame(
    row.names = dim_keys$rows,
    data,
    check.names = FALSE,
    fix.empty.names = FALSE
  )
}

#' Read multiple Zarr datatypes
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param item_names Vector of item names (in order)
#'
#' @return A named list
#'
#' @noRd
read_zarr_collection <- function(store, name, item_names) {
  items <- lapply(
    item_names,
    function(item_name) {
      new_name <- paste0(name, "/", item_name)
      encoding <- read_zarr_encoding(store, new_name)
      read_zarr_element(
        store = store,
        name = new_name,
        type = encoding$type,
        version = encoding$version
      )
    }
  )
  names(items) <- item_names
  items
}

#' Read Zarr element keys
#'
#' Read the keys of an element from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param type The encoding type of the element to read
#' @param version The encoding version of the element to read
#' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param ... Extra arguments passed to individual reading functions
#'
#' @return A character vector of keys
#'
#' @noRd
read_zarr_element_keys <- function(
  store,
  name,
  type = NULL,
  version = NULL,
  stop_on_error = FALSE,
  ...
) {
  if (!zarr_path_exists(store, name)) {
    return(NULL)
  }

  if (is.null(type)) {
    encoding_list <- read_zarr_encoding(store, name)
    type <- encoding_list$type
    version <- encoding_list$version
  }

  read_fun <- switch(
    type,
    "dataframe" = read_zarr_data_frame_keys,
    "dict" = read_zarr_mapping_keys,
    cli_abort(
      "No function for reading keys for Zarr encoding {.cls {type}} for element {.val {name}}"
    )
  )

  tryCatch(
    {
      read_fun(store = store, name = name, version = version, ...)
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

#' Read Zarr mapping keys
#'
#' Read keys for a mapping (dict) from a Zarr store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#'
#' @return A character vector of item names
#'
#' @noRd
read_zarr_mapping_keys <- function(store, name, version = "0.1.0") {
  version <- match.arg(version)

  items <- list.dirs(
    path = file.path(store, name),
    recursive = FALSE,
    full.names = FALSE
  )
  items[!items %in% ZARR_METADATA_FILES]
}

#' Read Zarr data frame keys
#'
#' Read the row names (index) and/or column names of a data frame from a Zarr
#' store
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param version Encoding version of the element to read
#' @param dim Dimension to read keys for: `"both"`, `"rows"`, or `"cols"`
#'
#' @return A character vector if `dim` is `"rows"` or `"cols"`, or a list with
#'   elements `"rows"` and `"cols"` if `dim` is `"both"`
#'
#' @noRd
read_zarr_data_frame_keys <- function(
  store,
  name,
  version = "0.2.0",
  dim = c("both", "rows", "cols")
) {
  version <- match.arg(version)
  dim <- match.arg(dim)

  attrs <- Rarr::read_zarr_attributes(file.path(store, name))
  index_name <- attrs[["_index"]]
  column_order <- attrs[["column-order"]]

  if (dim == "both") {
    list(
      rows = as.vector(read_zarr_element(store, file.path(name, index_name))),
      cols = as.character(column_order)
    )
  } else if (dim == "rows") {
    as.vector(read_zarr_element(store, file.path(name, index_name)))
  } else if (dim == "cols") {
    as.character(column_order)
  }
}
