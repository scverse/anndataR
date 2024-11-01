#' Write H5AD element
#'
#' Write an element to an H5AD file
#'
#' @param value The value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' #' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param ... Additional arguments passed to writing functions
#'
#' @noRd
#'
#' @details
#' `write_zarr_element()` should always be used instead of any of the specific
#' writing functions as it contains additional boilerplate to make sure
#' elements are written correctly.
write_zarr_element <- function(value, store, name, compression = c("none", "gzip", "lzf"), stop_on_error = FALSE, ...) { # nolint
  compression <- match.arg(compression)

  # Delete the path if it already exists
  # TODO: https://github.com/keller-mark/pizzarr/issues/69

  # Sparse matrices
  write_fun <-
    if (inherits(value, "sparseMatrix")) { # Sparse matrices
      write_zarr_sparse_array
    } else if (is.factor(value)) { # Categoricals
      write_zarr_categorical
    } else if (is.list(value)) { # Lists and data frames
      if (is.data.frame(value)) {
        write_zarr_data_frame
      } else {
        write_zarr_mapping
      }
    } else if (is.character(value)) { # Character values
      if (length(value) == 1 && !is.matrix(value)) {
        write_zarr_string_scalar
      } else {
        write_zarr_string_array
      }
    } else if (is.numeric(value) || inherits(value, "denseMatrix")) { # Numeric values
      if (length(value) == 1 && !is.matrix(value)) {
        write_zarr_numeric_scalar
      } else if (is.integer(value) && any(is.na(value))) {
        write_zarr_nullable_integer
      } else {
        write_zarr_dense_array
      }
    } else if (is.logical(value)) { # Logical values
      if (any(is.na(value))) {
        write_zarr_nullable_boolean
      } else if (length(value) == 1) {
        # Single Booleans should be written as numeric scalars
        write_zarr_numeric_scalar
      } else {
        write_zarr_dense_array
      }
    } else { # Fail if unknown
      stop("Writing '", class(value), "' objects to H5AD files is not supported")
    }


  tryCatch(
    {
      write_fun(value = value, store = store, name = name, compression = compression, ...)
    },
    error = function(e) {
      message <- paste0(
        "Could not write element '", name, "' of type '", class(value), "':\n",
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

#' Write H5AD encoding
#'
#' Write H5AD encoding attributes to an element in an H5AD file
#'
#' @noRd
#'
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param encoding The encoding type to set
#' @param version The encoding version to set
write_zarr_encoding <- function(store, name, encoding, version) {
  g <- pizzarr::zarr_open(store, path = name)
  attrs <- g$get_attrs()

  attrs$set_item("encoding-type", encoding)
  attrs$set_item("encoding-version", version)
}

#' Write H5AD dense array
#'
#' Write a dense array to an H5AD file
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
#'
#' @noRd
write_zarr_dense_array <- function(value, store, name, compression, version = "0.2.0") {
  version <- match.arg(version)

  zarr_write_compressed(store, name, value, compression)

  # Write attributes
  write_zarr_encoding(store, name, "array", version)
}

#' Write H5AD sparse array
#'
#' Write a sparse array to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_sparse_array <- function(value, store, name, compression, version = "0.1.0") {
  version <- match.arg(version)

  # check types
  stopifnot(inherits(value, "sparseMatrix"))

  if (inherits(value, "RsparseMatrix")) {
    type <- "csr_matrix"
    indices_attr <- "j"
  } else if (inherits(value, "CsparseMatrix")) {
    type <- "csc_matrix"
    indices_attr <- "i"
  } else {
    stop(
      "Unsupported matrix format in ", name, ".",
      "Supported formats are RsparseMatrix and CsparseMatrix",
      "(and objects that inherit from those)."
    )
  }

  # Write sparse matrix
  g <- pizzarr::zarr_open_group(store, path = name)
  zarr_write_compressed(store, paste0(name, "/indices"), attr(value, indices_attr), compression)
  zarr_write_compressed(store, paste0(name, "/indptr"), value@p, compression)
  zarr_write_compressed(store, paste0(name, "/data"), value@x, compression)

  # Add encoding
  write_zarr_encoding(store, name, type, version)

  # Write shape attribute
  g$get_attrs()$set_item("shape", dim(value))
}

#' Write H5AD nullable boolean
#'
#' Write a nullable boolean to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_nullable_boolean <- function(value, store, name, compression, version = "0.1.0") {
  # write mask and values
  g <- pizzarr::zarr_open_group(store, path = name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE

  zarr_write_compressed(store, paste0(name, "/values"), value_no_na, compression)
  zarr_write_compressed(store, paste0(name, "/mask"), is.na(value), compression)

  # Write attributes
  write_zarr_encoding(store, name, "nullable-boolean", version)
}

#' Write H5AD nullable integer
#'
#' Write a nullable integer to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_nullable_integer <- function(value, store, name, compression, version = "0.1.0") {
  # write mask and values
  g <- pizzarr::zarr_open_group(store, path = name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- -1L

  zarr_write_compressed(store, paste0(name, "/values"), value_no_na, compression)
  zarr_write_compressed(store, paste0(name, "/mask"), is.na(value), compression)

  # Write attributes
  write_zarr_encoding(store, name, "nullable-integer", version)
}

#' Write H5AD string array
#'
#' Write a string array to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_string_array <- function(value, store, name, compression, version = "0.2.0") {

  if (!is.null(dim(value))) {
    dims <- dim(value)
  } else {
    dims <- length(value)
  }

  object_codec <- pizzarr::VLenUtf8Codec$new()
  data <- array(data = value, dim = dims)
  a <- pizzarr::zarr_create_array(data, store = store, path = name, dtype = "|O", object_codec = object_codec, shape = dims)

  write_zarr_encoding(store, name, "string-array", version)
}

#' Write H5AD categorical
#'
#' Write a categorical to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_categorical <- function(value, store, name, compression, version = "0.2.0") {
  g <- pizzarr::zarr_open_group(store, path = name)
  zarr_write_compressed(store, paste0(name, "/categories"), levels(value), compression)
  zarr_write_compressed(store, paste0(name, "/codes"), as.integer(value), compression)
  zarr_write_compressed(store, paste0(name, "/ordered"), is.ordered(value), compression)

  write_zarr_encoding(store, name, "categorical", version)
}

#' Write H5AD string scalar
#'
#' Write a string scalar to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_string_scalar <- function(value, store, name, compression, version = "0.2.0") {
  # Write scalar
  object_codec = pizzarr::VLenUtf8Codec$new()
  a <- pizzarr::zarr_create_array(value, store = store, path = name, dtype = "|O", object_codec = object_codec, shape = list())

  # Write attributes
  write_zarr_encoding(store, name, "string", version)
}

#' Write H5AD numeric scalar
#'
#' Write a numeric scalar to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_numeric_scalar <- function(value, store, name, compression, version = "0.2.0") {
  # Write scalar
  zarr_write_compressed(store, name, value, compression)

  # Write attributes
  write_zarr_encoding(store, name, "numeric-scalar", version)
}

#' Write H5AD mapping
#'
#' Write a mapping to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_mapping <- function(value, store, name, compression, version = "0.1.0") {
  g <- pizzarr::zarr_open_group(store, path = name)

  # Write mapping elements
  for (key in names(value)) {
    write_zarr_element(value[[key]], store, paste0(name, "/", key), compression)
  }

  write_zarr_encoding(store, name, "dict", version)
}

#' Write H5AD data frame
#'
#' Write a data frame to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param index The index to write. Can either be a vector of length equal to
#' the number of rows in `values` or a single character string giving the name
#' of a column in `values`. If `NULL` then `rownames(value)` is used.
#' @param version Encoding version of the element to write
write_zarr_data_frame <- function(value, store, name, compression, index = NULL,
                                  version = "0.2.0") {
  g <- pizzarr::zarr_open_group(store, path = name)
  write_zarr_encoding(store, name, "dataframe", version)

  if (is.null(index)) {
    index_name <- "_index"
    index_value <- rownames(value)
  } else if (length(index) == nrow(value)) {
    index_name <- "_index"
    index_value <- index
  } else if (length(index) == 1 && index %in% colnames(value)) {
    index_name <- index
    index_value <- value[[index_name]]
    value[[index_name]] <- NULL
  } else {
    stop(
      "index must be a vector with length `nrow(value)` or a single character",
      "string giving the name of a column in `value`"
    )
  }

  # Write index
  write_zarr_data_frame_index(index_value, store, name, compression, index_name)

  # Write data frame columns
  for (col in colnames(value)) {
    write_zarr_element(value[[col]], store, paste0(name, "/", col), compression)
  }

  # Write additional data frame attributes
  col_order <- colnames(value)
  col_order <- col_order[col_order != index_name]
  # If there are no columns other than the index we set column order to an
  # empty numeric vector
  if (length(col_order) == 0) {
    col_order <- numeric()
  }

  g$get_attrs()$set_item("column-order", col_order)
}

#' Write H5AD data frame index
#'
#' Write an for a data frame to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write. Must be a vector to the same length as the data
#' frame.
#' @param file Path to a H5AD file or an open H5AD handle
#' @param name Name of the element within the H5AD file containing the data
#' frame
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param index_name Name of the data frame column storing the index
write_zarr_data_frame_index <- function(value, store, name, compression, index_name) {
  if (!zarr_path_exists(store, name)) {
    warning("The data frame '", name, "' does not exist in store. Creating it.")
    g <- pizzarr::zarr_open_group(store, path = name)
    write_zarr_encoding(store, name, "dataframe", "0.2.0")
  }

  encoding <- read_zarr_encoding(store, name)
  if (encoding$type != "dataframe") {
    stop("'", name, "' in '", store, "' is not a data frame")
  }

  # Write index columns
  write_zarr_element(value, store, paste0(name, "/", index_name))

  # Write data frame index attribute
  g <- pizzarr::zarr_open_group(store, path = name)
  g$get_attrs()$set_item("_index", index_name)
}

#' Write empty H5AD
#'
#' Write a new empty H5AD file
#'
#' @noRd
#'
#' @param file Path to the Zarr store to write
#' @param obs Data frame with observations
#' @param var Data frame with variables
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version The H5AD version to write
write_empty_zarr <- function(store, obs, var, compression, version = "0.1.0") {
  pizzarr::zarr_open_group(store, path = "/")
  write_zarr_encoding(store, "/", "anndata", "0.1.0")

  # write_zarr_element(data.frame(row.names = obs_names), store, "/obs", compression)
  # write_zarr_element(data.frame(row.names = var_names), store, "/var", compression)
  write_zarr_element(obs, store, "/obs", compression)
  write_zarr_element(var, store, "/var", compression)

  pizzarr::zarr_open_group(store, path = "layers")
  write_zarr_encoding(store, "/layers", "dict", "0.1.0")

  pizzarr::zarr_open_group(store, path = "obsm")
  write_zarr_encoding(store, "/obsm", "dict", "0.1.0")

  pizzarr::zarr_open_group(store, path = "obsp")
  write_zarr_encoding(store, "/obsp", "dict", "0.1.0")

  pizzarr::zarr_open_group(store, path = "uns")
  write_zarr_encoding(store, "/uns", "dict", "0.1.0")

  pizzarr::zarr_open_group(store, path = "varm")
  write_zarr_encoding(store, "/varm", "dict", "0.1.0")

  pizzarr::zarr_open_group(store, path = "varp")
  write_zarr_encoding(store, "/varp", "dict", "0.1.0")
}

#' HDF5 path exists
#'
#' Check that a path in HDF5 exists
#'
#' @noRd
#'
#' @param store Path to a Zarr store
#' @param target_path The path within the file to test for
#'
#' @return Whether the `path` exists in `file`
zarr_path_exists <- function(store, target_path) {
  store <- pizzarr::zarr_open(store, path = target_path)
  result <- tryCatch({
    store$get_item(target_path)
    return(TRUE)
  }, error = function(cond) {
    if(pizzarr::is_key_error(cond)) {
      return(FALSE)
    }
    stop(cond)
  })
  return(result)
}

#' HDF5 write compressed
#'
#' Write HDF5 dataset with chosen compression (can be none)
#'
#' @noRd
#'
#' @param file Path to a HDF5 file
#' @param name Name of the element within the H5AD file containing the data
#' frame
#' @param value Value to write. Must be a vector to the same length as the data
#' frame.
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#'
#' @return Whether the `path` exists in `file`
zarr_write_compressed <- function(store, name, value, compression = c("none", "gzip", "lzf")) {
  compression <- match.arg(compression)
  if (!is.null(dim(value))) {
    dims <- dim(value)
  } else {
    dims <- length(value)
  }

  object_codec <- NA
  if(is.integer(value)) {
    dtype <- "<i8"
  } else if(is.numeric(value)) {
    dtype <- "<f8"
  } else if(is.logical(value)) {
    dtype <- "|b1"
  } else if(is.character(value)) {
    dtype <- "|O"
    object_codec <- pizzarr::VLenUtf8Codec$new()
  } else {
    stop("Unsupported data type for writing to Zarr: ", class(value))
  }
  data <- array(data = value, dim = dims)
  pizzarr::zarr_create_array(data, store = store, path = name, shape = dims, dtype = dtype, object_codec = object_codec) # TODO: compression
}
