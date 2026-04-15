#' Write Zarr element
#'
#' Write an element to a Zarr store
#'
#' @param value The value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"` or `"gzip"`. Defaults to `"none"`.
#' @param stop_on_error Whether to stop on error or generate a warning instead
#' @param ... Additional arguments passed to writing functions
#'
#' @noRd
#'
#' @details
#' `write_zarr_element()` should always be used instead of any of the specific
#' writing functions as it contains additional boilerplate to make sure
#' elements are written correctly.
write_zarr_element <- function(
  value,
  store,
  name,
  compression = c(
    "none",
    "gzip",
    "blosc",
    "zstd",
    "lzma",
    "bz2",
    "zlib",
    "lz4"
  ),
  stop_on_error = FALSE,
  ...
) {
  compression <- match.arg(compression)

  # Sparse matrices
  write_fun <-
    if (is.null(value)) {
      write_zarr_null
    } else if (inherits(value, "sparseMatrix")) {
      # Sparse matrices
      write_zarr_sparse_array
    } else if (is.factor(value)) {
      # Categoricals
      write_zarr_categorical
    } else if (is.list(value)) {
      # Lists and data frames
      if (is.data.frame(value)) {
        write_zarr_data_frame
      } else {
        write_zarr_mapping
      }
    } else if (is.character(value)) {
      # Character values
      if (length(value) == 1 && !is.matrix(value)) {
        write_zarr_string_scalar
      } else {
        write_zarr_string_array
      }
    } else if (is.numeric(value) || inherits(value, "denseMatrix")) {
      # Numeric values
      if (length(value) == 1 && !is.matrix(value)) {
        write_zarr_numeric_scalar
      } else if (is.integer(value) && any(is.na(value))) {
        write_zarr_nullable_integer
      } else {
        write_zarr_dense_array
      }
    } else if (is.logical(value)) {
      # Logical values
      if (any(is.na(value))) {
        write_zarr_nullable_boolean
      } else if (length(value) == 1) {
        # Single Booleans should be written as numeric scalars
        write_zarr_numeric_scalar
      } else {
        write_zarr_dense_array
      }
    } else {
      # Fail if unknown
      cli_abort(c(
        "Writing {.cls {class(value)}} objects to Zarr is not supported",
        "i" = "Attempting to write to {.path {name}} in {.file {store}}"
      ))
    }

  # Delete the path if it already exists
  if (zarr_path_exists(store, name)) {
    unlink(file.path(store, name), recursive = TRUE)
  }

  tryCatch(
    {
      write_fun(
        value = value,
        store = store,
        name = name,
        compression = compression,
        ...
      )
    },
    error = function(e) {
      message <- paste0(
        "Could not write element '",
        name,
        "' of type '",
        class(value),
        "':\n",
        conditionMessage(e)
      )
      if (stop_on_error) {
        cli_abort(message)
      } else {
        cli_warn(message)
        NULL
      }
    }
  )
}

#' Write Zarr encoding
#'
#' Write Zarr encoding attributes to an element in a Zarr store
#'
#' @noRd
#'
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param encoding The encoding type to set
#' @param version The encoding version to set
write_zarr_encoding <- function(store, name, encoding, version) {
  Rarr::write_zarr_attributes(
    file.path(store, name),
    new.zattrs = list(`encoding-type` = encoding, `encoding-version` = version)
  )
}

#' Write Zarr null
#'
#' Write a null dataset to an Zarr file
#'
#' @param value Value to write, not used
#' @param store An open Zarr handle
#' @param name Name of the element within the Zarr store
#' @param compression Not used as there is no value
#' @param version Encoding version of the element to write
#'
#' @noRd
write_zarr_null <- function(
  value,
  store,
  name,
  compression,
  version = "0.1.0"
) {
  if (isFALSE(getOption("anndataR.write_null", "TRUE"))) {
    return(invisible(NULL))
  }

  Rarr::create_empty_zarr_array(
    file.path(store, name),
    dim = 0,
    chunk_dim = 0,
    data_type = "logical",
    zarr_version = 2L
  )

  write_zarr_encoding(store, name, "null", version)
}

#' Write Zarr dense array
#'
#' Write a dense array to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_dense_array <- function(
  value,
  store,
  name,
  compression,
  version = "0.2.0"
) {
  version <- match.arg(version)

  # matrices of type 'dgeMatrix' can simply be converted to a matrix
  if (inherits(value, "denseMatrix")) {
    value <- as.matrix(value)
  }

  zarr_write_compressed(
    store,
    name,
    value,
    compression
  )

  # Write attributes
  write_zarr_encoding(store, name, "array", version)
}

#' Write Zarr sparse array
#'
#' Write a sparse array to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_sparse_array <- function(
  value,
  store,
  name,
  compression,
  version = "0.1.0"
) {
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
    cli_abort(c(
      "Unsupported matrix format in {.path {name}}",
      "i" = "Supported matrices inherit from {.cls RsparseMatrix} or {.cls CsparseMatrix}"
    ))
  }

  # Write sparse matrix
  create_zarr_group(store, name)
  zarr_write_compressed(
    store,
    paste0(name, "/indices"),
    attr(value, indices_attr),
    compression
  )
  zarr_write_compressed(
    store,
    paste0(name, "/indptr"),
    value@p,
    compression
  )
  zarr_write_compressed(
    store,
    paste0(name, "/data"),
    value@x,
    compression
  )

  # Add encoding
  write_zarr_encoding(store, name, type, version)

  # Write shape attribute
  Rarr::write_zarr_attributes(file.path(store, name), list(shape = dim(value)))
}

#' Write Zarr nullable boolean
#'
#' Write a nullable boolean to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_nullable_boolean <- function(
  value,
  store,
  name,
  compression,
  version = "0.1.0"
) {
  # write mask and values
  create_zarr_group(store, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE

  zarr_write_compressed(
    store,
    paste0(name, "/values"),
    value_no_na,
    compression
  )
  zarr_write_compressed(
    store,
    paste0(name, "/mask"),
    is.na(value),
    compression
  )

  # Write attributes
  write_zarr_encoding(store, name, "nullable-boolean", version)
}

#' Write Zarr nullable integer
#'
#' Write a nullable integer to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_nullable_integer <- function(
  value,
  store,
  name,
  compression,
  version = "0.1.0"
) {
  # write mask and values
  create_zarr_group(store, name)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- -1L

  zarr_write_compressed(
    store,
    paste0(name, "/values"),
    value_no_na,
    compression
  )
  zarr_write_compressed(
    store,
    paste0(name, "/mask"),
    is.na(value),
    compression
  )

  # Write attributes
  write_zarr_encoding(store, name, "nullable-integer", version)
}

#' Write Zarr string array
#'
#' Write a string array to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_string_array <- function(
  value,
  store,
  name,
  compression,
  version = "0.2.0"
) {
  dims <- dim(value) %||% length(value)

  # replace NA to "NA" (as in rhdf5:::.h5postProcessDataset)
  # to read as "NA" -> NA later after Rarr:read_zarr_array
  value[is.na(value)] <- "NA"

  if (any(dims == 0)) {
    Rarr::create_empty_zarr_array(
      file.path(store, name),
      dim = dims,
      chunk_dim = dims,
      data_type = "<U",
      nchar = 1,
      compressor = .get_compressor(compression),
      zarr_version = 2L
    )
  } else {
    data <- array(data = value, dim = dims)
    Rarr::write_zarr_array(
      data,
      zarr_array_path = file.path(store, name),
      chunk_dim = dims,
      order = if (length(dims) > 1) "C" else "F",
      # TODO: string arrays require vlen-utf8 filter support
      # see https://github.com/Huber-group-EMBL/Rarr/issues/98
      data_type = "<U",
      nchar = max(nchar(value)),
      compressor = .get_compressor(compression),
      zarr_version = 2L
    )
  }

  write_zarr_encoding(store, name, "string-array", version)
}

#' Write Zarr categorical
#'
#' Write a categorical to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_categorical <- function(
  value,
  store,
  name,
  compression,
  version = "0.2.0"
) {
  create_zarr_group(store, name)

  categories <- levels(value)

  # Use zero-indexed values
  codes <- as.integer(value) - 1L

  # Set missing values to -1
  codes[is.na(codes)] <- -1L

  # write values to file
  write_zarr_string_array(
    categories,
    store,
    paste0(name, "/categories"),
    compression
  )
  write_zarr_dense_array(codes, store, paste0(name, "/codes"), compression)

  # Write encoding
  write_zarr_encoding(
    store = store,
    name = name,
    encoding = "categorical",
    version = version
  )

  # Write ordered attribute
  Rarr::write_zarr_attributes(
    file.path(store, name),
    new.zattrs = list("ordered" = is.ordered(value))
  )
}

#' Write Zarr string scalar
#'
#' Write a string scalar to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_string_scalar <- function(
  value,
  store,
  name,
  compression,
  version = "0.2.0"
) {
  # Write scalar
  value <- array(data = value, dim = 1)
  Rarr::write_zarr_array(
    value,
    zarr_array_path = file.path(store, name),
    chunk_dim = 1,
    compressor = .get_compressor(compression),
    zarr_version = 2L
  )

  # Write attributes
  write_zarr_encoding(store, name, "string", version)
}

#' Write Zarr numeric scalar
#'
#' Write a numeric scalar to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_numeric_scalar <- function(
  value,
  store,
  name,
  compression,
  version = "0.2.0"
) {
  # Write scalar
  zarr_write_compressed(store, name, value, compression)

  # Write attributes
  write_zarr_encoding(store, name, "numeric-scalar", version)
}

#' Write Zarr mapping
#'
#' Write a mapping to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
write_zarr_mapping <- function(
  value,
  store,
  name,
  compression,
  version = "0.1.0"
) {
  create_zarr_group(store, name)

  # Write mapping elements
  for (key in names(value)) {
    write_zarr_element(
      value[[key]],
      store,
      paste0(name, "/", key),
      compression
    )
  }

  write_zarr_encoding(store, name, "dict", version)
}

#' Write Zarr data frame
#'
#' Write a data frame to a Zarr store
#'
#' @noRd
#'
#' @param value Value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param index The index to write. Can either be a vector of length equal to
#' the number of rows in `values` or a single character string giving the name
#' of a column in `values`. If `NULL` then `rownames(value)` is used.
#' @param version Encoding version of the element to write
write_zarr_data_frame <- function(
  value,
  store,
  name,
  compression,
  index = NULL,
  version = "0.2.0"
) {
  create_zarr_group(store, name)
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
    cli_abort(paste(
      "{.arg index} must be a vector with length {.code nrow(value)} or",
      "a single character vector giving the name of a column in {.arg value}"
    ))
  }
  if (is.null(index_value)) {
    index_value <- seq_len(nrow(value)) - 1L
  }

  # Write data frame columns
  for (col in colnames(value)) {
    write_zarr_element(
      value[[col]],
      store,
      paste0(name, "/", col),
      compression
    )
  }

  # Write index
  write_zarr_element(
    index_value,
    store,
    paste0(name, "/", index_name),
    compression
  )

  # Write additional data frame attributes
  Rarr::write_zarr_attributes(
    zarr_path = file.path(store, name),
    new.zattrs = list("_index" = index_name)
  )

  # Write additional data frame attributes
  col_order <- colnames(value)
  col_order <- col_order[col_order != index_name]
  # If there are no columns other than the index we set column order to an
  # empty numeric vector
  if (length(col_order) == 0) {
    col_order <- numeric()
  } else {
    col_order <- array(col_order)
  }

  Rarr::write_zarr_attributes(
    zarr_path = file.path(store, name),
    new.zattrs = list(`column-order` = col_order)
  )
}

#' Write empty Zarr
#'
#' Write a new empty Zarr store
#'
#' @noRd
#'
#' @param store Path to the Zarr store to write
#' @param obs Data frame with observations
#' @param var Data frame with variables
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"` or `"gzip"`. Defaults to `"none"`.
#' @param version The anndata on-disk format version to write
write_empty_zarr <- function(
  store,
  obs,
  var,
  compression,
  version = "0.1.0"
) {
  create_zarr(store = store)
  write_zarr_encoding(store, "/", "anndata", "0.1.0")

  write_zarr_element(obs[, integer(0)], store, "/obs", compression)
  write_zarr_element(var[, integer(0)], store, "/var", compression)

  create_zarr_group(store, "layers")
  write_zarr_encoding(store, "/layers", "dict", "0.1.0")

  create_zarr_group(store, "obsm")
  write_zarr_encoding(store, "/obsm", "dict", "0.1.0")

  create_zarr_group(store, "obsp")
  write_zarr_encoding(store, "/obsp", "dict", "0.1.0")

  create_zarr_group(store, "uns")
  write_zarr_encoding(store, "/uns", "dict", "0.1.0")

  create_zarr_group(store, "varm")
  write_zarr_encoding(store, "/varm", "dict", "0.1.0")

  create_zarr_group(store, "varp")
  write_zarr_encoding(store, "/varp", "dict", "0.1.0")
}

#' Zarr write compressed
#'
#' Write Zarr dataset with chosen compression (can be none)
#'
#' @return  Returns (invisibly) `TRUE` if the array is successfully written.
#'
#' @noRd
#'
#' @param store Path to a Zarr store
#' @param name Name of the element within the Zarr store containing the data
#' frame
#' @param value Value to write. Must be a vector to the same length as the data
#' frame.
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"` or `"gzip"`. Defaults to `"none"`.
zarr_write_compressed <- function(
  store,
  name,
  value,
  compression = c("none", "gzip", "blosc", "zstd", "lzma", "bz2", "zlib", "lz4")
) {
  dims <- dim(value) %||% length(value)
  data <- array(data = value, dim = dims)
  Rarr::write_zarr_array(
    data,
    zarr_array_path = file.path(store, name),
    chunk_dim = dims,
    order = if (length(dims) > 1) "C" else "F",
    compressor = .get_compressor(compression),
    zarr_version = 2L
  )
}

#' Get Zarr compressor
#'
#' Convert a compression name to the corresponding Rarr compressor object.
#'
#' @param compression The compression algorithm name. One of `"none"`,
#'   `"gzip"`, `"blosc"`, `"zstd"`, `"lzma"`, `"bz2"`, `"zlib"`, `"lz4"`.
#'
#' @return A Rarr compressor object, or `NULL` for no compression.
#'
#' @noRd
.get_compressor <- function(compression) {
  switch(
    compression,
    "none" = NULL,
    "zstd" = Rarr::use_zstd(),
    "blosc" = Rarr::use_blosc(),
    "gzip" = Rarr::use_gzip(),
    "lzma" = Rarr::use_lzma(),
    "bz2" = Rarr::use_bz2(),
    "zlib" = Rarr::use_zlib(),
    "lz4" = Rarr::use_lz4()
  )
}
