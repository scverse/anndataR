#' Write Zarr element
#'
#' Write an element to a Zarr store
#'
#' @param value The value to write
#' @param store A Zarr store instance
#' @param name Name of the element within the Zarr store
#' @param compression The compression to use when writing the element.
#'   One of `"none"`, `"gzip"`, `"blosc"`, `"zstd"`, `"lzma"`, `"bz2"`,
#'   `"zlib"`, `"lz4"`.
#' @param zarr_format The format to write the element in, either 2 or 3 for the
#'   Zarr v2 or v3 formats. Always passed explicitly by the caller so that all
#'   elements of a store end up in the same format.
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
  zarr_format,
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
      } else if (anyNA(value) && length(dim(value)) <= 1) {
        write_zarr_nullable_string
      } else {
        write_zarr_string_array
      }
    } else if (is.numeric(value) || inherits(value, "denseMatrix")) {
      # Numeric values
      if (length(value) == 1 && !is.matrix(value)) {
        write_zarr_numeric_scalar
      } else if (is.integer(value) && anyNA(value)) {
        write_zarr_nullable_integer
      } else {
        write_zarr_dense_array
      }
    } else if (is.logical(value)) {
      # Logical values
      if (anyNA(value)) {
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
        zarr_format = zarr_format,
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
#' @inheritParams write_zarr_element
#' @param encoding The encoding type to set
#' @param version The encoding version to set
#'
#' @noRd
write_zarr_encoding <- function(store, name, encoding, version, zarr_format) {
  Rarr::write_zarr_attributes(
    file.path(store, name),
    new.zattrs = list(`encoding-type` = encoding, `encoding-version` = version),
    zarr_version = zarr_format
  )
}

#' Write Zarr null
#'
#' Write a null dataset to an Zarr file
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_null <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.1.0"
) {
  if (isFALSE(getOption("anndataR.write_null", "TRUE"))) {
    return(invisible(NULL))
  }
  # if dims is zero, fix chunk dim to 1, but raises warnings
  # https://github.com/Huber-group-EMBL/Rarr/issues/89
  suppressWarnings({
    Rarr::create_empty_zarr_array(
      file.path(store, name),
      dim = 0,
      chunk_dim = 1,
      data_type = "logical",
      zarr_version = zarr_format
    )
  })

  write_zarr_encoding(store, name, "null", version, zarr_format)
}

#' Write Zarr dense array
#'
#' Write a dense array to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_dense_array <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
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
    compression,
    zarr_format = zarr_format
  )

  # Write attributes
  write_zarr_encoding(store, name, "array", version, zarr_format)
}

#' Write Zarr sparse array
#'
#' Write a sparse array to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_sparse_array <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
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
  create_zarr_group(store, name, zarr_format)
  zarr_write_compressed(
    store,
    paste0(name, "/indices"),
    attr(value, indices_attr),
    compression,
    zarr_format = zarr_format
  )
  zarr_write_compressed(
    store,
    paste0(name, "/indptr"),
    value@p,
    compression,
    zarr_format = zarr_format
  )
  zarr_write_compressed(
    store,
    paste0(name, "/data"),
    value@x,
    compression,
    zarr_format = zarr_format
  )

  # Add encoding
  write_zarr_encoding(store, name, type, version, zarr_format)

  # Write shape attribute
  Rarr::write_zarr_attributes(
    file.path(store, name),
    list(shape = dim(value)),
    zarr_version = zarr_format
  )
}

#' Write Zarr nullable boolean
#'
#' Write a nullable boolean to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_nullable_boolean <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.1.0"
) {
  # write mask and values
  create_zarr_group(store, name, zarr_format)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE

  zarr_write_compressed(
    store,
    paste0(name, "/values"),
    value_no_na,
    compression,
    zarr_format = zarr_format
  )
  zarr_write_compressed(
    store,
    paste0(name, "/mask"),
    is.na(value),
    compression,
    zarr_format = zarr_format
  )

  # Write attributes
  write_zarr_encoding(store, name, "nullable-boolean", version, zarr_format)
}

#' Write Zarr nullable integer
#'
#' Write a nullable integer to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_nullable_integer <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.1.0"
) {
  # write mask and values
  create_zarr_group(store, name, zarr_format)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- -1L

  zarr_write_compressed(
    store,
    paste0(name, "/values"),
    value_no_na,
    compression,
    zarr_format = zarr_format
  )
  zarr_write_compressed(
    store,
    paste0(name, "/mask"),
    is.na(value),
    compression,
    zarr_format = zarr_format
  )

  # Write attributes
  write_zarr_encoding(store, name, "nullable-integer", version, zarr_format)
}

#' Patch Zarr array metadata to declare VLen-UTF8 strings
#'
#' TEMPORARY WORKAROUND. Rarr has no interface for writing VLen-UTF8 strings, so
#' the array is written with a placeholder data type and its metadata is then
#' rewritten by hand: a `vlen-utf8` filter for Zarr v2, and the `bytes` codec
#' swapped for `vlen-utf8` plus a `string` data type for Zarr v3.
#'
#' Editing metadata behind Rarr's back like this is fragile, and it means the
#' array on disk briefly disagrees with what Rarr thinks it wrote.
#'
#' TODO: Remove this function and pass the string type to Rarr directly once
#' https://github.com/Huber-group-EMBL/Rarr/issues/111 is resolved.
#'
#' @param store The location of the Zarr store
#' @param name Name of the array within the store
#' @param zarr_format The format the array was written in
#'
#' @noRd
patch_zarr_vlen_utf8 <- function(store, name, zarr_format) {
  metadata_path <- file.path(
    store,
    name,
    if (zarr_format == 2L) ".zarray" else "zarr.json"
  )
  check_requires("write_zarr", "jsonlite", where = "CRAN")
  metadata <- jsonlite::read_json(metadata_path)

  if (zarr_format == 2L) {
    metadata$filters <- c(metadata$filters, list(list(id = "vlen-utf8")))
  } else {
    metadata$data_type <- "string"
    # There should be only one array-bytes codec
    metadata$codecs <- lapply(
      metadata$codecs,
      function(codec) {
        if (codec$name == "bytes") {
          list(name = "vlen-utf8")
        } else {
          codec
        }
      }
    )
  }

  jsonlite::write_json(
    metadata,
    metadata_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
}

#' Write Zarr nullable string
#'
#' Write a nullable string to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_nullable_string <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.1.0"
) {
  # write mask and values
  create_zarr_group(store, name, zarr_format)
  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- ""

  write_zarr_string_array(
    value_no_na,
    store,
    paste0(name, "/values"),
    compression,
    zarr_format = zarr_format
  )
  zarr_write_compressed(
    store,
    paste0(name, "/mask"),
    is.na(value),
    compression,
    zarr_format = zarr_format
  )

  # Write attributes
  write_zarr_encoding(
    store,
    name,
    "nullable-string-array",
    version,
    zarr_format
  )

  Rarr::write_zarr_attributes(
    file.path(store, name),
    new.zattrs = list(`na-value` = "NA"),
    zarr_version = zarr_format
  )
}

#' Write Zarr string array
#'
#' Write a string array to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_string_array <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.2.0"
) {
  dims <- dim(value) %||% length(value)

  # if dims is zero, fix chunk dim to 1, but raises warnings
  # https://github.com/Huber-group-EMBL/Rarr/issues/89
  chunk_dims <- ifelse(dims == 0, 1, dims)

  # replace NA to "NA" (as in rhdf5:::.h5postProcessDataset)
  # to read as "NA" -> NA later after Rarr:read_zarr_array
  value[is.na(value)] <- "NA"

  # suppress chunk dim warnings
  suppressWarnings({
    Rarr::create_empty_zarr_array(
      file.path(store, name),
      dim = dims,
      chunk_dim = chunk_dims,
      order = if (length(dims) > 0) "C" else "F",
      # placeholder type, rewritten to a string type just below
      data_type = "|O",
      compressor = .get_compressor(compression),
      zarr_version = zarr_format
    )
  })
  # TODO: Remove this call, and write the array as strings directly, once
  # https://github.com/Huber-group-EMBL/Rarr/issues/111 is resolved
  patch_zarr_vlen_utf8(store, name, zarr_format)

  if (all(dims != 0)) {
    data <- array(data = value, dim = dims)
    Rarr::update_zarr_array(
      data,
      zarr_array_path = file.path(store, name),
      index = lapply(dims, seq_len)
    )
  }

  write_zarr_encoding(store, name, "string-array", version, zarr_format)
}

#' Write Zarr categorical
#'
#' Write a categorical to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_categorical <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.2.0"
) {
  create_zarr_group(store, name, zarr_format)

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
    compression,
    zarr_format = zarr_format
  )
  write_zarr_dense_array(
    codes,
    store,
    paste0(name, "/codes"),
    compression,
    zarr_format
  )

  # Write encoding
  write_zarr_encoding(
    store = store,
    name = name,
    encoding = "categorical",
    version = version,
    zarr_format = zarr_format
  )

  # Write ordered attribute
  Rarr::write_zarr_attributes(
    file.path(store, name),
    new.zattrs = list("ordered" = is.ordered(value)),
    zarr_version = zarr_format
  )
}

#' Write Zarr string scalar
#'
#' Write a string scalar to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_string_scalar <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.2.0"
) {
  # Write scalar
  value <- array(data = value, dim = 1)
  Rarr::create_empty_zarr_array(
    zarr_array_path = file.path(store, name),
    # placeholder type, rewritten to a string type just below
    data_type = "|O",
    dim = 1,
    chunk_dim = 1,
    compressor = .get_compressor(compression),
    zarr_version = zarr_format
  )
  # TODO: Remove this call, and write the scalar as a string directly, once
  # https://github.com/Huber-group-EMBL/Rarr/issues/111 is resolved
  patch_zarr_vlen_utf8(store, name, zarr_format)
  Rarr::update_zarr_array(
    value,
    zarr_array_path = file.path(store, name),
    index = list(1)
  )

  # Write attributes
  write_zarr_encoding(store, name, "string", version, zarr_format)
}

#' Write Zarr numeric scalar
#'
#' Write a numeric scalar to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_numeric_scalar <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.2.0"
) {
  # Write scalar
  zarr_write_compressed(store, name, value, compression, zarr_format)

  # Write attributes
  write_zarr_encoding(store, name, "numeric-scalar", version, zarr_format)
}

#' Write Zarr mapping
#'
#' Write a mapping to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
write_zarr_mapping <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  version = "0.1.0"
) {
  create_zarr_group(store, name, zarr_format)

  # Write mapping elements
  for (key in names(value)) {
    write_zarr_element(
      value[[key]],
      store,
      paste0(name, "/", key),
      compression,
      zarr_format = zarr_format
    )
  }

  write_zarr_encoding(store, name, "dict", version, zarr_format)
}

#' Write Zarr data frame
#'
#' Write a data frame to a Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#' @param index The index to write. Can either be a vector of length equal to
#' the number of rows in `values` or a single character string giving the name
#' of a column in `values`. If `NULL` then `rownames(value)` is used.
#'
#' @noRd
write_zarr_data_frame <- function(
  value,
  store,
  name,
  compression,
  zarr_format,
  index = NULL,
  version = "0.2.0"
) {
  create_zarr_group(store, name, zarr_format)
  write_zarr_encoding(store, name, "dataframe", version, zarr_format)

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
      compression,
      zarr_format
    )
  }

  # Write index
  write_zarr_element(
    index_value,
    store,
    paste0(name, "/", index_name),
    compression,
    zarr_format
  )

  # Write additional data frame attributes
  Rarr::write_zarr_attributes(
    zarr_path = file.path(store, name),
    new.zattrs = list("_index" = index_name),
    zarr_version = zarr_format
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
    new.zattrs = list(`column-order` = col_order),
    zarr_version = zarr_format
  )
}

#' Write empty Zarr
#'
#' Write a new empty Zarr store
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#' @param obs Data frame with observations
#' @param var Data frame with variables
#'
#' @noRd
write_empty_zarr <- function(
  store,
  obs,
  var,
  compression,
  zarr_format,
  version = "0.1.0"
) {
  create_zarr(store = store, zarr_format)
  write_zarr_encoding(store, "/", "anndata", "0.1.0", zarr_format)

  write_zarr_element(
    obs[, integer(0)],
    store,
    "/obs",
    compression,
    zarr_format
  )
  write_zarr_element(
    var[, integer(0)],
    store,
    "/var",
    compression,
    zarr_format
  )

  create_zarr_group(store, "layers", zarr_format)
  write_zarr_encoding(store, "/layers", "dict", "0.1.0", zarr_format)

  create_zarr_group(store, "obsm", zarr_format)
  write_zarr_encoding(store, "/obsm", "dict", "0.1.0", zarr_format)

  create_zarr_group(store, "obsp", zarr_format)
  write_zarr_encoding(store, "/obsp", "dict", "0.1.0", zarr_format)

  create_zarr_group(store, "uns", zarr_format)
  write_zarr_encoding(store, "/uns", "dict", "0.1.0", zarr_format)

  create_zarr_group(store, "varm", zarr_format)
  write_zarr_encoding(store, "/varm", "dict", "0.1.0", zarr_format)

  create_zarr_group(store, "varp", zarr_format)
  write_zarr_encoding(store, "/varp", "dict", "0.1.0", zarr_format)
}

#' Zarr write compressed
#'
#' Write Zarr dataset with chosen compression (can be none)
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
#'
#' @return  Returns (invisibly) `TRUE` if the array is successfully written.
zarr_write_compressed <- function(
  store,
  name,
  value,
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
  zarr_format
) {
  dims <- dim(value) %||% length(value)
  data <- array(data = value, dim = dims)
  Rarr::write_zarr_array(
    data,
    zarr_array_path = file.path(store, name),
    chunk_dim = dims,
    order = if (length(dims) > 0) "C" else "F",
    compressor = .get_compressor(compression),
    zarr_version = zarr_format
  )
}

#' Get Zarr compressor
#'
#' Convert a compression name to the corresponding Rarr compressor object.
#'
#' @inheritParams write_zarr_element
#' @inheritParams write_zarr_encoding version
#'
#' @noRd
#'
#' @return A Rarr compressor object, or `NULL` for no compression.
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
