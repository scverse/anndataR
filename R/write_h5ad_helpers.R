# should take a look at
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html
# again

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
#' `write_h5ad_element()` should always be used instead of any of the specific
#' writing functions as it contains additional boilerplate to make sure
#' elements are written correctly.
# nolint start: cyclocomp_linter
write_h5ad_element <- function(
  value,
  file,
  name,
  compression = c("none", "gzip", "lzf"),
  stop_on_error = TRUE,
  rhdf5 = FALSE,
  ...
) {
  if (rhdf5) {
    rhdf5_write_h5ad_element(
      value = value,
      file = file,
      name = name,
      compression = compression,
      stop_on_error = stop_on_error,
      ...
    )

    return(invisible(NULL))
  }

  compression <- match.arg(compression)

  # Sparse matrices
  write_fun <-
    if (inherits(value, "sparseMatrix")) {
      # Sparse matrices
      write_h5ad_sparse_array
    } else if (is.factor(value)) {
      # Categoricals
      write_h5ad_categorical
    } else if (is.list(value)) {
      # Lists and data frames
      if (is.data.frame(value)) {
        write_h5ad_data_frame
      } else {
        write_h5ad_mapping
      }
    } else if (is.character(value)) {
      # Character values
      if (length(value) == 1 && !is.matrix(value)) {
        write_h5ad_string_scalar
      } else {
        write_h5ad_string_array
      }
    } else if (is.numeric(value) || inherits(value, "denseMatrix")) {
      # Numeric values
      if (length(value) == 1 && !is.matrix(value)) {
        write_h5ad_numeric_scalar
      } else if (
        is.integer(value) && any(is.na(value)) && length(dim(value)) <= 1
      ) {
        write_h5ad_nullable_integer
      } else {
        write_h5ad_dense_array
      }
    } else if (is.logical(value)) {
      # Logical values
      if (any(is.na(value))) {
        write_h5ad_nullable_boolean
      } else if (length(value) == 1) {
        # Single Booleans should be written as numeric scalars
        write_h5ad_numeric_scalar
      } else {
        write_h5ad_dense_array
      }
    } else {
      # Fail if unknown
      cli_abort(c(
        "Writing {.cls {class(value)}} objects to H5AD is not supported",
        "i" = "Attempting to write to {.path {name}} in {.file {file}}"
      ))
    }

  # Delete the path if it already exists
  # TODO: do this here?
  if (hdf5_path_exists(file, name)) {
    hdf5r::h5unlink(file, name)
  }

  tryCatch(
    {
      write_fun(
        value = value,
        file = file,
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
# nolint end: cyclocomp_linter

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
write_h5ad_encoding <- function(file, name, encoding, version) {
  # interpreted from:
  # https://github.com/ycli1995/hdf5r.Extra/blob/a09078234a9457d74c019dabae8619393826b581/R/hdf5-internal.R#L76

  hdf5_create_attribute(
    file,
    name,
    "encoding-type",
    encoding,
    is_scalar = TRUE
  )
  hdf5_create_attribute(
    file,
    name,
    "encoding-version",
    version,
    is_scalar = TRUE
  )

  file
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
write_h5ad_dense_array <- function(
  value,
  file,
  name,
  compression,
  version = "0.2.0"
) {
  version <- match.arg(version)

  # matrices of type 'dgeMatrix' can simply be converted to a matrix
  if (inherits(value, "denseMatrix")) {
    value <- as.matrix(value)
  }

  if (is.matrix(value) && any(is.na(value))) {
    # is.na(value) <- NaN gets ignored
    na_indices <- is.na(value)
    value[na_indices] <- NaN
  }

  if (!is.vector(value)) {
    if (is.matrix(value)) {
      value <- t(value)
    } else if (is.array(value)) {
      value <- aperm(value)
    }
  }

  # Guess data type
  dtype <- NULL

  if (is.logical(value)) {
    dtype <- hdf5r::H5T_LOGICAL$new(include_NA = FALSE)
  }

  # Write dense array
  hdf5_create_dataset(
    file = file,
    name = name,
    value = value,
    compression = compression,
    dtype = dtype
  )

  # Write encoding
  write_h5ad_encoding(file, name, "array", version)
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
write_h5ad_sparse_array <- function(
  value,
  file,
  name,
  compression,
  version = "0.1.0"
) {
  version <- match.arg(version)

  # check types
  if (!(inherits(value, "sparseMatrix"))) {
    cli_abort(
      "{.arg value} must be a {.cls sparseMatrix} but has class {.cls {class(value)}}"
    )
  }

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
  file$create_group(name)
  write_h5ad_dense_array(
    attr(value, indices_attr),
    file,
    paste0(name, "/indices"),
    compression
  )
  write_h5ad_dense_array(value@p, file, paste0(name, "/indptr"), compression)
  write_h5ad_dense_array(value@x, file, paste0(name, "/data"), compression)
  write_h5ad_encoding(file, name, type, version)

  # Write shape attribute
  hdf5_create_attribute(
    file,
    name,
    "shape",
    dim(value),
    is_scalar = FALSE
  )
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
write_h5ad_nullable_boolean <- function(
  value,
  file,
  name,
  compression,
  version = "0.1.0"
) {
  file$create_group(name)

  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE

  write_h5ad_dense_array(
    value_no_na,
    file,
    paste0(name, "/values"),
    compression
  )

  # write mask manually
  # nolint start: commented_code_linter
  # write_h5ad_dense_array(is.na(value), file, paste0(name, "/mask"), compression)
  # nolint end: commented_code_linter

  # NOTE: `mask_dtype` will be written as a H5T_STD_U8LE, but h5py writes this as a H5T_STD_I8LE
  mask_dtype <- hdf5r::H5T_LOGICAL$new(include_NA = FALSE)

  hdf5_create_dataset(
    file = file,
    name = paste0(name, "/mask"),
    value = is.na(value),
    dtype = mask_dtype
  )
  write_h5ad_encoding(file, paste0(name, "/mask"), "array", "0.2.0")

  # set encoding
  write_h5ad_encoding(file, name, "nullable-boolean", version)
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
write_h5ad_nullable_integer <- function(
  value,
  file,
  name,
  compression,
  version = "0.1.0"
) {
  file$create_group(name)

  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- 1L

  write_h5ad_dense_array(
    value_no_na,
    file,
    paste0(name, "/values"),
    compression
  )
  # nolint start: commented_code_linter
  # write_h5ad_dense_array(is.na(value), file, paste0(name, "/mask"), compression)
  # nolint end: commented_code_linter

  # write mask manually
  mask_value <- is.na(value)
  mask_dims <- length(mask_value)
  mask_dtype <- hdf5r::H5T_LOGICAL$new(include_NA = FALSE)
  mask_space <- hdf5r::guess_space(
    mask_value,
    dtype = mask_dtype,
    chunked = FALSE
  )
  mask_gzip_level <- if (compression == "none") 0 else 9
  mask_name <- paste0(name, "/mask")
  file$create_dataset(
    name = mask_name,
    dims = mask_dims,
    gzip_level = mask_gzip_level,
    robj = mask_value,
    chunk_dims = NULL,
    space = mask_space,
    dtype = mask_dtype
  )
  write_h5ad_encoding(file, mask_name, "array", "0.2.0")
  write_h5ad_encoding(file, name, "nullable-integer", version)
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
write_h5ad_string_array <- function(
  value,
  file,
  name,
  compression,
  version = "0.2.0"
) {
  # TODO: add variable length string and encoding?
  hdf5_create_dataset(file, name, value, compression)

  write_h5ad_encoding(file, name, "string-array", version)
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
write_h5ad_categorical <- function(
  value,
  file,
  name,
  compression,
  version = "0.2.0"
) {
  file$create_group(name)

  categories <- levels(value)

  # Use zero-indexed values
  codes <- as.integer(value) - 1L

  # Set missing values to -1
  codes[is.na(codes)] <- -1L

  # write values to file
  write_h5ad_string_array(
    categories,
    file,
    paste0(name, "/categories"),
    compression
  )
  write_h5ad_dense_array(codes, file, paste0(name, "/codes"), compression)

  # Write encoding
  write_h5ad_encoding(
    file = file,
    name = name,
    encoding = "categorical",
    version = version
  )

  # Write ordered attribute
  hdf5_create_attribute(
    file,
    name,
    "ordered",
    is.ordered(value),
    is_scalar = TRUE
  )
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
write_h5ad_string_scalar <- function(
  value,
  file,
  name,
  compression,
  version = "0.2.0"
) {
  # TODO: add compression, variable length string and encoding!
  hdf5_create_dataset(file, name, value, compression, scalar = TRUE)

  # Write encoding
  write_h5ad_encoding(file, name, "string", version)
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
write_h5ad_numeric_scalar <- function(
  value,
  file,
  name,
  compression,
  version = "0.2.0"
) {
  # Write scalar
  hdf5_create_dataset(file, name, value, compression, scalar = TRUE)

  # Write encoding
  write_h5ad_encoding(file, name, "numeric-scalar", version)
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
write_h5ad_mapping <- function(
  value,
  file,
  name,
  compression,
  version = "0.1.0"
) {
  file$create_group(name)

  # Write mapping elements
  for (key in names(value)) {
    write_h5ad_element(value[[key]], file, paste0(name, "/", key), compression)
  }

  write_h5ad_encoding(file, name, "dict", version)
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
write_h5ad_data_frame <- function(
  value,
  file,
  name,
  compression,
  index = NULL,
  version = "0.2.0"
) {
  file$create_group(name)
  write_h5ad_encoding(file, name, "dataframe", version)

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
    write_h5ad_element(value[[col]], file, paste0(name, "/", col), compression)
  }

  # write index
  write_h5ad_element(
    index_value,
    file,
    paste0(name, "/", index_name),
    compression
  )

  # Write additional data frame attributes
  hdf5_create_attribute(
    file,
    name,
    "_index",
    index_name,
    is_scalar = TRUE
  )
  hdf5_create_attribute(
    file,
    name,
    "column-order",
    colnames(value),
    is_scalar = FALSE
  )
}

#' Write empty H5AD
#'
#' Write a new empty H5AD file
#'
#' @noRd
#'
#' @param file Path to the H5AD file to write
#' @param obs Data frame with observations
#' @param var Data frame with variables
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version The H5AD version to write
write_empty_h5ad <- function(file, obs, var, compression, version = "0.1.0") {
  write_h5ad_encoding(file, "/", "anndata", "0.1.0")

  write_h5ad_element(obs[, integer(0)], file, "/obs", compression)
  write_h5ad_element(var[, integer(0)], file, "/var", compression)

  file$create_group("layers")
  write_h5ad_encoding(file, "/layers", "dict", "0.1.0")

  file$create_group("obsm")
  write_h5ad_encoding(file, "/obsm", "dict", "0.1.0")

  file$create_group("obsp")
  write_h5ad_encoding(file, "/obsp", "dict", "0.1.0")

  file$create_group("uns")
  write_h5ad_encoding(file, "/uns", "dict", "0.1.0")

  file$create_group("varm")
  write_h5ad_encoding(file, "/varm", "dict", "0.1.0")

  file$create_group("varp")
  write_h5ad_encoding(file, "/varp", "dict", "0.1.0")

  invisible(NULL)
}
