# should take a look at
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html
# again

#' Write H5AD element
#'
#' Write an element to an H5AD file
#'
#' @param value The value to write
#' @param file An open H5AD handle
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
  stop_on_error = FALSE,
  ...
) {
  compression <- match.arg(compression)

  # Sparse matrices
  write_fun <-
    if (is.null(value)) {
      write_h5ad_null
    } else if (inherits(value, "sparseMatrix")) {
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
    rhdf5::h5delete(file, name)
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
  hdf5_write_attribute(
    file,
    name,
    "encoding-type",
    encoding,
    is_scalar = TRUE
  )

  hdf5_write_attribute(
    file,
    name,
    "encoding-version",
    version,
    is_scalar = TRUE
  )
}

#' Write H5AD null
#'
#' Write a null dataset to an H5AD file
#'
#' @param value Value to write, not used
#' @param file An open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression Not used as there is no value
#' @param version Encoding version of the element to write
#'
#' @noRd
write_h5ad_null <- function(value, file, name, compression, version = "0.1.0") {
  if (isFALSE(getOption("anndataR.write_null", "TRUE"))) {
    return(invisible(NULL))
  }

  h5s <- rhdf5::H5Screate("H5S_NULL")
  on.exit(rhdf5::H5Sclose(h5s), add = TRUE)

  h5d <- rhdf5::H5Dcreate(
    file,
    name = name,
    dtype_id = "H5T_IEEE_F32LE",
    h5space = h5s
  )
  on.exit(rhdf5::H5Dclose(h5d), add = TRUE)

  write_h5ad_encoding(file, name, "null", version)
}

#' Write H5AD dense array
#'
#' Write a dense array to an H5AD file
#'
#' @param value Value to write
#' @param file An open H5AD handle
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

  H5type <- if (is.integer(value)) {
    "H5T_STD_I64LE"
  } else {
    NULL
  }

  # Write dense array
  if (is.logical(value)) {
    hdf5_write_boolean_dataset(
      file = file,
      name = name,
      value = value,
      compression = compression
    )
  } else {
    hdf5_write_dataset(
      file = file,
      name = name,
      value = value,
      H5type = H5type,
      compression = compression
    )
  }

  write_h5ad_encoding(file, name, "array", version)
}

#' Write H5AD sparse array
#'
#' Write a sparse array to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file An open H5AD handle
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
  rhdf5::h5createGroup(file, name)
  hdf5_write_dataset(
    file = file,
    name = paste0(name, "/indices"),
    value = attr(value, indices_attr),
    compression = compression
  )
  hdf5_write_dataset(
    file = file,
    name = paste0(name, "/indptr"),
    value = value@p,
    compression = compression
  )
  hdf5_write_dataset(
    file = file,
    name = paste0(name, "/data"),
    value = value@x,
    compression = compression
  )
  write_h5ad_encoding(file, name, type, version)

  # Write shape attribute
  hdf5_write_attribute(
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
#' @param file An open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
# nolint start: object_length_linter
write_h5ad_nullable_boolean <- function(
  value,
  file,
  name,
  compression,
  version = "0.1.0"
) {
  # nolint end: object_length_linter
  rhdf5::h5createGroup(file, name)

  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- FALSE

  write_h5ad_dense_array(
    value_no_na,
    file,
    paste0(name, "/values"),
    compression
  )

  write_h5ad_dense_array(
    is.na(value),
    file,
    paste0(name, "/mask"),
    compression
  )

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
#' @param file An open H5AD handle
#' @param name Name of the element within the H5AD file
#' @param compression The compression to use when writing the element. Can be
#' one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.
#' @param version Encoding version of the element to write
# nolint start: object_length_linter
write_h5ad_nullable_integer <- function(
  value,
  file,
  name,
  compression,
  version = "0.1.0"
) {
  # nolint end: object_length_linter
  rhdf5::h5createGroup(file, name)

  value_no_na <- value
  value_no_na[is.na(value_no_na)] <- 0L

  write_h5ad_dense_array(
    value_no_na,
    file,
    paste0(name, "/values"),
    compression
  )

  write_h5ad_dense_array(
    is.na(value),
    file,
    paste0(name, "/mask"),
    compression
  )

  write_h5ad_encoding(file, name, "nullable-integer", version)
}

#' Write H5AD string array
#'
#' Write a string array to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file An open H5AD handle
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
  if (!is.vector(value)) {
    if (is.matrix(value)) {
      value <- t(value)
    } else if (is.array(value)) {
      value <- aperm(value)
    }
  }

  hdf5_write_dataset(
    file = file,
    name = name,
    value = value,
    compression = compression
  )

  write_h5ad_encoding(file, name, "string-array", version)
}

#' Write H5AD categorical
#'
#' Write a categorical to an H5AD file
#'
#' @noRd
#'
#' @param value Value to write
#' @param file An open H5AD handle
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
  rhdf5::h5createGroup(file, name)

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
  hdf5_write_attribute(
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
#' @param file An open H5AD handle
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
  hdf5_write_scalar(
    file = file,
    name = name,
    value = value
  )

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
#' @param file An open H5AD handle
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
  if (is.logical(value)) {
    hdf5_write_boolean_dataset(
      file = file,
      name = name,
      value = value,
      is_scalar = TRUE,
      compression = compression
    )
  } else {
    hdf5_write_scalar(
      file = file,
      name = name,
      value = value
    )
  }

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
#' @param file An open H5AD handle
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
  rhdf5::h5createGroup(file, name)

  # Write mapping elements
  for (key in names(value)) {
    write_h5ad_element(
      value[[key]],
      file,
      paste0(name, "/", key),
      compression
    )
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
#' @param file An open H5AD handle
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
  rhdf5::h5createGroup(file, name)
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
    write_h5ad_element(
      value[[col]],
      file,
      paste0(name, "/", col),
      compression
    )
  }

  # write index
  write_h5ad_element(
    index_value,
    file,
    paste0(name, "/", index_name),
    compression
  )

  # Write additional data frame attributes
  hdf5_write_attribute(
    file,
    name,
    "_index",
    index_name,
    is_scalar = TRUE
  )

  col_order <- colnames(value)
  col_order <- col_order[col_order != index_name]
  # If there are no columns other than the index we set column order to an
  # empty numeric vector
  if (length(col_order) == 0) {
    col_order <- numeric()
  }

  hdf5_write_attribute(
    file,
    name,
    "column-order",
    col_order,
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
write_empty_h5ad <- function(
  file,
  obs,
  var,
  compression,
  version = "0.1.0"
) {
  write_h5ad_encoding(file, "/", "anndata", "0.1.0")

  write_h5ad_element(obs[, integer(0)], file, "/obs", compression)
  write_h5ad_element(var[, integer(0)], file, "/var", compression)

  rhdf5::h5createGroup(file, "layers")
  write_h5ad_encoding(file, "/layers", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "obsm")
  write_h5ad_encoding(file, "/obsm", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "obsp")
  write_h5ad_encoding(file, "/obsp", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "uns")
  write_h5ad_encoding(file, "/uns", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "varm")
  write_h5ad_encoding(file, "/varm", "dict", "0.1.0")

  rhdf5::h5createGroup(file, "varp")
  write_h5ad_encoding(file, "/varp", "dict", "0.1.0")

  invisible(NULL)
}
