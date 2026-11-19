#' Maximum number of non-zero values in a sparse matrix from {Matrix}
#'
#' The `p` slot of a `dgCMatrix` or `dgRMatrix` is an integer vector, so these
#' classes can't hold more non-zero values than this. Larger matrices are stored
#' as an `SVT_SparseMatrix`, which has no such limit.
#'
#' This is a function so that tests can lower the limit, as matrices that are
#' really this large need tens of gigabytes of memory.
#'
#' @return The maximum number of non-zero values
#' @noRd
sparse_max_nnz <- function() {
  .Machine$integer.max
}

#' Number of non-zero values to process at once for large sparse matrices
#'
#' @return The target number of non-zero values in a block
#' @noRd
sparse_block_nnz <- function() {
  2^28
}

#' Split the major axis of a sparse matrix into blocks
#'
#' Only a block with a single row or column can be larger than
#' `sparse_block_nnz()`, which still fits in a `dgCMatrix` because a row or
#' column holds at most `.Machine$integer.max` values.
#'
#' @param pointers Index pointers, 0-based, as a double so that they can be
#'   larger than `.Machine$integer.max`
#'
#' @return A list of length 2 integer vectors, with the first and the last row
#'   or column of each block
#' @noRd
sparse_block_bounds <- function(pointers) {
  n_major <- length(pointers) - 1L
  block_nnz <- sparse_block_nnz()

  bounds <- list()
  first <- 1L
  while (first <= n_major) {
    last <- findInterval(pointers[[first]] + block_nnz, pointers) - 1L
    last <- min(max(last, first), n_major)
    bounds[[length(bounds) + 1L]] <- c(first, last)
    first <- last + 1L
  }
  bounds
}

#' Read a sparse matrix given a way to read parts of its arrays
#'
#' Matrices that fit in a `dgCMatrix` or `dgRMatrix` are read in one go. Larger
#' ones are read in blocks along the major axis into an `SVT_SparseMatrix`, so
#' that neither the reader nor R needs to hold the arrays as one long vector.
#'
#' @param read_slice A function with arguments `what` (`"data"` or `"indices"`),
#'   `start` (1-based) and `count`, which returns that slice as a vector, or the
#'   whole array if `start` and `count` are `NULL`
#' @param indptr Index pointers, 0-based, as a numeric so that they can be larger
#'   than `.Machine$integer.max`
#' @param shape Matrix dimensions
#' @param type Either `"csc_matrix"` or `"csr_matrix"`
#'
#' @return A `dgCMatrix`, `dgRMatrix` or `SVT_SparseMatrix`
#' @noRd
read_sparse_matrix <- function(read_slice, indptr, shape, type) {
  nnz <- if (length(indptr) > 0L) indptr[[length(indptr)]] else 0

  if (nnz <= sparse_max_nnz()) {
    return(construct_sparse_matrix(
      data = read_slice("data"),
      indices = read_slice("indices"),
      indptr = indptr,
      shape = shape,
      type = type
    ))
  }

  check_requires(
    "Reading a sparse matrix with more than 2^31 - 1 non-zero values",
    c("SparseArray", "BiocGenerics"),
    where = "Bioc"
  )
  if (is.unsorted(indptr)) {
    cli_abort("The index pointers of a sparse matrix must be increasing")
  }

  # An SVT_SparseMatrix is stored by column. Build the major axis as columns:
  # each block is a valid dgCMatrix, and cbind() reuses the columns of its input.
  n_minor <- if (type == "csc_matrix") shape[[1]] else shape[[2]]
  blocks <- lapply(sparse_block_bounds(indptr), function(bound) {
    first <- bound[[1L]]
    last <- bound[[2L]]
    offset <- indptr[[first]]
    count <- indptr[[last + 1L]] - offset

    read_block <- function(what) {
      if (count > 0) read_slice(what, start = offset + 1, count = count)
    }

    block <- construct_sparse_matrix(
      data = read_block("data"),
      indices = read_block("indices"),
      indptr = indptr[first:(last + 1L)] - offset,
      shape = c(n_minor, last - first + 1L),
      type = "csc_matrix"
    )
    as(block, "SVT_SparseMatrix")
  })
  svt <- do.call(BiocGenerics::cbind, blocks)

  if (type == "csc_matrix") svt else t(svt)
}

#' Realize a matrix so that it can be written
#'
#' @param value The value being written
#'
#' @return If `value` is a `DelayedArray` or `SparseMatrix` a `matrix`,
#'   `CsparseMatrix`, or an `SVT_SparseMatrix` if it has too many non-zero values
#'   for a `CsparseMatrix`. Any other `value` is returned as it is.
#' @noRd
realize_matrix_for_write <- function(value) {
  if (inherits(value, "DelayedArray")) {
    if (!DelayedArray::is_sparse(value)) {
      return(as.matrix(value))
    }
    value <- as(value, "SVT_SparseMatrix")
  }

  if (inherits(value, "SparseMatrix")) {
    value <- if (SparseArray::nzcount(value) <= sparse_max_nnz()) {
      as(value, "CsparseMatrix")
    } else {
      as(value, "SVT_SparseMatrix")
    }
  }

  value
}

#' Prepare writing an `SVT_SparseMatrix` in blocks
#'
#' @param value An `SVT_SparseMatrix`
#'
#' @return A list with `value` as a double matrix, its `indptr` as a double and
#'   the `bounds` of the blocks to write
#' @noRd
svt_write_layout <- function(value) {
  check_requires(
    "Writing a sparse matrix with more than 2^31 - 1 non-zero values",
    c("BiocGenerics", "MatrixGenerics"),
    where = "Bioc"
  )
  if (BiocGenerics::type(value) != "double") {
    BiocGenerics::type(value) <- "double"
  }
  indptr <- c(
    0,
    cumsum(MatrixGenerics::colSums(SparseArray::is_nonzero(value)))
  )

  list(
    value = value,
    indptr = unname(indptr),
    bounds = sparse_block_bounds(indptr)
  )
}

#' Get a block of columns of an `SVT_SparseMatrix` as a `CsparseMatrix`
#'
#' @param value An `SVT_SparseMatrix`
#' @param bound The first and last column of the block
#'
#' @return A `CsparseMatrix`
#' @noRd
svt_block <- function(value, bound) {
  as(value[, bound[[1L]]:bound[[2L]], drop = FALSE], "CsparseMatrix")
}

#' Convert 64-bit integers to and from doubles with the same bytes
#'
#' {Rarr} can only read and write 64-bit integers that fit in an R integer. A
#' 64-bit array read as `float64` keeps its bytes, whatever its codecs do with
#' them, and this converts those doubles to the value of the integer.
#'
#' Values can't be larger than `2^53`, the largest integer a double can hold
#' exactly. As an index pointer that would be far more than fits in memory.
#'
#' @param x A numeric vector of whole numbers (`int64_to_double_bits()`), or of
#'   doubles read as `float64` (`double_bits_to_int64()`)
#'
#' @return A numeric vector
#' @noRd
int64_to_double_bits <- function(x) {
  bytes <- vapply(
    256^(0:7),
    function(base) floor(x / base) %% 256,
    numeric(length(x))
  )
  dim(bytes) <- c(length(x), 8L)
  readBin(
    as.raw(t(bytes)),
    "double",
    n = length(x),
    size = 8L,
    endian = "little"
  )
}

#' @rdname int64_to_double_bits
#' @noRd
double_bits_to_int64 <- function(x) {
  bytes <- as.numeric(writeBin(x, raw(), size = 8L, endian = "little"))
  colSums(matrix(bytes, nrow = 8L) * 256^(0:7))
}

#' Read a Zarr array of 64-bit integers as doubles
#'
#' TODO: Replace with `Rarr::read_zarr_array()` when it can return 64-bit
#' integers that don't fit in an R integer. For now it makes those `NA`, see
#' https://github.com/Huber-group-EMBL/Rarr/issues/289.
#'
#' @param path Path to the Zarr array
#'
#' @return A numeric vector
#' @noRd
zarr_read_int64 <- function(path) {
  values <- withCallingHandlers(
    as.double(Rarr::read_zarr_array(path)),
    warning = function(w) {
      if (grepl("Integer overflow", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!anyNA(values)) {
    return(values)
  }

  metadata_file <- zarr_array_metadata_file(path)
  metadata <- readChar(metadata_file, file.size(metadata_file))

  float_metadata <- sub(
    '("dtype"\\s*:\\s*"[<>])[iu]8"',
    '\\1f8"',
    sub('("data_type"\\s*:\\s*")u?int64"', '\\1float64"', metadata)
  )
  if (identical(float_metadata, metadata)) {
    return(values)
  }

  # The array is read through a copy with different metadata. An indptr has 8
  # bytes for each row or column, so this copy is small.
  copy_dir <- tempfile("anndataR-int64")
  dir.create(copy_dir)
  on.exit(unlink(copy_dir, recursive = TRUE), add = TRUE)
  file.copy(path, copy_dir, recursive = TRUE)
  copy <- file.path(copy_dir, basename(path))
  writeChar(float_metadata, zarr_array_metadata_file(copy), eos = NULL)

  double_bits_to_int64(as.vector(Rarr::read_zarr_array(copy)))
}

#' Write a Zarr array of 64-bit integers
#'
#' TODO: Replace with `Rarr::write_zarr_array()` when it can write 64-bit
#' integers that don't fit in an R integer. For now it fails on doubles, see
#' https://github.com/Huber-group-EMBL/Rarr/issues/289.
#'
#' @param store Path to the Zarr store
#' @param name Name of the array in the store
#' @param value A numeric vector of whole numbers
#' @inheritParams write_zarr_element
#'
#' @noRd
zarr_write_int64 <- function(store, name, value, compression, zarr_format) {
  zarr_write_compressed(
    store,
    name,
    int64_to_double_bits(value),
    compression,
    zarr_format = zarr_format
  )

  metadata_file <- zarr_array_metadata_file(file.path(store, name))
  metadata <- readChar(metadata_file, file.size(metadata_file))
  metadata <- sub(
    '("dtype"\\s*:\\s*"[<>])f8"',
    '\\1i8"',
    sub('("data_type"\\s*:\\s*")float64"', '\\1int64"', metadata)
  )
  writeChar(metadata, metadata_file, eos = NULL)

  invisible()
}

#' Check that a backed read can use the index pointers of a Zarr sparse matrix
#'
#' A backed read goes through {ZarrArray}, which reads `indptr` with
#' `Rarr::read_zarr_array()`. That makes 64-bit integers that don't fit in an R
#' integer `NA`, and the read then fails with an error that doesn't say why. This
#' fails with an informative error instead. Only the last index pointer is read,
#' as it is the largest.
#'
#' TODO: Remove when {Rarr} can read these integers and bump the version in
#' `DESCRIPTION`, see https://github.com/Huber-group-EMBL/Rarr/issues/289.
#'
#' @param store Path to the Zarr store
#' @param name Name of the sparse matrix in the store
#' @param n_pointers Length of `indptr`
#' @param call The environment to report the error from
#'
#' @noRd
zarr_check_backed_indptr <- function(
  store,
  name,
  n_pointers,
  call = rlang::caller_env()
) {
  last_pointer <- withCallingHandlers(
    Rarr::read_zarr_array(
      file.path(store, name, "indptr"),
      index = list(n_pointers)
    ),
    warning = function(w) {
      if (grepl("Integer overflow", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!anyNA(last_pointer)) {
    return(invisible(TRUE))
  }

  cli_abort(
    c(
      paste(
        "Sparse matrices with more than 2^31 - 1 non-zero values can't be",
        "read from Zarr with {.code backed = TRUE} yet, as {.pkg Rarr} can't",
        "read their 64-bit {.field indptr}",
        "({.url https://github.com/Huber-group-EMBL/Rarr/issues/289})"
      ),
      "i" = paste(
        "Read the matrix into memory with {.code backed = FALSE}, or use H5AD",
        "instead, see",
        "{.url https://anndataR.scverse.org/articles/large_matrices.html}"
      )
    ),
    call = call
  )
}

#' Get the path of the metadata file of a Zarr array
#'
#' @param path Path to the Zarr array
#'
#' @return The path of `zarr.json` for Zarr v3 and `.zarray` for Zarr v2
#' @noRd
zarr_array_metadata_file <- function(path) {
  metadata_file <- file.path(path, "zarr.json")
  if (file.exists(metadata_file)) metadata_file else file.path(path, ".zarray")
}

#' Convert a large Python sparse matrix to R
#'
#' `reticulate::py_to_r()` builds a `dgCMatrix` or `dgRMatrix`, so it can't
#' convert a matrix with too many non-zero values for those. Those matrices are
#' converted to an `SVT_SparseMatrix` instead.
#'
#' @param x A Python object
#'
#' @return The R version of `x`
#' @noRd
py_to_r_matrix <- function(x) {
  if (!inherits(x, "python.builtin.object")) {
    return(py_to_r(x))
  }

  bi <- reticulate::import_builtins(convert = FALSE)
  is_sparse <- reticulate::py_to_r(bi$hasattr(x, "nnz")) &&
    reticulate::py_to_r(bi$hasattr(x, "tocsr"))
  # A Python int above .Machine$integer.max overflows in reticulate, a float doesn't
  if (
    !is_sparse ||
      reticulate::py_to_r(bi$float(bi$getattr(x, "nnz"))) <= sparse_max_nnz()
  ) {
    return(py_to_r(x))
  }

  format <- reticulate::py_to_r(bi$getattr(x, "format"))
  if (!format %in% c("csr", "csc")) {
    x <- bi$getattr(x, "tocsr")()
    format <- "csr"
  }

  read_slice <- function(what, start = NULL, count = NULL) {
    values <- bi$getattr(x, what)
    if (!is.null(start)) {
      values <- values$`__getitem__`(
        bi$slice(bi$int(start - 1), bi$int(start - 1 + count))
      )
    }
    as.vector(reticulate::py_to_r(values))
  }

  read_sparse_matrix(
    read_slice = read_slice,
    indptr = as.vector(reticulate::py_to_r(
      bi$getattr(x, "indptr")$astype("float64")
    )),
    shape = as.integer(unlist(reticulate::py_to_r(bi$getattr(x, "shape")))),
    type = paste0(format, "_matrix")
  )
}
