#' Maximum number of non-zero values in a sparse matrix from {Matrix}
#'
#' The `p` slot of a `dgCMatrix` or `dgRMatrix` is an integer vector, so these
#' classes can't hold more non-zero values than this. Larger matrices are stored
#' as an `SVT_SparseMatrix`, which has no such limit.
#'
#' This is a function so that tests can lower the limit, as matrices that are
#' really this large need tens of gigabytes of memory.
#'
#' TODO: Revisit when {Matrix} supports a 64-bit `p` slot, see
#' https://r-forge.r-project.org/tracker/index.php?func=detail&aid=6611&group_id=61&atid=297
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
      indptr = indptr[seq(first, last + 1L)] - offset,
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
    value <- as_csparse_or_svt(value)
  }

  value
}

#' Convert a `SparseMatrix` to a `CsparseMatrix` if it fits in one
#'
#' @param x A `SparseMatrix`
#'
#' @return A `CsparseMatrix`, or an `SVT_SparseMatrix` if `x` has too many
#'   non-zero values for a `CsparseMatrix`
#' @noRd
as_csparse_or_svt <- function(x) {
  if (SparseArray::nzcount(x) <= sparse_max_nnz()) {
    as(x, "CsparseMatrix")
  } else {
    as(x, "SVT_SparseMatrix")
  }
}

#' Prepare writing an `SVT_SparseMatrix` in blocks
#'
#' @param value An `SVT_SparseMatrix`
#'
#' @return A list with `value`, its `indptr` as a double, its number of non-zero
#'   values `nnz` and the `bounds` of the blocks to write
#' @noRd
svt_write_layout <- function(value) {
  check_requires(
    "Writing a sparse matrix with more than 2^31 - 1 non-zero values",
    "MatrixGenerics",
    where = "Bioc"
  )
  indptr <- c(
    0,
    cumsum(MatrixGenerics::colSums(SparseArray::is_nonzero(value)))
  )

  list(
    value = value,
    indptr = unname(indptr),
    nnz = indptr[[length(indptr)]],
    bounds = sparse_block_bounds(indptr)
  )
}

#' Write the `indices` and `data` of an `SVT_SparseMatrix` in blocks
#'
#' The counterpart of `read_slice` in `read_sparse_matrix()`. Each block of
#' columns is converted to a `CsparseMatrix`, so only one block is held twice.
#' The `data` is written as doubles, whatever the type of the matrix.
#'
#' @param layout The result of `svt_write_layout()`
#' @param write_slice A function with arguments `what` (`"indices"` or
#'   `"data"`), `start` (1-based) and `values`, which writes `values` to that
#'   array from position `start`
#'
#' @noRd
svt_write_blocks <- function(layout, write_slice) {
  for (bound in layout$bounds) {
    columns <- bound[[1L]]:bound[[2L]]
    block <- as(layout$value[, columns, drop = FALSE], "CsparseMatrix")
    if (length(block@x) > 0L) {
      start <- layout$indptr[[bound[[1L]]]] + 1
      write_slice("indices", start, block@i)
      # a logical matrix gives an lgCMatrix
      write_slice("data", start, as.double(block@x))
    }
  }
  invisible()
}

#' Convert 64-bit integers to doubles with the same bytes
#'
#' {Rarr} can't write doubles to a 64-bit integer array. Written as `float64`,
#' these doubles have the bytes of the integers, so changing the data type in
#' the metadata afterwards gives the 64-bit integers.
#'
#' Values can't be larger than `2^53`, the largest integer a double can hold
#' exactly. As an index pointer that would be far more than fits in memory.
#'
#' @param x A numeric vector of whole numbers
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
