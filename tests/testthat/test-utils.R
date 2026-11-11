col_matrix <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = 0.1)
row_matrix <- as(col_matrix, "RsparseMatrix")

# we expect to receive a transposed dgRMatrix from to_py_matrix
test_that("to_py_matrix works with dgCMatrix", {
  result <- to_py_matrix(col_matrix)
  expect_s4_class(result, "dgRMatrix")
  expect_equal(result, Matrix::t(row_matrix))
})

test_that("to_py_matrix works with dgRMatrix", {
  result <- to_py_matrix(row_matrix)
  expect_s4_class(result, "dgRMatrix")
  expect_equal(result, Matrix::t(row_matrix))
})

# we expect to receive a transposed dgCMatrix from to_R_matrix
test_that("to_R_matrix works with dgRMatrix", {
  result <- to_R_matrix(row_matrix)
  expect_s4_class(result, "dgCMatrix")
  expect_equal(result, Matrix::t(col_matrix))
})

test_that("to_R_matrix works with dgCMatrix", {
  result <- to_R_matrix(col_matrix)
  expect_s4_class(result, "dgCMatrix")
  expect_equal(result, Matrix::t(col_matrix))
})

test_that("materialize_matrix returns in-memory matrices unchanged", {
  expect_identical(materialize_matrix(col_matrix), col_matrix)
  expect_identical(materialize_matrix(row_matrix), row_matrix)
})

test_that("materialize_matrix keeps sparse DelayedMatrix objects sparse", {
  skip_if_not_installed("DelayedArray")

  result <- materialize_matrix(DelayedArray::DelayedArray(row_matrix))
  expect_s4_class(result, "dgCMatrix")
  expect_equal(as.matrix(result), as.matrix(col_matrix), ignore_attr = TRUE)
})

test_that("materialize_matrix realizes dense DelayedMatrix objects", {
  skip_if_not_installed("DelayedArray")

  dense_matrix <- as.matrix(col_matrix)
  result <- materialize_matrix(DelayedArray::DelayedArray(dense_matrix))
  expect_identical(result, dense_matrix)
})

test_that("construct_sparse_matrix errors on out-of-range indices", {
  # 1-based row indices for a 3 x 4 CSC matrix
  expect_error(
    construct_sparse_matrix(
      data = rep(1, 4),
      indices = c(1L, 2L, 3L, 1L),
      indptr = c(0L, 1L, 2L, 3L, 4L),
      shape = c(3L, 4L),
      type = "csc_matrix"
    ),
    "must be 0-based"
  )

  # 1-based column indices for a 3 x 4 CSR matrix
  expect_error(
    construct_sparse_matrix(
      data = rep(1, 3),
      indices = c(1L, 2L, 4L),
      indptr = c(0L, 1L, 2L, 3L),
      shape = c(3L, 4L),
      type = "csr_matrix"
    ),
    "must be 0-based"
  )
})
