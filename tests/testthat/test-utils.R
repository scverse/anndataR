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
