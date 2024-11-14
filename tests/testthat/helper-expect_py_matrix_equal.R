expect_py_matrix_equal <- function(a, b, rtol = 1e-6, atol = 1e-6) {
  requireNamespace("testthat")
  requireNamespace("reticulate")

  bi <- reticulate::import_builtins()
  np <- reticulate::import("numpy")
  scipy <- reticulate::import("scipy")

  testthat::expect_equal(bi$type(a), bi$type(b)) # does this always work?

  testthat::expect_equal(a$dtype, b$dtype)

  testthat::expect_equal(reticulate::py_to_r(a$shape), reticulate::py_to_r(b$shape))

  a_dense <-
    if (scipy$sparse$issparse(a)) {
      a$toarray()
    } else {
      a
    }
  b_dense <-
    if (scipy$sparse$issparse(b)) {
      b$toarray()
    } else {
      b
    }

  testthat::expect_null(
    np$testing$assert_allclose(a_dense, b_dense, rtol = rtol, atol = atol)
  )
}
