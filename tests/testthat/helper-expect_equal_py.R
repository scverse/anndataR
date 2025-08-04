expect_equal_py <- function(a, b) {
  requireNamespace("testthat")
  requireNamespace("reticulate")

  bi <- reticulate::import_builtins()

  # ↓ does this always work?
  testthat::expect_equal(bi$str(bi$type(a)), bi$str(bi$type(b)))

  if (inherits(a, "pandas.core.frame.DataFrame")) {
    pd <- reticulate::import("pandas")
    testthat::expect_null(
      pd$testing$assert_frame_equal(
        a,
        b,
        check_dtype = FALSE,
        check_exact = FALSE
      )
    )
  } else if (
    inherits(a, "np.ndarray") || inherits(a, "scipy.sparse.base.spmatrix")
  ) {
    scipy <- reticulate::import("scipy")
    np <- reticulate::import("numpy")

    testthat::expect_equal(a$dtype, b$dtype)

    testthat::expect_equal(
      py_to_r_ifneedbe(a$shape),
      py_to_r_ifneedbe(b$shape)
    )

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
      np$testing$assert_allclose(a_dense, b_dense)
    )
  }
}

py_to_r_ifneedbe <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    requireNamespace("reticulate")
    reticulate::py_to_r(x)
  } else {
    x
  }
}
