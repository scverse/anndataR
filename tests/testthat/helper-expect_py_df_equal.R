expect_py_df_equal <- function(a, b) {
  requireNamespace("testthat")
  requireNamespace("reticulate")

  bi <- reticulate::import_builtins()
  pd <- reticulate::import("pandas")

  testthat::expect_equal(bi$type(a), bi$type(b)) # does this always work?

  testthat::expect_null(
    pd$testing$assert_frame_equal(
      a,
      b,
      check_dtype = FALSE,
      check_exact = FALSE
    )
  )
}
