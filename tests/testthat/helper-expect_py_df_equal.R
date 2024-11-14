expect_py_df_equal <- function(a, b) {
  pd <- reticulate::import("pandas")

  expect_null(
    pd$testing$assert_frame_equal(
      a,
      b,
      check_dtype = FALSE,
      check_exact = FALSE
    )
  )
}
