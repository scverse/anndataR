test_that("has_row_names works on matrix", {
  expect_false(has_row_names(matrix(1:26, ncol = 1)))
  expect_true(has_row_names(matrix(
    1:26,
    ncol = 1,
    dimnames = list(letters, NULL)
  )))
})

test_that("has_row_names works on data.frame", {
  expect_false(has_row_names(data.frame(x = 1:26)))
  expect_true(has_row_names(data.frame(x = 1:26, row.names = letters)))
})
