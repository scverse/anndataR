file <- system.file("extdata", "example.h5ad", package = "anndataR")

# Should contain only `type` and `version`
test_that("reading encoding", {
  encoding <- read_h5ad_encoding(file, "obs")
  expect_equal(names(encoding), c("type", "version"))
})

test_that("reading dense matrix", {
  matrix2d <- read_h5ad_dense_array(file, "obsm/dense")
  # TODO: test dimensions, test a random coordinate to ensure the dimensions didn't flip
  expect_equal(dim(matrix2d), c(640, 10))
  expect_equal(matrix2d[43, 6], 47)

  matrix1d <- read_h5ad_dense_array(file, "obs/dummy_num")
  expect_equal(matrix1d, rep(42.42, 640))
})

test_that("reading sparse matrix", {
  matrix2d_csr <- read_h5ad_sparse_array(file, "obsm/sparse_csr", type = "csr")
  matrix2d_csc <- read_h5ad_sparse_array(file, "obsm/sparse_csc", type = "csc")
  expect_equal(dim(matrix2d_csr), c(640, 10))
  expect_equal(dim(matrix2d_csc), c(640, 10))

  expect_equal(matrix2d_csr[43, 6], 47)
  expect_equal(matrix2d_csc[43, 6], 47)

  matrix1d_csr <- read_h5ad_sparse_array(file, "obsm/sparse_csr_1", type = "csr")
  matrix1d_csc <- read_h5ad_sparse_array(file, "obsm/sparse_csc_1", type = "csc")

  expect_equal(matrix1d_csr[15], 14)
  expect_equal(matrix1d_csc[15], 14)
})

test_that("reading string scalar", {
  scalar <- read_h5ad_string_scalar(file, "uns/dummy_string_scalar")
  expect_equal(scalar, "foo")
})

test_that("reading numeric scalar", {
  scalar <- read_h5ad_numeric_scalar(file, "uns/dummy_int_scalar")
  expect_equal(scalar, 1)
})

test_that("reading string array", {
  sarray <- read_h5ad_string_array(file, "uns/dummy_string_array")
  expect_equal(sarray[43, 6], "row542")
})

test_that("reading mapping", {
  mapping <- read_h5ad_mapping(file, "uns")
  expect_type(mapping, "list")
  expect_type(names(mapping), "character")
})

test_that("reading dataframe", {
  df <- read_h5ad_data_frame(file, "obs")
  expect_true(is.data.frame(df))
})
