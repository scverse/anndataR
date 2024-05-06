skip_if_not_installed("pizzarr")

file <- system.file("extdata", "example.zarr", package = "anndataR")
store <- pizzarr::DirectoryStore$new(file)

test_that("reading encoding works", {
  encoding <- read_zarr_encoding(store, "obs")
  expect_equal(names(encoding), c("type", "version"))
})

test_that("reading dense matrices works", {
  mat <- read_zarr_dense_array(store, "layers/dense_counts")
  expect_true(is.matrix(mat))
  expect_type(mat, "integer")
  expect_equal(dim(mat), c(50, 100))

  mat <- read_zarr_dense_array(file, "layers/dense_X")
  expect_true(is.matrix(mat))
  expect_type(mat, "double")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading sparse matrices works", {
  mat <- read_zarr_sparse_array(file, "layers/csc_counts", type = "csc")
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(50, 100))

  mat <- read_zarr_sparse_array(file, "layers/counts", type = "csr")
  expect_s4_class(mat, "dgRMatrix")
  expect_equal(dim(mat), c(50, 100))
})
