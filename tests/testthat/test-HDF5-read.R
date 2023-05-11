file <- system.file("extdata", "example.h5ad", package = "anndataR")

# Should contain only `type` and `version`
test_that("reading encoding works", {
  encoding <- read_h5ad_encoding(file, "obs")
  expect_equal(names(encoding), c("type", "version"))
})

test_that("reading dense matrices works", {
  mat <- read_h5ad_dense_array(file, "layers/dense_counts")
  expect_true(is.matrix(mat))
  expect_type(mat, "integer")
  expect_equal(dim(mat), c(50, 100))

  mat <- read_h5ad_dense_array(file, "layers/dense_X")
  expect_true(is.matrix(mat))
  expect_type(mat, "double")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading sparse matrices works", {
  mat <- read_h5ad_sparse_array(file, "layers/csc_counts", type = "csc")
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(50, 100))
  
  mat <- read_h5ad_sparse_array(file, "layers/counts", type = "csr")
  expect_s4_class(mat, "dgRMatrix")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading 1D numeric arrays works", {
  array1D <- read_h5ad_dense_array(file, "obs/Int")
  expect_vector(array1D, ptype = integer(), size = 50)
  
  array1D <- read_h5ad_dense_array(file, "obs/Float")
  expect_vector(array1D, ptype = double(), size = 50)
})

# test_that("reading 1D sparse numeric arrays works", {
#   array1D <- read_h5ad_dense_array(file, "obs/dummy_num")
#   expect_equal(matrix1d, rep(42.42, 640))
# })

test_that("reading 1D nullable arrays works", {
  array1D <- read_h5ad_nullable_integer(file, "obs/IntNA")
  expect_vector(array1D, ptype = integer(), size = 50)
  expect_true(any(is.na(array1D)))
  
  array1D <- read_h5ad_dense_array(file, "obs/FloatNA")
  expect_vector(array1D, ptype = double(), size = 50)
  expect_true(any(is.na(array1D)))
  
  array1D <- read_h5ad_nullable_boolean(file, "obs/Bool")
  expect_vector(array1D, ptype = logical(), size = 50)
  expect_false(any(is.na(array1D)))
  
  array1D <- read_h5ad_nullable_boolean(file, "obs/BoolNA")
  expect_vector(array1D, ptype = logical(), size = 50)
  expect_true(any(is.na(array1D)))
})

test_that("reading string scalars works", {
  scalar <- read_h5ad_string_scalar(file, "uns/StringScalar")
  expect_equal(scalar, "A string")
})

test_that("reading numeric scalars works", {
  scalar <- read_h5ad_numeric_scalar(file, "uns/IntScalar")
  expect_equal(scalar, 1)
})

test_that("reading string arrays works", {
  array <- read_h5ad_string_array(file, "uns/String")
  expect_vector(array, ptype = character(), size = 10)
  expect_equal(array[3], "String 2")
  
  array <- read_h5ad_string_array(file, "uns/String2D")
  expect_true(is.matrix(array))
  expect_type(array, "character")
  expect_equal(dim(array), c(10, 20))
})

test_that("reading mappings works", {
  mapping <- read_h5ad_mapping(file, "uns")
  expect_type(mapping, "list")
  expect_type(names(mapping), "character")
})

test_that("reading dataframes works", {
  df <- read_h5ad_data_frame(file, "obs")
  expect_s3_class(df, "data.frame")
  expect_equal(
    colnames(df),
    c("Float", "FloatNA", "Int", "IntNA", "Bool", "BoolNA", "n_genes_by_counts",
      "log1p_n_genes_by_counts", "total_counts", "log1p_total_counts", "leiden")
  )
})
