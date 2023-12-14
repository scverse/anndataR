skip_if_not_installed("hdf5r")

file <- hdf5r::H5File$new(system.file("extdata", "example.h5ad", package = "anndataR"), mode = "r")
on.exit(file$close())

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

test_that("reading recarrays works", {
  array_list <- read_h5ad_rec_array(
    file, "uns/rank_genes_groups/logfoldchanges"
  )
  expect_true(is.list(array_list))
  expect_equal(names(array_list), c("0", "1", "2", "3", "4", "5"))
  for (array in array_list) {
    expect_true(is.vector(array))
    expect_type(array, "double")
    expect_equal(length(array), 100)
  }
})

test_that("reading 1D numeric arrays works", {
  array_1d <- read_h5ad_dense_array(file, "obs/Int")
  expect_vector(array_1d, ptype = integer(), size = 50)

  array_1d <- read_h5ad_dense_array(file, "obs/Float")
  expect_vector(array_1d, ptype = double(), size = 50)
})

test_that("reading 1D sparse numeric arrays works", {
  array_1d <- read_h5ad_sparse_array(file, "uns/Sparse1D", type = "csc")
  expect_s4_class(array_1d, "dgCMatrix")
  expect_equal(dim(array_1d), c(1, 6))
})

test_that("reading 1D nullable arrays works", {
  array_1d <- read_h5ad_nullable_integer(file, "obs/IntNA")
  expect_vector(array_1d, ptype = integer(), size = 50)
  expect_true(any(is.na(array_1d)))

  array_1d <- read_h5ad_dense_array(file, "obs/FloatNA")
  expect_vector(array_1d, ptype = double(), size = 50)
  expect_true(any(is.na(array_1d)))

  array_1d <- read_h5ad_nullable_boolean(file, "obs/Bool")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_false(any(is.na(array_1d)))

  array_1d <- read_h5ad_nullable_boolean(file, "obs/BoolNA")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_true(any(is.na(array_1d)))
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
  expect_equal(dim(array), c(5, 10))
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
    c(
      "Float", "FloatNA", "Int", "IntNA", "Bool", "BoolNA",
      "n_genes_by_counts", "log1p_n_genes_by_counts", "total_counts",
      "log1p_total_counts", "leiden"
    )
  )
})

test_that("reading H5AD as SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_h5ad(file, to = "SingleCellExperiment")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading H5AD as Seurat works", {
  skip_if_not_installed("SeuratObject")

  # TODO: remove this suppression when the to_seurat, from_seurat functions are updated.
  seurat <- suppressWarnings(read_h5ad(file, to = "Seurat"))
  expect_s4_class(seurat, "Seurat")
})
