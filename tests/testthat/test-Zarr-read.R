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

  mat <- read_zarr_dense_array(store, "layers/dense_X")
  expect_true(is.matrix(mat))
  expect_type(mat, "double")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading sparse matrices works", {
  mat <- read_zarr_sparse_array(store, "layers/csc_counts", type = "csc")
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(50, 100))

  mat <- read_zarr_sparse_array(store, "layers/counts", type = "csr")
  expect_s4_class(mat, "dgRMatrix")
  expect_equal(dim(mat), c(50, 100))
})


test_that("reading recarrays works", {
  f <- function() read_zarr_rec_array(store, "uns/rank_genes_groups/logfoldchanges")
  expect_error(f())
})

test_that("reading 1D numeric arrays works", {
  array_1d <- read_zarr_dense_array(store, "obs/Int")
  expect_vector(array_1d, ptype = integer(), size = 50)

  array_1d <- read_zarr_dense_array(store, "obs/Float")
  expect_vector(array_1d, ptype = double(), size = 50)
})

test_that("reading 1D sparse numeric arrays works", {
  array_1d <- read_zarr_sparse_array(store, "uns/Sparse1D", type = "csc")
  expect_s4_class(array_1d, "dgCMatrix")
  expect_equal(dim(array_1d), c(1, 6))
})

test_that("reading 1D nullable arrays works", {
  array_1d <- read_zarr_nullable_integer(store, "obs/IntNA")
  expect_vector(array_1d, ptype = integer(), size = 50)
  expect_true(any(is.na(array_1d)))

  array_1d <- read_zarr_dense_array(store, "obs/FloatNA")
  expect_vector(array_1d, ptype = double(), size = 50)
  expect_true(any(is.na(array_1d)))

  array_1d <- read_zarr_nullable_boolean(store, "obs/BoolNA")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_true(any(is.na(array_1d)))
})

test_that("reading 1D nullable arrays works (Nullable boolean)", {
  skip("TODO: non NA booleans dont have mask arrays, should they ?")
  array_1d <- read_zarr_nullable_boolean(store, "obs/Bool")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_false(any(is.na(array_1d)))
})

test_that("reading string scalars works", {
  scalar <- read_zarr_string_scalar(store, "uns/StringScalar")
  expect_equal(scalar, "A string")
})

test_that("reading numeric scalars works", {
  scalar <- read_zarr_numeric_scalar(store, "uns/IntScalar")
  expect_equal(scalar, 1)
})

test_that("reading string arrays works", {
  array <- read_zarr_string_array(store, "uns/String")
  expect_vector(array, ptype = character(), size = 10)
  expect_equal(array[3], "String 2")

  array <- read_zarr_string_array(store, "uns/String2D")
  expect_true(is.matrix(array))
  expect_type(array, "character")
  expect_equal(dim(array), c(5, 10))
})

test_that("reading mappings works", {
  mapping <- read_zarr_mapping(store, "uns")
  expect_type(mapping, "list")
  expect_type(names(mapping), "character")
})

test_that("reading dataframes works", {
  df <- read_zarr_data_frame(store, "obs", include_index = TRUE)
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

test_that("reading Zarr as SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_zarr(store, to = "SingleCellExperiment")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading Zarr as Seurat works", {
  skip_if_not_installed("SeuratObject")

  # TODO: remove this suppression when the to_seurat, from_seurat functions are updated.
  seurat <- suppressWarnings(read_zarr(store, to = "Seurat"))
  expect_s4_class(seurat, "Seurat")
})
