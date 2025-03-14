skip_if_not_installed("rhdf5")
skip_if_not_installed("pizzarr")

# h5ad file
file <- hdf5r::H5File$new(system.file("extdata", "example.h5ad", package = "anndataR"), mode = "r")

# zarr file
zarr_dir <- system.file("extdata", "example.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(zarr_dir, exdir = td)
zarr_dir <- file.path(td, "example.zarr")
store <- pizzarr::DirectoryStore$new(zarr_dir)

test_that("reading dense matrices is same for h5ad and zarr", {
  mat_h5ad <- read_h5ad_dense_array(file, "layers/dense_counts")
  mat_zarr <- read_zarr_dense_array(store, "layers/dense_counts")
  expect_equal(mat_h5ad, mat_zarr)

  mat_h5ad <- read_h5ad_dense_array(file, "layers/dense_X")
  mat_zarr <- read_zarr_dense_array(store, "layers/dense_X")
  expect_equal(mat_h5ad, mat_zarr)
})

test_that("reading sparse matrices is same for h5ad and zarr", {
  mat_h5ad <- read_h5ad_sparse_array(file, "layers/csc_counts", type = "csc")
  mat_zarr <- read_zarr_sparse_array(store, "layers/csc_counts", type = "csc")
  expect_equal(mat_h5ad, mat_zarr)

  mat_h5ad <- read_h5ad_sparse_array(file, "layers/counts", type = "csr")
  mat_zarr <- read_zarr_sparse_array(store, "layers/counts", type = "csr")
  expect_equal(mat_h5ad, mat_zarr)
})

test_that("reading recarrays works", {
  skip("read_zarr_rec_array is not implemented yet")
  array_list <- read_zarr_rec_array(
    file, "uns/rank_genes_groups/logfoldchanges"
  )
  expect_true(is.list(array_list))
  expect_equal(names(array_list), c("0", "1", "2", "3", "4", "5"))
  for (array in array_list) {
    expect_true(is.array(array))
    expect_type(array, "double")
    expect_equal(dim(array), 100)
  }
})

test_that("reading 1D numeric arrays is same for h5ad and zarr", {
  array_1d_h5ad <- read_h5ad_dense_array(file, "obs/Int")
  array_1d_zarr <- read_zarr_dense_array(store, "obs/Int")
  expect_equal(array_1d_h5ad, array_1d_zarr)

  array_1d_h5ad <- read_h5ad_dense_array(file, "obs/Float")
  array_1d_zarr <- read_zarr_dense_array(store, "obs/Float")
  expect_equal(array_1d_h5ad, array_1d_zarr)
})

test_that("reading 1D sparse numeric arrays is same for h5ad and zarr", {
  array_1d_h5ad <- read_h5ad_sparse_array(file, "uns/Sparse1D", type = "csc")
  array_1d_zarr <- read_zarr_sparse_array(store, "uns/Sparse1D", type = "csc")
  expect_equal(array_1d_h5ad, array_1d_zarr)
})

test_that("reading 1D nullable arrays is same for h5ad and zarr", {
  array_1d_h5ad <- read_h5ad_nullable_integer(file, "obs/IntNA")
  array_1d_zarr <- read_zarr_nullable_integer(store, "obs/IntNA")
  expect_equal(array_1d_h5ad, array_1d_zarr)

  array_1d_h5ad <- read_h5ad_dense_array(file, "obs/FloatNA")
  array_1d_zarr <- read_zarr_dense_array(store, "obs/FloatNA")
  expect_equal(array_1d_h5ad, array_1d_zarr)

  # TODO: check this test, zarr Bools are stored as dense array hence no mask is given
  array_1d_h5ad <- read_h5ad_nullable_boolean(file, "obs/Bool")
  array_1d_zarr <- read_zarr_dense_array(store, "obs/Bool") # TODO: read_zarr_nullable_boolean should be used instead ?
  expect_equal(array_1d_h5ad, array_1d_zarr)

  array_1d_h5ad <- read_h5ad_nullable_boolean(file, "obs/BoolNA")
  array_1d_zarr <- read_zarr_nullable_boolean(store, "obs/BoolNA")
  expect_equal(array_1d_h5ad, array_1d_zarr)
})

test_that("reading string scalars is same for h5ad and zarr", {
  scalar_h5ad <- read_h5ad_string_scalar(file, "uns/StringScalar")
  scalar_zarr <- read_zarr_string_scalar(store, "uns/StringScalar")
  expect_equal(scalar_h5ad, scalar_zarr)
})

test_that("reading numeric scalars is same for h5ad and zarr", {
  scalar_h5ad <- read_h5ad_numeric_scalar(file, "uns/IntScalar")
  scalar_zarr <- read_zarr_numeric_scalar(store, "uns/IntScalar")
  expect_equal(scalar_h5ad, scalar_zarr)
})

test_that("reading string arrays is same for h5ad and zarr", {
  array_h5ad <- read_h5ad_string_array(file, "uns/String")
  array_zarr <- read_zarr_string_array(store, "uns/String")
  expect_equal(array_h5ad, array_zarr)

  array_h5ad <- read_h5ad_string_array(file, "uns/String2D")
  array_zarr <- read_zarr_string_array(store, "uns/String2D")
  expect_equal(array_h5ad, array_zarr)
})

test_that("reading mappings is same for h5ad and zarr", {
  skip("read_zarr_mapping returns list in different order")
  mapping_h5ad <- read_h5ad_mapping(file, "uns")
  mapping_zarr <- read_zarr_mapping(store, "uns")

  expect_equal(mapping_h5ad, mapping_zarr)
})

test_that("reading dataframes works", {
  df_h5ad <- read_h5ad_data_frame(file, "obs")
  df_zarr <- read_zarr_data_frame(store, "obs", include_index = TRUE)

  expect_equal(df_h5ad, df_zarr)
})

test_that("reading H5AD as SingleCellExperiment is same for h5ad and zarr", {
  skip_if_not_installed("SingleCellExperiment")

  sce_h5ad <- read_h5ad(file, to = "SingleCellExperiment")
  sce_zarr <- read_zarr(store, to = "SingleCellExperiment")

  expect_equal(sce_h5ad, sce_zarr)
})
