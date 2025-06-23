skip_if_not_installed("rhdf5")

requireNamespace("vctrs")

filename <- system.file("extdata", "example.h5ad", package = "anndataR")
file <- rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY", native = TRUE)

test_that("reading encoding works", {
  encoding <- rhdf5_read_h5ad_encoding(file, "obs")
  expect_equal(names(encoding), c("type", "version"))
})

test_that("reading dense matrices works", {
  mat <- rhdf5_read_h5ad_dense_array(file, "layers/dense_counts")
  expect_true(is.matrix(mat))
  expect_type(mat, "integer")
  expect_equal(dim(mat), c(50, 100))

  mat <- rhdf5_read_h5ad_dense_array(file, "layers/dense_X")
  expect_true(is.matrix(mat))
  expect_type(mat, "double")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading sparse matrices works", {
  mat <- rhdf5_read_h5ad_sparse_array(file, "layers/csc_counts", type = "csc")
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(50, 100))

  mat <- rhdf5_read_h5ad_sparse_array(file, "layers/counts", type = "csr")
  expect_s4_class(mat, "dgRMatrix")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading recarrays works", {
  array_list <- rhdf5_read_h5ad_rec_array(
    file,
    "uns/rank_genes_groups/logfoldchanges"
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
  array_1d <- rhdf5_read_h5ad_dense_array(file, "obs/Int")
  expect_equal(array_1d, array(0L:49L))

  array_1d <- rhdf5_read_h5ad_dense_array(file, "obs/Float")
  expect_equal(array_1d, array(rep(42.42, 50)))
})

test_that("reading 1D sparse numeric arrays works", {
  array_1d <- rhdf5_read_h5ad_sparse_array(file, "uns/Sparse1D", type = "csc")
  expect_s4_class(array_1d, "dgCMatrix")
  expect_equal(dim(array_1d), c(1, 6))
})

test_that("reading 1D nullable arrays works", {
  array_1d <- rhdf5_read_h5ad_nullable_integer(file, "obs/IntNA")
  expect_vector(array_1d, ptype = integer(), size = 50)
  expect_true(any(is.na(array_1d)))

  array_1d <- rhdf5_read_h5ad_dense_array(file, "obs/FloatNA")
  expected <- array(rep(42.42, 50))
  expected[1] <- NA
  expect_equal(array_1d, expected)

  array_1d <- rhdf5_read_h5ad_nullable_boolean(file, "obs/Bool")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_false(any(is.na(array_1d)))

  array_1d <- rhdf5_read_h5ad_nullable_boolean(file, "obs/BoolNA")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_true(any(is.na(array_1d)))
})

test_that("reading string scalars works", {
  scalar <- rhdf5_read_h5ad_string_scalar(file, "uns/StringScalar")
  expect_equal(scalar, "A string")
})

test_that("reading numeric scalars works", {
  scalar <- rhdf5_read_h5ad_numeric_scalar(file, "uns/IntScalar")
  expect_equal(scalar, 1)
})

test_that("reading string arrays works", {
  array <- rhdf5_read_h5ad_string_array(file, "uns/String")
  expect_equal(array, array(paste0("String ", 0L:9L)))

  array <- rhdf5_read_h5ad_string_array(file, "uns/String2D")
  expect_true(is.matrix(array))
  expect_type(array, "character")
  expect_equal(dim(array), c(5, 10))
})

test_that("reading mappings works", {
  mapping <- rhdf5_read_h5ad_mapping(file, "uns")
  expect_type(mapping, "list")
  expect_type(names(mapping), "character")
})

test_that("reading dataframes works", {
  df <- rhdf5_read_h5ad_data_frame(file, "obs")
  expect_s3_class(df, "data.frame")
  expect_equal(
    colnames(df),
    c(
      "Float",
      "FloatNA",
      "Int",
      "IntNA",
      "Bool",
      "BoolNA",
      "n_genes_by_counts",
      "log1p_n_genes_by_counts",
      "total_counts",
      "log1p_total_counts",
      "leiden"
    )
  )
})

rhdf5::H5Fclose(file)

test_that("reading H5AD as SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_h5ad(filename, as = "SingleCellExperiment", rhdf5 = TRUE)
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading H5AD as Seurat works", {
  skip_if_not_installed("Seurat")

  seurat <- read_h5ad(filename, as = "Seurat", rhdf5 = TRUE)
  expect_s4_class(seurat, "Seurat")
})

test_that("deprecated to argument in read_h5ad() works", {
  expect_warning(
    mem_ad <- read_h5ad(filename, to = "InMemoryAnnData", rhdf5 = TRUE)
  )
  expect_true(inherits(mem_ad, "InMemoryAnnData"))
})
