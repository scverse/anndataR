skip_if_not_installed("rhdf5")

h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")
if (file.exists(h5ad_file)) {
  file.remove(h5ad_file)
}

rhdf5::h5createFile(file = h5ad_file)

test_that("Writing H5AD dense arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  expect_silent(write_h5ad_element(array, h5ad_file, "dense_array", compression = "none"))
  expect_true(hdf5_path_exists(h5ad_file, "/dense_array"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "dense_array")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing H5AD sparse arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  csc_array <- as(array, "CsparseMatrix")
  expect_silent(
    write_h5ad_element(csc_array, h5ad_file, "csc_array", compression = "none")
  )
  expect_true(hdf5_path_exists(h5ad_file, "/csc_array"))
  expect_true(hdf5_path_exists(h5ad_file, "/csc_array/data"))
  expect_true(hdf5_path_exists(h5ad_file, "/csc_array/indices"))
  expect_true(hdf5_path_exists(h5ad_file, "/csc_array/indptr"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "csc_array")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csc_matrix")

  csr_array <- as(array, "RsparseMatrix")
  expect_silent(write_h5ad_element(csr_array, h5ad_file, "csr_array", compression = "none"))
  expect_true(hdf5_path_exists(h5ad_file, "/csr_array"))
  expect_true(hdf5_path_exists(h5ad_file, "/csr_array/data"))
  expect_true(hdf5_path_exists(h5ad_file, "/csr_array/indices"))
  expect_true(hdf5_path_exists(h5ad_file, "/csr_array/indptr"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "csr_array")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csr_matrix")
})

test_that("Writing H5AD nullable booleans works", {
  nullable <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  nullable[5] <- NA

  expect_silent(write_h5ad_element(nullable, h5ad_file, "nullable_bool"))
  expect_true(hdf5_path_exists(h5ad_file, "/nullable_bool"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "nullable_bool")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-boolean")
})

test_that("Writing H5AD nullable integers works", {
  nullable <- as.integer(1:5)
  nullable[5] <- NA

  expect_silent(write_h5ad_element(nullable, h5ad_file, "nullable_int"))
  expect_true(hdf5_path_exists(h5ad_file, "/nullable_int"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "nullable_int")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-integer")
})

test_that("Writing H5AD string arrays works", {
  string <- LETTERS[1:5]

  expect_silent(write_h5ad_element(string, h5ad_file, "string_array"))
  expect_true(hdf5_path_exists(h5ad_file, "/string_array"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "string_array")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")

  string2d <- matrix(LETTERS[1:20], nrow = 5, ncol = 4)

  expect_silent(write_h5ad_element(string2d, h5ad_file, "string_array2D"))
  expect_true(hdf5_path_exists(h5ad_file, "/string_array2D"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "string_array2D")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")
})

test_that("Writing H5AD categoricals works", {
  categorical <- factor(LETTERS[1:5])

  expect_no_error(write_h5ad_element(categorical, h5ad_file, "categorical"))
  expect_true(hdf5_path_exists(h5ad_file, "/categorical"))
  expect_true(hdf5_path_exists(h5ad_file, "/categorical/categories"))
  expect_true(hdf5_path_exists(h5ad_file, "/categorical/codes"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "categorical")
  expect_equal(names(attrs), c("encoding-type", "encoding-version", "ordered"))
  expect_equal(attrs[["encoding-type"]], "categorical")
})

test_that("Writing H5AD string scalars works", {
  string <- "A"

  expect_silent(write_h5ad_element(string, h5ad_file, "string_scalar"))
  expect_true(hdf5_path_exists(h5ad_file, "/string_scalar"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "string_scalar")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string")
})

test_that("Writing H5AD numeric scalars works", {
  number <- 1.0

  expect_silent(write_h5ad_element(number, h5ad_file, "numeric_scalar"))
  expect_true(hdf5_path_exists(h5ad_file, "/numeric_scalar"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "numeric_scalar")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "numeric-scalar")
})

test_that("Writing H5AD mappings works", {
  mapping <- list(
    array = matrix(rnorm(20), nrow = 5, ncol = 4),
    sparse = as(matrix(rnorm(20), nrow = 5, ncol = 4), "CsparseMatrix"),
    string = LETTERS[1:5],
    numeric = rnorm(5),
    scalar = 2
  )

  expect_silent(write_h5ad_element(mapping, h5ad_file, "mapping", compression = "none"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/array"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/sparse"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/sparse/data"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/sparse/indices"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/sparse/indptr"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/string"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/numeric"))
  expect_true(hdf5_path_exists(h5ad_file, "/mapping/scalar"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "mapping")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dict")
})

test_that("Writing H5AD data frames works", {
  df <- data.frame(
    Letters = letters[1:5],
    Numbers = 1:5
  )

  expect_silent(write_h5ad_element(df, h5ad_file, "dataframe"))
  expect_true(hdf5_path_exists(h5ad_file, "/dataframe"))
  expect_true(hdf5_path_exists(h5ad_file, "/dataframe/Letters"))
  expect_true(hdf5_path_exists(h5ad_file, "/dataframe/Numbers"))
  expect_true(hdf5_path_exists(h5ad_file, "/dataframe/_index"))
  attrs <- rhdf5::h5readAttributes(h5ad_file, "dataframe")
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dataframe")
  expect_true(all(c("_index", "column-order") %in% names(attrs)))
  expect_equal(attrs[["_index"]], "_index")
  expect_identical(as.vector(attrs[["column-order"]]), c("Letters", "Numbers"))
})

test_that("writing H5AD from SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  file <- withr::local_file("SingleCellExperiment.h5ad")

  sce <- generate_dataset(format = "SingleCellExperiment")
  write_h5ad(sce, file)
  expect_true(file.exists(file))
})

test_that("writing H5AD from Seurat works", {
  skip_if_not_installed("SeuratObject")
  skip("while Seurat converter is failing")

  file <- withr::local_file("Seurat.h5ad")

  seurat <- generate_dataset(format = "Seurat")
  write_h5ad(seurat, file)
  expect_true(file.exists(file))
})

test_that("writing gzip compressed files works", {
  dummy <- generate_dataset(100, 200)
  non_random_X <- matrix(5, 100, 200) # nolint

  adata <- AnnData(
    X = non_random_X,
    obs = dummy$obs[, c(), drop = FALSE],
    var = dummy$var[, c(), drop = FALSE]
  )

  h5ad_file_none <- tempfile(pattern = "hdf5_write_none_", fileext = ".h5ad")
  h5ad_file_gzip <- tempfile(pattern = "hdf5_write_gzip_", fileext = ".h5ad")

  write_h5ad(adata, h5ad_file_none, compression = "none")
  write_h5ad(adata, h5ad_file_gzip, compression = "gzip")

  expect_true(file.info(h5ad_file_none)$size > file.info(h5ad_file_gzip)$size)
})

test_that("writing lzf compressed files works", {
  dummy <- generate_dataset(100, 200)
  non_random_X <- matrix(5, 100, 200) # nolint

  adata <- AnnData(
    X = non_random_X,
    obs = dummy$obs[, c(), drop = FALSE],
    var = dummy$var[, c(), drop = FALSE]
  )

  h5ad_file_none <- tempfile(pattern = "hdf5_write_none_", fileext = ".h5ad")
  h5ad_file_lzf <- tempfile(pattern = "hdf5_write_lzf_", fileext = ".h5ad")

  write_h5ad(adata, h5ad_file_none, compression = "none")
  write_h5ad(adata, h5ad_file_lzf, compression = "lzf")

  expect_true(file.info(h5ad_file_none)$size > file.info(h5ad_file_lzf)$size)
})
