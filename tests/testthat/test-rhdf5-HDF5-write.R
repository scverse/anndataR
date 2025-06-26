skip_if_not_installed("rhdf5")

requireNamespace("vctrs")

file <- tempfile(pattern = "hdf5_rhdf5_write_", fileext = ".h5ad")
if (file.exists(file)) {
  file.remove(file)
}

file <- rhdf5::H5Fcreate(file)

test_that("Writing H5AD dense arrays works", {
  value <- matrix(rnorm(20), nrow = 5, ncol = 4)

  expect_silent(
    rhdf5_write_h5ad_element(value, file, "dense_array", compression = "none")
  )
  expect_true(rhdf5_hdf5_path_exists(file, "/dense_array"))
  attrs <- rhdf5::h5readAttributes(file, "dense_array", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing H5AD dense 3D arrays works", {
  value <- array(rnorm(60), dim = c(5, 4, 3))

  expect_silent(
    rhdf5_write_h5ad_element(
      value,
      file,
      "dense_3d_array",
      compression = "none"
    )
  )
  expect_true(rhdf5_hdf5_path_exists(file, "/dense_3d_array"))
  attrs <- rhdf5::h5readAttributes(file, "dense_3d_array", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing H5AD sparse arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  csc_array <- as(array, "CsparseMatrix")
  expect_silent(
    rhdf5_write_h5ad_element(csc_array, file, "csc_array", compression = "none")
  )
  expect_true(rhdf5_hdf5_path_exists(file, "/csc_array"))
  expect_true(rhdf5_hdf5_path_exists(file, "/csc_array/data"))
  expect_true(rhdf5_hdf5_path_exists(file, "/csc_array/indices"))
  expect_true(rhdf5_hdf5_path_exists(file, "/csc_array/indptr"))
  attrs <- rhdf5::h5readAttributes(file, "csc_array", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csc_matrix")

  csr_array <- as(array, "RsparseMatrix")
  expect_silent(
    rhdf5_write_h5ad_element(csr_array, file, "csr_array", compression = "none")
  )
  expect_true(rhdf5_hdf5_path_exists(file, "/csr_array"))
  expect_true(rhdf5_hdf5_path_exists(file, "/csr_array/data"))
  expect_true(rhdf5_hdf5_path_exists(file, "/csr_array/indices"))
  expect_true(rhdf5_hdf5_path_exists(file, "/csr_array/indptr"))
  attrs <- rhdf5::h5readAttributes(file, "csr_array", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csr_matrix")
})

test_that("Writing dgeMatrix", {
  value <- matrix(rnorm(20), nrow = 5, ncol = 4) |>
    as("dMatrix") |>
    as("generalMatrix") |>
    as("unpackedMatrix")

  expect_silent(
    rhdf5_write_h5ad_element(value, file, "dgematrix", compression = "none")
  )
  expect_true(rhdf5_hdf5_path_exists(file, "/dgematrix"))
  attrs <- rhdf5::h5readAttributes(file, "dgematrix", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing H5AD nullable booleans works", {
  nullable <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  nullable[5] <- NA

  expect_silent(rhdf5_write_h5ad_element(nullable, file, "nullable_bool"))
  expect_true(rhdf5_hdf5_path_exists(file, "/nullable_bool"))
  attrs <- rhdf5::h5readAttributes(file, "nullable_bool", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-boolean")
})

test_that("Writing H5AD nullable integers works", {
  nullable <- as.integer(1:5)
  nullable[5] <- NA

  expect_silent(rhdf5_write_h5ad_element(nullable, file, "nullable_int"))
  expect_true(rhdf5_hdf5_path_exists(file, "/nullable_int"))
  attrs <- rhdf5::h5readAttributes(file, "nullable_int", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-integer")
})

test_that("Writing H5AD string arrays works", {
  string <- LETTERS[1:5]

  expect_silent(rhdf5_write_h5ad_element(string, file, "string_array"))
  expect_true(rhdf5_hdf5_path_exists(file, "/string_array"))
  attrs <- rhdf5::h5readAttributes(file, "string_array", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")

  string2d <- matrix(LETTERS[1:20], nrow = 5, ncol = 4)

  expect_silent(rhdf5_write_h5ad_element(string2d, file, "string_array2D"))
  expect_true(rhdf5_hdf5_path_exists(file, "/string_array2D"))
  attrs <- rhdf5::h5readAttributes(file, "string_array2D", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")
})

test_that("Writing H5AD categoricals works", {
  categorical <- factor(LETTERS[1:5])

  expect_no_error(rhdf5_write_h5ad_element(categorical, file, "categorical"))
  expect_true(rhdf5_hdf5_path_exists(file, "/categorical"))
  expect_true(rhdf5_hdf5_path_exists(file, "/categorical/categories"))
  expect_true(rhdf5_hdf5_path_exists(file, "/categorical/codes"))
  attrs <- rhdf5::h5readAttributes(file, "categorical", native = FALSE)
  expect_equal(names(attrs), c("encoding-type", "encoding-version", "ordered"))
  expect_equal(attrs[["encoding-type"]], "categorical")
})

test_that("Writing H5AD string scalars works", {
  string <- "A"

  expect_silent(rhdf5_write_h5ad_element(string, file, "string_scalar"))
  expect_true(rhdf5_hdf5_path_exists(file, "/string_scalar"))
  attrs <- rhdf5::h5readAttributes(file, "string_scalar", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string")
})

test_that("Writing H5AD numeric scalars works", {
  number <- 1.0

  expect_silent(rhdf5_write_h5ad_element(number, file, "numeric_scalar"))
  expect_true(rhdf5_hdf5_path_exists(file, "/numeric_scalar"))
  attrs <- rhdf5::h5readAttributes(file, "numeric_scalar", native = FALSE)
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

  expect_silent(
    rhdf5_write_h5ad_element(mapping, file, "mapping", compression = "none")
  )
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/array"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/sparse"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/sparse/data"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/sparse/indices"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/sparse/indptr"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/string"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/numeric"))
  expect_true(rhdf5_hdf5_path_exists(file, "/mapping/scalar"))
  attrs <- rhdf5::h5readAttributes(file, "mapping", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dict")
})

test_that("Writing H5AD data frames works", {
  df <- data.frame(
    Letters = letters[1:5],
    Numbers = 1:5
  )

  expect_silent(rhdf5_write_h5ad_element(df, file, "dataframe"))
  expect_true(rhdf5_hdf5_path_exists(file, "/dataframe"))
  expect_true(rhdf5_hdf5_path_exists(file, "/dataframe/Letters"))
  expect_true(rhdf5_hdf5_path_exists(file, "/dataframe/Numbers"))
  expect_true(rhdf5_hdf5_path_exists(file, "/dataframe/_index"))
  attrs <- rhdf5::h5readAttributes(file, "dataframe", native = FALSE)
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dataframe")
  expect_true(all(c("_index", "column-order") %in% names(attrs)))
  expect_equal(attrs[["_index"]], "_index")
  expect_identical(as.vector(attrs[["column-order"]]), c("Letters", "Numbers"))
})

test_that("writing H5AD from SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  file <- withr::local_file(tempfile(fileext = ".h5ad"))

  sce <- generate_dataset(format = "SingleCellExperiment")
  write_h5ad(sce, file, rhdf5 = TRUE)
  expect_true(file.exists(file))
})

test_that("writing H5AD from Seurat works", {
  skip_if_not_installed("SeuratObject")

  file <- withr::local_file(tempfile(fileext = ".h5ad"))

  seurat <- generate_dataset(format = "Seurat")
  write_h5ad(seurat, file, rhdf5 = TRUE)
  expect_true(file.exists(file))
})

test_that("writing gzip compressed files works", {
  dummy <- generate_dataset(100, 200, example = TRUE)
  non_random_X <- matrix(5, 100, 200) # nolint

  adata <- AnnData(
    X = non_random_X,
    obs = dummy$obs,
    var = dummy$var
  )

  h5ad_file_none <- tempfile(pattern = "hdf5_rhdf5_write_none_", fileext = ".h5ad")
  h5ad_file_gzip <- tempfile(pattern = "hdf5_rhdf5_write_gzip_", fileext = ".h5ad")

  write_h5ad(adata, h5ad_file_none, compression = "none", rhdf5 = TRUE)
  write_h5ad(adata, h5ad_file_gzip, compression = "gzip", rhdf5 = TRUE)

  expect_true(file.info(h5ad_file_none)$size > file.info(h5ad_file_gzip)$size)
})

# TODO: re-enable these tests
# nolint start
# test_that("writing lzf compressed files works", {
#   dummy <- generate_dataset(100, 200, example = FALSE)
#   non_random_X <- matrix(5, 100, 200) # nolint
#
#   adata <- AnnData(
#     X = non_random_X,
#     obs = dummy$obs,
#     var = dummy$var
#   )
#
#   h5ad_file_none <- tempfile(pattern = "hdf5_rhdf5_write_none_", fileext = ".h5ad")
#   h5ad_file_lzf <- tempfile(pattern = "hdf5_rhdf5_write_lzf_", fileext = ".h5ad")
#
#   write_h5ad(adata, h5ad_file_none, compression = "none", rhdf5 = TRUE)
#   write_h5ad(adata, h5ad_file_lzf, compression = "lzf", rhdf5 = TRUE)
#
#   expect_true(file.info(h5ad_file_none)$size > file.info(h5ad_file_lzf)$size)
# })
# nolint end

test_that("write_h5ad() gives consistent hashes", {
  dummy <- generate_dataset(100, 200, example = TRUE, format = "AnnData")

  file <- withr::local_file(tempfile(fileext = ".h5ad"))

  write_h5ad(dummy, file, rhdf5 = TRUE)
  hash1 <- rlang::hash_file(file)

  write_h5ad(dummy, file, rhdf5 = TRUE)
  hash2 <- rlang::hash_file(file)

  expect_identical(hash1, hash2)
})
