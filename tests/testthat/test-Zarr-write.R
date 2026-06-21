skip_if_not_installed("Rarr")

store <- tempfile(fileext = ".zarr")
if (dir.exists(store)) {
  unlink(store, recursive = TRUE)
}

create_zarr(store = store)

test_that("Writing Zarr dense arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  expect_silent(write_zarr_element(
    array,
    store,
    "dense_array",
    compression = "none"
  ))
  expect_true(zarr_path_exists(store, "dense_array"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "dense_array"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing Zarr dense 3D arrays works", {
  value <- array(rnorm(60), dim = c(5, 4, 3))

  expect_silent(
    write_zarr_element(
      value,
      store,
      "dense_3d_array"
    )
  )
  expect_true(zarr_path_exists(store, "dense_3d_array"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "dense_3d_array"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing Zarr sparse arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  csc_array <- as(array, "CsparseMatrix")
  expect_silent(write_zarr_element(
    csc_array,
    store,
    "csc_array",
    compression = "none"
  ))
  expect_true(zarr_path_exists(store, "csc_array"))
  expect_true(zarr_path_exists(store, "csc_array/data"))
  expect_true(zarr_path_exists(store, "csc_array/indices"))
  expect_true(zarr_path_exists(store, "csc_array/indptr"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "csc_array"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csc_matrix")

  csr_array <- as(array, "RsparseMatrix")
  expect_silent(write_zarr_element(
    csr_array,
    store,
    "csr_array",
    compression = "none"
  ))
  expect_true(zarr_path_exists(store, "csr_array"))
  expect_true(zarr_path_exists(store, "csr_array/data"))
  expect_true(zarr_path_exists(store, "csr_array/indices"))
  expect_true(zarr_path_exists(store, "csr_array/indptr"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "csr_array"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csr_matrix")
})

test_that("Writing dgeMatrix", {
  value <- matrix(rnorm(20), nrow = 5, ncol = 4) |>
    as("dMatrix") |>
    as("generalMatrix") |>
    as("unpackedMatrix")

  expect_silent(
    write_zarr_element(value, store, "dgematrix")
  )
  expect_true(zarr_path_exists(store, "dgematrix"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "dgematrix"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing Zarr nullable booleans works", {
  nullable <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  nullable[5] <- NA

  expect_silent(write_zarr_element(nullable, store, "nullable_bool"))
  expect_true(zarr_path_exists(store, "nullable_bool"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "nullable_bool"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-boolean")
})

test_that("Writing Zarr nullable integers works", {
  nullable <- as.integer(1:5)
  nullable[5] <- NA

  expect_silent(write_zarr_element(nullable, store, "nullable_int"))
  expect_true(zarr_path_exists(store, "nullable_int"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "nullable_int"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-integer")
})

test_that("Writing Zarr string arrays works", {
  string <- LETTERS[1:5]

  write_zarr_element(string, store, "string_array")
  expect_true(zarr_path_exists(store, "string_array"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "string_array"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")

  string2d <- matrix(LETTERS[1:20], nrow = 5, ncol = 4)

  expect_silent(write_zarr_element(string2d, store, "string_array2D"))
  expect_true(zarr_path_exists(store, "string_array2D"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "string_array2D"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")
})

test_that("Writing Zarr categoricals works", {
  categorical <- factor(LETTERS[1:5])

  expect_no_error(write_zarr_element(categorical, store, "categorical"))
  expect_true(zarr_path_exists(store, "categorical"))
  expect_true(zarr_path_exists(store, "categorical/categories"))
  expect_true(zarr_path_exists(store, "categorical/codes"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "categorical"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "categorical")
})

test_that("Writing Zarr string scalars works", {
  string <- "A"

  expect_silent(write_zarr_element(string, store, "string_scalar"))
  expect_true(zarr_path_exists(store, "string_scalar"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "string_scalar"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string")
})

test_that("Writing Zarr numeric scalars works", {
  number <- 1.0

  expect_silent(write_zarr_element(number, store, "numeric_scalar"))
  expect_true(zarr_path_exists(store, "numeric_scalar"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "numeric_scalar"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "numeric-scalar")
})

test_that("Writing Zarr mappings works", {
  mapping <- list(
    array = matrix(rnorm(20), nrow = 5, ncol = 4),
    sparse = as(matrix(rnorm(20), nrow = 5, ncol = 4), "CsparseMatrix"),
    string = LETTERS[1:5],
    numeric = rnorm(5),
    scalar = 2
  )

  expect_silent(write_zarr_element(
    mapping,
    store,
    "mapping",
    compression = "none"
  ))
  expect_true(zarr_path_exists(store, "mapping"))
  expect_true(zarr_path_exists(store, "mapping/array"))
  expect_true(zarr_path_exists(store, "mapping/sparse"))
  expect_true(zarr_path_exists(store, "mapping/sparse/data"))
  expect_true(zarr_path_exists(store, "mapping/sparse/indices"))
  expect_true(zarr_path_exists(store, "mapping/sparse/indptr"))
  expect_true(zarr_path_exists(store, "mapping/string"))
  expect_true(zarr_path_exists(store, "mapping/numeric"))
  expect_true(zarr_path_exists(store, "mapping/scalar"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "mapping"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dict")
})

test_that("Writing Zarr data frames works", {
  df <- data.frame(
    Letters = letters[1:5],
    Numbers = 1:5
  )

  expect_silent(write_zarr_element(df, store, "dataframe"))
  expect_true(zarr_path_exists(store, "dataframe"))
  expect_true(zarr_path_exists(store, "dataframe/Letters"))
  expect_true(zarr_path_exists(store, "dataframe/Numbers"))
  expect_true(zarr_path_exists(store, "dataframe/_index"))
  attrs <- Rarr::read_zarr_attributes(file.path(store, "dataframe"))
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dataframe")
  expect_true(all(c("_index", "column-order") %in% names(attrs)))
  expect_equal(attrs[["_index"]], "_index")
  expect_identical(
    as.character(attrs[["column-order"]]),
    c("Letters", "Numbers")
  )
})

test_that("writing Zarr from SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")
  store <- tempfile(fileext = ".zarr")
  sce <- generate_dataset(format = "SingleCellExperiment")
  write_zarr(sce, store)
  expect_true(dir.exists(store))
})

test_that("writing Zarr from Seurat works", {
  skip_if_not_installed("SeuratObject")
  store <- tempfile(fileext = ".zarr")
  sce <- generate_dataset(format = "Seurat")
  write_zarr(sce, store)
  expect_true(dir.exists(store))
})

dir_size <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  sum(file.info(files)$size, na.rm = TRUE)
}

test_that("writing compressed files works for Zarr", {
  dummy <- generate_dataset(100, 200)
  non_random_X <- matrix(5, 100, 200) # nolint

  adata <- AnnData(
    X = non_random_X,
    obs = dummy$obs,
    var = dummy$var
  )

  store_none <- tempfile(fileext = ".zarr")
  store_compressed <- tempfile(fileext = ".zarr")

  write_zarr(adata, store_none, compression = "none")

  comp_list <- c("gzip", "blosc", "zstd", "lzma", "bz2", "zlib", "lz4")
  for (comp in comp_list) {
    write_zarr(adata, store_compressed, compression = comp)
    unlink(store_compressed, recursive = TRUE)
    expect_true(dir_size(store_none) > dir_size(store_compressed))
  }
})
