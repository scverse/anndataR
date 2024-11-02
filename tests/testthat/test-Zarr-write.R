skip_if_not_installed("pizzarr")

store <- pizzarr::MemoryStore$new()

test_that("Writing Zarr dense arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  expect_silent(write_zarr_element(array, store, "dense_array", compression = "none"))
  expect_true(zarr_path_exists(store, "dense_array"))
  g <- pizzarr::zarr_open(store, path = "dense_array")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "array")
})

test_that("Writing Zarr sparse arrays works", {
  array <- matrix(rnorm(20), nrow = 5, ncol = 4)

  csc_array <- as(array, "CsparseMatrix")
  expect_silent(write_zarr_element(csc_array, store, "csc_array", compression = "none"))
  expect_true(zarr_path_exists(store, "csc_array"))
  expect_true(zarr_path_exists(store, "csc_array/data"))
  expect_true(zarr_path_exists(store, "csc_array/indices"))
  expect_true(zarr_path_exists(store, "csc_array/indptr"))
  g <- pizzarr::zarr_open(store, path = "csc_array")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csc_matrix")

  csr_array <- as(array, "RsparseMatrix")
  expect_silent(write_zarr_element(csr_array, store, "csr_array", compression = "none"))
  expect_true(zarr_path_exists(store, "csr_array"))
  expect_true(zarr_path_exists(store, "csr_array/data"))
  expect_true(zarr_path_exists(store, "csr_array/indices"))
  expect_true(zarr_path_exists(store, "csr_array/indptr"))
  g <- pizzarr::zarr_open(store, path = "csr_array")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "csr_matrix")
})

test_that("Writing Zarr nullable booleans works", {
  nullable <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  nullable[5] <- NA

  expect_silent(write_zarr_element(nullable, store, "nullable_bool"))
  expect_true(zarr_path_exists(store, "nullable_bool"))
  g <- pizzarr::zarr_open(store, path = "nullable_bool")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-boolean")
})

test_that("Writing Zarr nullable integers works", {
  nullable <- as.integer(1:5)
  nullable[5] <- NA

  expect_silent(write_zarr_element(nullable, store, "nullable_int"))
  expect_true(zarr_path_exists(store, "nullable_int"))
  g <- pizzarr::zarr_open(store, path = "nullable_int")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "nullable-integer")
})

test_that("Writing Zarr string arrays works", {
  string <- LETTERS[1:5]

  write_zarr_element(string, store, "string_array")
  expect_true(zarr_path_exists(store, "string_array"))
  g <- pizzarr::zarr_open(store, path = "string_array")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")

  string2d <- matrix(LETTERS[1:20], nrow = 5, ncol = 4)

  expect_silent(write_zarr_element(string2d, store, "string_array2D"))
  expect_true(zarr_path_exists(store, "string_array2D"))
  g <- pizzarr::zarr_open(store, path = "string_array2D")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string-array")
})

test_that("Writing Zarr categoricals works", {
  categorical <- factor(LETTERS[1:5])

  expect_no_error(write_zarr_element(categorical, store, "categorical"))
  expect_true(zarr_path_exists(store, "categorical"))
  expect_true(zarr_path_exists(store, "categorical/categories"))
  expect_true(zarr_path_exists(store, "categorical/codes"))
  expect_true(zarr_path_exists(store, "categorical/ordered"))
  g <- pizzarr::zarr_open(store, path = "categorical")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "categorical")
})

test_that("Writing Zarr string scalars works", {
  string <- "A"

  expect_silent(write_zarr_element(string, store, "string_scalar"))
  expect_true(zarr_path_exists(store, "string_scalar"))
  g <- pizzarr::zarr_open(store, path = "string_scalar")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "string")
})

test_that("Writing Zarr numeric scalars works", {
  number <- 1.0

  expect_silent(write_zarr_element(number, store, "numeric_scalar"))
  expect_true(zarr_path_exists(store, "numeric_scalar"))
  g <- pizzarr::zarr_open(store, path = "numeric_scalar")
  attrs <- g$get_attrs()$to_list()
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

  expect_silent(write_zarr_element(mapping, store, "mapping", compression = "none"))
  expect_true(zarr_path_exists(store, "mapping"))
  expect_true(zarr_path_exists(store, "mapping/array"))
  expect_true(zarr_path_exists(store, "mapping/sparse"))
  expect_true(zarr_path_exists(store, "mapping/sparse/data"))
  expect_true(zarr_path_exists(store, "mapping/sparse/indices"))
  expect_true(zarr_path_exists(store, "mapping/sparse/indptr"))
  expect_true(zarr_path_exists(store, "mapping/string"))
  expect_true(zarr_path_exists(store, "mapping/numeric"))
  expect_true(zarr_path_exists(store, "mapping/scalar"))
  g <- pizzarr::zarr_open(store, path = "mapping")
  attrs <- g$get_attrs()$to_list()
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
  g <- pizzarr::zarr_open_group(store, path = "dataframe")
  attrs <- g$get_attrs()$to_list()
  expect_true(all(c("encoding-type", "encoding-version") %in% names(attrs)))
  expect_equal(attrs[["encoding-type"]], "dataframe")
  expect_true(all(c("_index", "column-order") %in% names(attrs)))
  expect_equal(attrs[["_index"]], "_index")
  expect_identical(as.character(attrs[["column-order"]]), c("Letters", "Numbers"))
})

test_that("writing Zarr from SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  store <- pizzarr::MemoryStore$new()

  sce <- generate_dataset(format = "SingleCellExperiment")
  write_zarr(sce, store)
  #expect_true(file.exists(file))
  # TODO: expect things
})

test_that("writing Zarr from Seurat works", {
  skip_if_not_installed("SeuratObject")
  skip("while Seurat converter is failing")

  store <- pizzarr::MemoryStore$new()

  seurat <- generate_dataset(format = "Seurat")
  write_zarr(seurat, store)
  expect_true(file.exists(file))
})

test_that("writing gzip compressed files works for Zarr", {
  dummy <- generate_dataset(100, 200)
  non_random_X <- matrix(5, 100, 200) # nolint

  adata <- AnnData(
    X = non_random_X,
    obs = dummy$obs,
    var = dummy$var
  )

  store_none <- pizzarr::MemoryStore$new()
  store_gzip <- pizzarr::MemoryStore$new()

  write_zarr(adata, store_none, compression = "none")
  write_zarr(adata, store_gzip, compression = "gzip")

  #expect_true(file.info(h5ad_file_none)$size > file.info(h5ad_file_gzip)$size)
  # TODO: expect things
})

test_that("writing lzf compressed files works for Zarr", {
  dummy <- generate_dataset(100, 200)
  non_random_X <- matrix(5, 100, 200) # nolint

  adata <- AnnData(
    X = non_random_X,
    obs = dummy$obs,
    var = dummy$var
  )

  store_none <- pizzarr::MemoryStore$new()
  store_lzf <- pizzarr::MemoryStore$new()

  write_zarr(adata, store_none, compression = "none")
  write_zarr(adata, store_lzf, compression = "lzf")

  #expect_true(file.info(h5ad_file_none)$size > file.info(h5ad_file_lzf)$size)
  # TODO: expect things
})
