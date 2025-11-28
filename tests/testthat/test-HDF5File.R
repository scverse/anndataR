skip_if_not_installed("rhdf5")

file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("creating an HDF5File workds", {
  h5file <- HDF5File$new(file)
  expect_r6_class(h5file, "HDF5File")
})

h5file <- HDF5File$new(file)

test_that("opening an HDF5File works", {
  h5file$open()
  expect_true(h5file$is_open)
})

test_that("opening an open HDF5File works", {
  h5file$open()
  expect_true(h5file$is_open)
})

test_that("opening an HDF5File in read-only mode works", {
  h5file$open(readonly = TRUE)
  expect_true(h5file$is_open)
  expect_true(h5file$is_readonly)
})

test_that("switching HDF5File mode works", {
  h5file$open(readonly = FALSE)
  expect_false(h5file$is_readonly)
})

test_that("closing an HDF5File works", {
  h5file$close()
  expect_false(h5file$is_open)
})

test_that("closing a closed HDF5File works", {
  h5file$close()
  expect_false(h5file$is_open)
})

test_that("opening and deferred closing an HDF5File works", {
  expect_false(h5file$is_open)
  test_func <- function(h5file) {
    h5file$open_and_defer_close()
    expect_true(h5file$is_open)
  }
  test_func(h5file)
  expect_false(h5file$is_open)
})

test_that("nested opening and deferred closing an HDF5File works", {
  expect_false(h5file$is_open)

  # Test when the outer function opens the file
  test_func_outer <- function(h5file) {
    h5file$open_and_defer_close()
    expect_true(h5file$is_open)
    test_func_inner(h5file)
  }

  test_func_inner <- function(h5file) {
    expect_true(h5file$is_open)
  }

  test_func_outer(h5file)
  expect_false(h5file$is_open)

  # Test when the inner function opens the file
  test_func_outer <- function(h5file) {
    expect_false(h5file$is_open)
    test_func_inner(h5file)
    expect_false(h5file$is_open)
  }

  test_func_inner <- function(h5file) {
    h5file$open_and_defer_close()
    expect_true(h5file$is_open)
  }

  test_func_outer(h5file)
  expect_false(h5file$is_open)
})
