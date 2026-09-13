skip_if_not_installed("rhdf5")

test_that("hdf5_ensure_group_path creates all missing intermediate groups", {
  file <- withr::local_tempfile(fileext = ".h5")
  rhdf5::h5createFile(file)
  hdf5_file <- HDF5File$new(file)

  hdf5_ensure_group_path(hdf5_file, "mod/rna")

  expect_true(hdf5_path_exists(hdf5_file, "mod"))
  expect_true(hdf5_path_exists(hdf5_file, "mod/rna"))
})

test_that("hdf5_ensure_group_path is a no-op when the group already exists", {
  file <- withr::local_tempfile(fileext = ".h5")
  rhdf5::h5createFile(file)
  hdf5_file <- HDF5File$new(file)

  hdf5_create_group(hdf5_file, "mod")
  hdf5_create_group(hdf5_file, "mod/rna")

  expect_no_error(hdf5_ensure_group_path(hdf5_file, "mod/rna"))
})

test_that("hdf5_ensure_group_path only creates missing segments", {
  file <- withr::local_tempfile(fileext = ".h5")
  rhdf5::h5createFile(file)
  hdf5_file <- HDF5File$new(file)

  hdf5_create_group(hdf5_file, "mod")

  hdf5_ensure_group_path(hdf5_file, "mod/rna")

  expect_true(hdf5_path_exists(hdf5_file, "mod/rna"))
})

test_that("hdf5_ensure_group_path handles a single-segment path", {
  file <- withr::local_tempfile(fileext = ".h5")
  rhdf5::h5createFile(file)
  hdf5_file <- HDF5File$new(file)

  hdf5_ensure_group_path(hdf5_file, "mod")

  expect_true(hdf5_path_exists(hdf5_file, "mod"))
})
