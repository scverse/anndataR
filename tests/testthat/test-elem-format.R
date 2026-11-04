make_zarr_dir <- function(meta_file) {
  store <- withr::local_tempfile(
    fileext = ".zarr",
    .local_envir = parent.frame()
  )
  dir.create(store)
  file.create(file.path(store, meta_file))
  store
}

test_that("detect_elem_format() detects an existing H5AD file", {
  skip_if_not_installed("rhdf5")

  path <- withr::local_tempfile(fileext = ".anything")
  rhdf5::h5createFile(path)

  expect_identical(detect_elem_format(path), "h5ad")
  expect_identical(detect_elem_format(path, "h5ad"), "h5ad")
})

test_that("detect_elem_format() detects existing Zarr v2 and v3 stores", {
  v2 <- make_zarr_dir(".zgroup")
  v3 <- make_zarr_dir("zarr.json")

  expect_identical(detect_elem_format(v2), "zarr")
  expect_identical(detect_elem_format(v3), "zarr")
  expect_identical(detect_elem_format(paste0(v3, "/")), "zarr")
  expect_identical(detect_elem_format(v3, "zarr"), "zarr")
})

test_that("detect_elem_format() uses the extension for new paths", {
  dir <- withr::local_tempdir()

  for (ext in c("h5ad", "H5AD", "h5", "H5", "hdf5", "HDF5", "h5mu", "H5MU")) {
    expect_identical(
      detect_elem_format(file.path(dir, paste0("new.", ext))),
      "h5ad",
      info = ext
    )
  }
  for (ext in c("zarr", "ZARR", "Zarr")) {
    expect_identical(
      detect_elem_format(file.path(dir, paste0("new.", ext))),
      "zarr",
      info = ext
    )
    expect_identical(
      detect_elem_format(file.path(dir, paste0("new.", ext, "/"))),
      "zarr",
      info = ext
    )
  }
})

test_that("detect_elem_format() errors when the format is unclear", {
  dir <- withr::local_tempdir()

  # New path with an unknown extension
  expect_error(
    detect_elem_format(file.path(dir, "new.txt")),
    "Could not detect the format"
  )
  expect_error(
    detect_elem_format(file.path(dir, "new")),
    "Could not detect the format"
  )

  # Existing directory without Zarr metadata
  empty_dir <- file.path(dir, "empty.zarr")
  dir.create(empty_dir)
  expect_error(detect_elem_format(empty_dir), "Could not detect the format")
})

test_that("detect_elem_format() errors for an existing non-HDF5 file", {
  skip_if_not_installed("rhdf5")

  path <- withr::local_tempfile(fileext = ".h5ad")
  writeLines("not hdf5", path)

  expect_error(detect_elem_format(path), "Could not detect the format")
})

test_that("detect_elem_format() returns an explicit format for new paths", {
  dir <- withr::local_tempdir()

  expect_identical(detect_elem_format(file.path(dir, "new"), "h5ad"), "h5ad")
  expect_identical(
    detect_elem_format(file.path(dir, "new.txt"), "zarr"),
    "zarr"
  )
  expect_identical(
    detect_elem_format(file.path(dir, "new.h5ad"), "zarr"),
    "zarr"
  )
})

test_that("detect_elem_format() errors when an explicit format doesn't match", {
  skip_if_not_installed("rhdf5")

  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  rhdf5::h5createFile(h5ad)
  zarr <- make_zarr_dir("zarr.json")

  expect_error(detect_elem_format(h5ad, "zarr"), "is \"h5ad\"")
  expect_error(detect_elem_format(zarr, "h5ad"), "is \"zarr\"")
})

test_that("detect_elem_format() validates its arguments", {
  expect_error(detect_elem_format(NULL), "single non-empty string")
  expect_error(detect_elem_format(NA_character_), "single non-empty string")
  expect_error(detect_elem_format(""), "single non-empty string")
  expect_error(
    detect_elem_format(c("a.h5ad", "b.h5ad")),
    "single non-empty string"
  )
  expect_error(detect_elem_format("a.h5ad", "loom"), "must be one of")
})
