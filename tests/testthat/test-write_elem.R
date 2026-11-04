skip_if_not_installed("rhdf5")

adata <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")

# Every top-level slot and one level below
elem_names <- c(
  "X",
  "obs",
  "var",
  "layers",
  "obsm",
  "varm",
  "obsp",
  "varp",
  "uns",
  unlist(lapply(c("layers", "obsm", "uns"), function(slot) {
    paste0(slot, "/", names(adata[[slot]]))
  }))
)

ref_h5ad <- withr::local_tempfile(fileext = ".h5ad")
write_h5ad(adata, ref_h5ad)

# Write each element to its own group so nested names don't clash with the
# slot they are in
elem_target <- function(name) {
  paste0("elem/", gsub("/", "__", name, fixed = TRUE))
}

# TODO: Scalar NA strings are written to Zarr as "NA", remove this once fixed
drop_scalar_na_string <- function(value) {
  if (is.list(value) && !is.data.frame(value)) {
    value$scalar_character_with_nas <- NULL
  }
  value
}

test_that("write_elem() roundtrips every element for H5AD", {
  path <- withr::local_tempfile(fileext = ".h5ad")

  for (name in elem_names) {
    value <- read_elem(ref_h5ad, name)
    write_elem(path, elem_target(name), value)
    expect_equal(read_elem(path, elem_target(name)), value, info = name)
  }
})

for (zarr_format in c(2L, 3L)) {
  test_that(
    paste0("write_elem() roundtrips every element for Zarr v", zarr_format),
    {
      skip_if_not_installed("Rarr")
      withr::local_options(anndataR.zarr_format = zarr_format)

      path <- withr::local_tempfile(fileext = ".zarr")

      for (name in setdiff(elem_names, "uns/scalar_character_with_nas")) {
        value <- drop_scalar_na_string(read_elem(ref_h5ad, name))
        write_elem(path, elem_target(name), value)
        expect_equal(read_elem(path, elem_target(name)), value, info = name)
      }
      expect_identical(get_zarr_format(path), zarr_format)
    }
  )
}

test_that("write_elem() keeps other elements of an existing H5AD file", {
  path <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, path)
  before <- read_h5ad(path)

  write_elem(path, "uns/new", 1:3)

  after <- read_h5ad(path)
  expect_equal(as.vector(after$uns$new), 1:3)
  expect_equal(after$X, before$X)
  expect_equal(after$obs, before$obs)
  expect_equal(after$layers, before$layers)
  expect_equal(after$uns[names(before$uns)], before$uns)
})

test_that("write_elem() keeps other elements of an existing Zarr store", {
  skip_if_not_installed("Rarr")

  for (zarr_format in c(2L, 3L)) {
    info <- paste0("Zarr v", zarr_format)
    path <- withr::local_tempfile(fileext = ".zarr")
    write_zarr(adata, path, zarr_format = zarr_format)
    before <- read_zarr(path)

    # The store format is used even if the option is different
    withr::with_options(list(anndataR.zarr_format = 5L - zarr_format), {
      write_elem(path, "uns/new", 1:3)
    })

    expect_identical(
      zarr_node_format(path, "uns/new"),
      zarr_format,
      info = info
    )
    expect_equal(as.vector(read_elem(path, "uns/new")), 1:3, info = info)
    after <- read_zarr(path)
    expect_equal(after$X, before$X, info = info)
    expect_equal(after$obs, before$obs, info = info)
    expect_equal(after$layers, before$layers, info = info)
    expect_equal(after$uns[names(before$uns)], before$uns, info = info)
  }
})

test_that("write_elem() creates parent groups", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_elem(h5ad, "/a/b/c", 1:3)
  expect_equal(as.vector(read_elem(h5ad, "a/b/c")), 1:3)

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_elem(zarr, "/a/b/c", 1:3)
  expect_equal(as.vector(read_elem(zarr, "a/b/c")), 1:3)
  expect_true(zarr_path_exists(zarr, "a/b"))
})

test_that("write_elem() only replaces existing elements with overwrite = TRUE", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_elem(h5ad, "x", 1:3)
  expect_error(write_elem(h5ad, "x", 4:6), "already exists")
  expect_equal(as.vector(read_elem(h5ad, "x")), 1:3)
  write_elem(h5ad, "x", letters, overwrite = TRUE)
  expect_equal(as.vector(read_elem(h5ad, "x")), letters)

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_elem(zarr, "x", 1:3)
  expect_error(write_elem(zarr, "x", 4:6), "already exists")
  expect_equal(as.vector(read_elem(zarr, "x")), 1:3)
  write_elem(zarr, "x", letters, overwrite = TRUE)
  expect_equal(as.vector(read_elem(zarr, "x")), letters)
})

test_that("write_elem() follows the anndataR.write_null option", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  withr::with_options(list(anndataR.write_null = TRUE), {
    write_elem(h5ad, "null", NULL)
  })
  expect_null(read_elem(h5ad, "null"))

  withr::with_options(list(anndataR.write_null = FALSE), {
    write_elem(h5ad, "not_written", NULL)
  })
  expect_error(read_elem(h5ad, "not_written"), "not found")
})

test_that("write_elem() passes options in ... to the writers", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  mat <- matrix(as.numeric(1:200), nrow = 10)

  write_elem(h5ad, "gzip", mat, compression = "gzip", chunk_size = NULL)
  expect_equal(read_elem(h5ad, "gzip"), mat)
  expect_error(
    write_elem(h5ad, "bad", mat, compression = "bad"),
    "should be one of"
  )

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_elem(zarr, "gzip", mat, compression = "gzip")
  expect_equal(read_elem(zarr, "gzip"), mat)
})

test_that("write_elem() errors for unsupported arguments in ...", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  expect_error(write_elem(h5ad, "x", 1:3, compresion = "gzip"), "compresion")
  expect_error(
    write_elem(h5ad, "x", 1:3, FALSE, "auto", "gzip"),
    "must be named"
  )
  expect_error(write_elem(h5ad, "x", 1:3, zarr_format = 2), "zarr_format")
  expect_error(
    write_elem(h5ad, "x", 1:3, stop_on_error = FALSE),
    "stop_on_error"
  )

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  expect_error(write_elem(zarr, "x", 1:3, chunk_size = 100), "chunk_size")
  expect_error(write_elem(zarr, "x", 1:3, zarr_format = 2), "zarr_format")
})

test_that("write_elem() errors for the root and unsupported values", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  expect_error(write_elem(h5ad, "/", 1:3), "root group")
  expect_error(write_elem(h5ad, "", 1:3), "root group")
  expect_error(write_elem(h5ad, "x", adata), "not supported")
  expect_error(write_elem(h5ad, "x", function() NULL), "not supported")
})

test_that("write_elem() validates its arguments", {
  dir <- withr::local_tempdir()
  h5ad <- file.path(dir, "test.h5ad")

  expect_error(write_elem(file.path(dir, "test.txt"), "x", 1:3), "format")
  expect_error(write_elem(h5ad, NULL, 1:3), "single string")
  expect_error(write_elem(h5ad, "x", 1:3, overwrite = NA), "overwrite")
  expect_error(write_elem(h5ad, "x", 1:3, format = "loom"), "must be one of")
  expect_false(file.exists(h5ad))
})

test_that("write_elem() returns the path invisibly", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  expect_invisible(write_elem(h5ad, "x", 1:3))
  expect_identical(write_elem(h5ad, "y", 1:3), h5ad)
})
