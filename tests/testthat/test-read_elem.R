skip_if_not_installed("rhdf5")

adata <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")

slots <- c("X", "obs", "var", "layers", "obsm", "varm", "obsp", "varp", "uns")

# Read every top-level slot and one level below, and compare with reading the
# element through the internal reader
expect_read_elem_matches <- function(path, read_internal) {
  for (slot in slots) {
    expect_equal(read_elem(path, slot), read_internal(slot), info = slot)
  }
  for (slot in c("layers", "obsm", "uns")) {
    for (key in names(adata[[slot]])) {
      name <- paste0(slot, "/", key)
      expect_equal(read_elem(path, name), read_internal(name), info = name)
    }
  }
}

test_that("read_elem() reads every element of an H5AD file", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)
  hdf5_file <- HDF5File$new(h5ad)

  expect_read_elem_matches(h5ad, function(name) {
    read_h5ad_element(hdf5_file, name, stop_on_error = TRUE)
  })
})

for (zarr_format in c(2L, 3L)) {
  test_that(
    paste0(
      "read_elem() reads every element of a Zarr v",
      zarr_format,
      " store"
    ),
    {
      skip_if_not_installed("Rarr")

      store <- withr::local_tempfile(fileext = ".zarr")
      write_zarr(adata, store, zarr_format = zarr_format)

      expect_read_elem_matches(store, function(name) {
        read_zarr_element(store, name, stop_on_error = TRUE)
      })
    }
  )
}

test_that("read_elem() returns values matching the AnnData", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)

  expect_equal(read_elem(h5ad, "obs"), adata$obs)
  expect_equal(read_elem(h5ad, "var"), adata$var)
  # Matrices are returned without dimnames, as stored
  expect_equal(read_elem(h5ad, "X"), unname(adata$X))
  expect_null(dimnames(read_elem(h5ad, "X")))
})

test_that("read_elem() accepts a leading slash and an explicit format", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)

  expect_equal(read_elem(h5ad, "/obs"), read_elem(h5ad, "obs"))
  expect_equal(read_elem(h5ad, "obs", format = "h5ad"), read_elem(h5ad, "obs"))
})

test_that("read_elem() with backed = TRUE returns a DelayedArray for H5AD", {
  skip_if_not_installed("HDF5Array")

  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)

  X <- read_elem(h5ad, "X", backed = TRUE)
  expect_s4_class(X, "DelayedArray")
  expect_equal(as.matrix(X), read_elem(h5ad, "X"))
})

test_that("read_elem() errors for missing paths and elements", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)

  expect_error(
    read_elem(file.path(withr::local_tempdir(), "missing.h5ad"), "obs"),
    "does not exist"
  )
  expect_error(read_elem(h5ad, "missing"), "not found")
  expect_error(read_elem(h5ad, "obsm/missing"), "not found")
})

test_that("read_elem() errors for missing Zarr elements and backed Zarr", {
  skip_if_not_installed("Rarr")

  store <- withr::local_tempfile(fileext = ".zarr")
  write_zarr(adata, store)

  expect_error(read_elem(store, "missing"), "not found")
  expect_error(read_elem(store, "X", backed = TRUE), "not yet supported")
})

test_that("read_elem() errors for unsupported encodings", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)

  # Root anndata group
  expect_error(read_elem(h5ad, "/"), "No function for reading")

  # Group without encoding attributes, e.g. a .h5mu modality container
  rhdf5::h5createGroup(h5ad, "no_encoding")
  expect_error(read_elem(h5ad, "no_encoding"), "Encoding attributes not found")

  # raw group
  rhdf5::h5createGroup(h5ad, "raw")
  hdf5_file <- HDF5File$new(h5ad)
  write_h5ad_encoding(hdf5_file, "raw", "raw", "0.1.0")
  expect_error(read_elem(h5ad, "raw"), "No function for reading")
})

test_that("read_elem() validates its arguments", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, h5ad)

  expect_error(read_elem(h5ad, NULL), "single string")
  expect_error(read_elem(h5ad, c("obs", "var")), "single string")
  expect_error(read_elem(h5ad, NA_character_), "single string")
  expect_error(read_elem(h5ad, "obs", backed = NA), "backed")
  expect_error(read_elem(h5ad, "obs", format = "loom"), "must be one of")
  expect_error(read_elem(h5ad, "obs", format = "zarr"), "is \"h5ad\"")
  expect_error(read_elem(h5ad, "X", foo = 1), "foo")
})
