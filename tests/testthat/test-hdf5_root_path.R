test_that("hdf5_root_path is a no-op for the default root", {
  expect_equal(hdf5_root_path("/", "X"), "X")
  expect_equal(hdf5_root_path("", "X"), "X")
  expect_equal(hdf5_root_path("/", "/obs"), "obs")
  expect_equal(hdf5_root_path("/", "/"), "/")
})

test_that("hdf5_root_path prefixes a non-default root", {
  expect_equal(hdf5_root_path("mod/rna", "X"), "mod/rna/X")
  expect_equal(hdf5_root_path("mod/rna", "layers/foo"), "mod/rna/layers/foo")
})

test_that("hdf5_root_path normalizes leading/trailing slashes", {
  expect_equal(hdf5_root_path("/mod/rna/", "X"), "mod/rna/X")
  expect_equal(hdf5_root_path("mod/rna", "/X"), "mod/rna/X")
  expect_equal(hdf5_root_path("/mod/rna/", "/X"), "mod/rna/X")
})

test_that("hdf5_root_path resolves the root itself", {
  expect_equal(hdf5_root_path("mod/rna", "/"), "mod/rna")
  expect_equal(hdf5_root_path("mod/rna", ""), "mod/rna")
})
