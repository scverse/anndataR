test_that("hdf5_group_path is a no-op for the default group", {
  expect_equal(hdf5_group_path("/", "X"), "X")
  expect_equal(hdf5_group_path("", "X"), "X")
  expect_equal(hdf5_group_path("/", "/obs"), "obs")
  expect_equal(hdf5_group_path("/", "/"), "/")
})

test_that("hdf5_group_path prefixes a non-default group", {
  expect_equal(hdf5_group_path("mod/rna", "X"), "mod/rna/X")
  expect_equal(hdf5_group_path("mod/rna", "layers/foo"), "mod/rna/layers/foo")
})

test_that("hdf5_group_path normalizes leading/trailing slashes", {
  expect_equal(hdf5_group_path("/mod/rna/", "X"), "mod/rna/X")
  expect_equal(hdf5_group_path("mod/rna", "/X"), "mod/rna/X")
  expect_equal(hdf5_group_path("/mod/rna/", "/X"), "mod/rna/X")
})

test_that("hdf5_group_path resolves the group itself", {
  expect_equal(hdf5_group_path("mod/rna", "/"), "mod/rna")
  expect_equal(hdf5_group_path("mod/rna", ""), "mod/rna")
})
