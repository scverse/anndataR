test_that("generating dummy data works", {
  dummy <- dummy_data()
  expect_type(dummy, "list")
  expect_identical(
    names(dummy),
    c("X", "obs", "obs_names", "var", "var_names", "layers")
  )
  expect_identical(dim(dummy$X), c(10L, 20L))
})

test_that("generating dummy SingleCellExperiment works", {
  dummy <- dummy_data(output = "SingleCellExperiment")
  expect_s4_class(dummy, "SingleCellExperiment")
})

test_that("generating dummy Seurat works", {
  dummy <- dummy_data(output = "Seurat")
  expect_s4_class(dummy, "Seurat")
})

test_that("generating dummy InMemoryAnnData works", {
  dummy <- dummy_data(output = "InMemoryAnnData")
  expect_true(methods::is(dummy,"InMemoryAnnData"))
})

test_that("generating dummy HDF5AnnData works", {
  dummy <- dummy_data(output = "HDF5AnnData")
  expect_true(methods::is(dummy,"HDF5AnnData"))
})
