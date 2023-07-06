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
  dummy <- dummy_SingleCellExperiment()
  expect_s4_class(dummy, "SingleCellExperiment")
})

test_that("generating dummy Seurat works", {
  dummy <- dummy_Seurat()
  expect_s4_class(dummy, "Seurat")
})
