test_that("generating dummy data works", {
  dataset <- generate_dataset()
  expect_type(dataset, "list")
  expect_setequal(
    names(dataset),
    c("X", "obs", "obsp", "obsm", "obs_names", "var", "varp", "varm", "var_names", "layers", "uns")
  )
  expect_identical(dim(dataset$X), c(10L, 20L))
})

test_that("generating dummy SingleCellExperiment works", {
  dummy <- generate_dataset(format = "SingleCellExperiment")
  expect_s4_class(dummy, "SingleCellExperiment")
})

suppressPackageStartupMessages(library(SeuratObject))

test_that("generating dummy Seurat works", {
  dummy <- generate_dataset(format = "Seurat")
  expect_s4_class(dummy, "Seurat")
})
