file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("reading H5AD as SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_h5ad(file, as = "SingleCellExperiment", mode = "r", rhdf5 = TRUE)
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading H5AD as Seurat works", {
  skip_if_not_installed("SeuratObject")
  # TODO: remove this suppression when the to_seurat, from_seurat functions are updated.
  seurat <- suppressWarnings(read_h5ad(file, as = "Seurat", rhdf5 = TRUE))
  expect_s4_class(seurat, "Seurat")
})

test_that("reading H5AD as InMemoryAnnData works", {
  adata <- read_h5ad(file, as = "InMemoryAnnData", mode = "r", rhdf5 = TRUE)
  expect_equal(class(adata), c("InMemoryAnnData", "AbstractAnnData", "R6"))
})
