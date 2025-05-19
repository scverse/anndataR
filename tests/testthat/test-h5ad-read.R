file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("reading H5AD as SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_h5ad(file, to = "SingleCellExperiment", mode = "r")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading H5AD as Seurat works", {
  skip_if_not_installed("SeuratObject")
  # TODO: remove this suppression when the to_seurat, from_seurat functions are updated.
  seurat <- suppressWarnings(read_h5ad(file, to = "Seurat"))
  expect_s4_class(seurat, "Seurat")
})

test_that("reading H5AD as InMemoryAnnData works", {
  adata <- read_h5ad(file, to = "InMemoryAnnData", mode = "r")
  expect_equal(class(adata), c("InMemoryAnnData", "AbstractAnnData", "R6"))
})
