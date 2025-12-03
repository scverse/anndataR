skip_if_not_installed("rhdf5")
skip_if_not_installed("HDF5Array")

file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("reading H5AD as SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_h5ad(file, as = "SingleCellExperiment", mode = "r")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading H5AD as Seurat works", {
  skip_if_not_installed("SeuratObject")
  # TODO: remove this suppression when the to_seurat, from_seurat functions are updated.
  seurat <- suppressWarnings(read_h5ad(file, as = "Seurat"))
  expect_s4_class(seurat, "Seurat")
})

test_that("reading H5AD as InMemoryAnnData works", {
  adata <- read_h5ad(file, as = "InMemoryAnnData", mode = "r")
  expect_equal(class(adata), c("InMemoryAnnData", "AbstractAnnData", "R6"))
})

test_that("reading H5AD as HDF5AnnData works", {
  adata <- read_h5ad(file, as = "HDF5AnnData", mode = "r")
  expect_equal(class(adata), c("HDF5AnnData", "AbstractAnnData", "R6"))
})

test_that("reading H5AD as backed HDF5AnnData works", {
  adata <- read_h5ad(file, as = "HDF5AnnData", backed = TRUE)

  expect_equal(class(adata), c("HDF5AnnData", "AbstractAnnData", "R6"))

  expect_s4_class(adata$X, "DelayedArray")

  for (layer in adata$layers) {
    expect_s4_class(layer, "DelayedArray")
  }

  for (obsm in adata$obsm) {
    expect_s4_class(obsm, "DelayedArray")
  }

  for (varm in adata$varm) {
    expect_s4_class(varm, "DelayedArray")
  }

  for (obsp in adata$obsp) {
    expect_s4_class(obsp, "DelayedArray")
  }

  for (varp in adata$varp) {
    expect_s4_class(varp, "DelayedArray")
  }
})

test_that("reading H5AD as backed SingleCellExperiment works", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- read_h5ad(file, as = "SingleCellExperiment", backed = TRUE)

  expect_s4_class(sce, "SingleCellExperiment")

  for (assay_name in SummarizedExperiment::assayNames(sce)) {
    expect_s4_class(SummarizedExperiment::assay(sce, assay_name), "DelayedArray")
  }

  for (reduced_dim_name in SingleCellExperiment::reducedDimNames(sce)) {
    reduced_dim <- SingleCellExperiment::reducedDim(sce, reduced_dim_name)

    if (inherits(reduced_dim, "LinearEmbeddingMatrix")) {
      expect_s4_class(SingleCellExperiment::sampleFactors(reduced_dim), "DelayedArray")
      expect_s4_class(SingleCellExperiment::featureLoadings(reduced_dim), "DelayedArray")
    } else {
      expect_s4_class(reduced_dim, "DelayedArray")
    }
  }

  for (row_pair_name in SingleCellExperiment::rowPairNames(sce)) {
    # SingleCellExperiment does not support DelayedArray for rowPairs
    expect_s4_class(SingleCellExperiment::rowPair(sce, row_pair_name), "SelfHits")
  }

  for (col_pair_name in SingleCellExperiment::colPairNames(sce)) {
    # SingleCellExperiment does not support DelayedArray for colPairs
    expect_s4_class(SingleCellExperiment::colPair(sce, col_pair_name), "SelfHits")
  }
})

test_that("reading H5AD as backed Seurat works", {
  skip_if_not_installed("SeuratObject")

  seurat <- read_h5ad(file, as = "Seurat", backed = TRUE)

  expect_s4_class(seurat, "Seurat")

  for (layer in SeuratObject::Layers(seurat)) {
    expect_s4_class(SeuratObject::LayerData(seurat, layer), "DelayedArray")
  }

  for (reduction in SeuratObject::Reductions(seurat)) {
    # Seurat does not support DelayedArray for embeddings and loadings
    expect_true(is.matrix(SeuratObject::Embeddings(seurat, reduction)))
    expect_true(is.matrix(SeuratObject::Loadings(seurat, reduction)))
  }

  for (graph in SeuratObject::Graphs(seurat)) {
    # Seurat does not support DelayedArray for graphs
    expect_true(inherits(seurat[[graph]], "dgCMatrix"))
  }
})
