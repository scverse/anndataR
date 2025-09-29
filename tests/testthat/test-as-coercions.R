test_that("as() from SingleCellExperiment to InMemoryAnnData warns and coerces", {
  skip_if_not_installed("SingleCellExperiment")

  counts <- matrix(as.numeric(1:6), nrow = 3, ncol = 2)
  colnames(counts) <- paste0("cell", 1:2)
  rownames(counts) <- paste0("gene", 1:3)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  expect_warning(
    ad <- methods::as(sce, "InMemoryAnnData"),
    "Prefer `as_AnnData\\("
  )
  expect_true(inherits(ad, "InMemoryAnnData"))
  expect_equal(ad$n_obs(), ncol(sce))
  expect_equal(ad$n_vars(), nrow(sce))

  expect_warning(
    sce_roundtrip <- methods::as(ad, "SingleCellExperiment"),
    "as_SingleCellExperiment"
  )
  expect_s4_class(sce_roundtrip, "SingleCellExperiment")
  expect_equal(dim(sce_roundtrip), dim(sce))
})

test_that("as() refuses HDF5AnnData without extra arguments", {
  skip_if_not_installed("SingleCellExperiment")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1))
  )

  expect_error(
    methods::as(sce, "HDF5AnnData"),
    "Use `as_AnnData"
  )
})

test_that("as() converts AnnDataView to InMemoryAnnData with warning", {
  ad <- AnnData(X = matrix(1:4, nrow = 2))
  view <- ad[1, ]
  expect_warning(
    view_materialised <- methods::as(view, "InMemoryAnnData"),
    "adata\\$as_InMemoryAnnData"
  )
  expect_true(inherits(view_materialised, "InMemoryAnnData"))
  expect_equal(view_materialised$n_obs(), 1L)
  expect_equal(view_materialised$n_vars(), ad$n_vars())

  expect_error(
    methods::as(ad, "HDF5AnnData"),
    "adata\\$as_HDF5AnnData"
  )
})

test_that("as() converts AnnData to ReticulateAnnData when available", {
  skip_if_no_anndata_py()

  ad <- AnnData(X = matrix(1:4, nrow = 2))

  expect_warning(
    ad_reticulate <- methods::as(ad, "ReticulateAnnData"),
    "adata\\$as_ReticulateAnnData"
  )
  expect_true(inherits(ad_reticulate, "ReticulateAnnData"))

  expect_warning(
    ad_roundtrip <- methods::as(ad_reticulate, "InMemoryAnnData"),
    "adata\\$as_InMemoryAnnData"
  )
  expect_true(inherits(ad_roundtrip, "InMemoryAnnData"))
  expect_equal(ad_roundtrip$n_obs(), ad$n_obs())
  expect_equal(ad_roundtrip$n_vars(), ad$n_vars())
})

test_that("as() between Seurat and AnnData warns appropriately", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  counts <- Matrix::Matrix(matrix(1:4, nrow = 2), sparse = TRUE)
  colnames(counts) <- paste0("cell", 1:2)
  rownames(counts) <- paste0("gene", 1:2)

  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_warning(
    ad <- methods::as(seurat_obj, "InMemoryAnnData"),
    "Prefer `as_AnnData\\("
  )
  expect_true(inherits(ad, "InMemoryAnnData"))

  expect_warning(
    seurat_roundtrip <- methods::as(ad, "Seurat"),
    "as_Seurat"
  )
  expect_true(inherits(seurat_roundtrip, "Seurat"))

  expect_error(
    methods::as(seurat_obj, "HDF5AnnData"),
    "Use `as_AnnData"
  )
})
