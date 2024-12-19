dummy <- generate_dataset(10L, 20L)

test_that("to_Seurat with inmemoryanndata", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var
  )
  # running to_seurat when ad0$X is null probably doesn't make any sense
  ad0 <- AnnData(
    X = matrix(1L, nrow = 5, ncol = 10),
    obs = data.frame(row.names = letters[1:5]),
    var = data.frame(row.names = LETTERS[1:10])
  )

  # TODO: remove suppressWarnings when to_Seurat gets updated
  seu <- suppressWarnings(ad$to_Seurat())
  seu0 <- suppressWarnings(ad0$to_Seurat())

  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)

  # trackstatus: class=Seurat, feature=test_get_var_names, status=done
  expect_equal(rownames(seu), rownames(dummy$var))
  expect_equal(rownames(seu0), LETTERS[1:10])

  # trackstatus: class=Seurat, feature=test_get_obs_names, status=done
  expect_equal(colnames(seu), rownames(dummy$obs))
  expect_equal(colnames(seu0), letters[1:5])

  # check whether all obs keys are found in the seu metadata
  # seurat will have extra metadata columns
  # trackstatus: class=Seurat, feature=test_get_obs, status=done
  for (obs_key in colnames(dummy$obs)) {
    expect_true(obs_key %in% colnames(seu@meta.data))
    expect_equal(seu@meta.data[[obs_key]], dummy$obs[[obs_key]])
  }
  seu0@meta.data

  # check whether all var keys are found in the seu assay metadata
  # trackstatus: class=Seurat, feature=test_get_var, status=done
  active_assay <- seu@assays[[seu@active.assay]]
  for (var_key in colnames(dummy$var)) {
    expect_true(var_key %in% colnames(active_assay@meta.data))
    expect_equal(active_assay@meta.data[[var_key]], dummy$var[[var_key]])
  }
})

test_that("to_Seurat() fails gracefully", {
  expect_error(to_Seurat(), regexp = "adata.*is missing")
  expect_error(to_Seurat("foo"), regexp = "AbstractAnnData.*not TRUE")
})

test_that("from_Seurat() works", {
  skip_if_not_installed("Seurat")

  library(Seurat)

  suppressWarnings({
    counts <- matrix(rbinom(20000, 1000, .001), nrow = 100)
    obj <- CreateSeuratObject(counts = counts)
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, npcs = 10L)
    obj <- FindNeighbors(obj)
    obj <- RunUMAP(obj, dims = 1:10)
  })

  ad <- from_Seurat(obj)

  # trackstatus: class=Seurat, feature=test_set_X, status=done
  # trackstatus: class=Seurat, feature=test_set_layers, status=done
  expect_equal(ad$n_obs(), 200L)
  expect_equal(ad$n_vars(), 100L)

  expect_equal(ad$layers_keys(), c("counts", "data", "scale.data"))

  # trackstatus: class=Seurat, feature=test_set_obsm, status=done
  expect_equal(ad$obsm_keys(), c("X_pca", "X_umap"))

  # trackstatus: class=Seurat, feature=test_set_varm, status=done
  expect_equal(ad$varm_keys(), "PCs")
})
