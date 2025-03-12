skip_if_not_installed("Seurat")
library(Seurat)

known_issues <- read_known_issues()

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

active_assay <- obj@assays[[obj@active.assay]]

ad <- from_Seurat(obj)

test_that("from_SCE retains number of observations and features", {
  expect_equal(ad$n_obs(), 200L)
  expect_equal(ad$n_vars(), 100L)

  # trackstatus: class=Seurat, feature=test_set_var_names, status=done
  expect_equal(rownames(ad$var), rownames(obj))
  # trackstatus: class=Seurat, feature=test_set_obs_names, status=done
  expect_equal(rownames(ad$obs), colnames(obj))
})

# trackstatus: class=Seurat, feature=test_set_obs, status=done
for (obs_key in colnames(obj@meta.data)) {
  test_that(paste0("from_Seurat retains obs key: ", obs_key), {
    msg <- message_if_known(
      backend = "from_Seurat",
      slot = c("obs"),
      dtype = obs_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obs_key %in% colnames(ad$obs))
    expect_equal(
      ad$obs[[obs_key]],
      obj@meta.data[[obs_key]],
      info = paste0("obs_key: ", obs_key)
    )
  })
}

# trackstatus: class=Seurat, feature=test_set_var, status=done
for (var_key in colnames(active_assay@meta.data)) {
  test_that(paste0("from_Seurat retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "from_Seurat",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(var_key %in% colnames(ad$var))
    expect_equal(ad$var[[var_key]], active_assay@meta.data[[var_key]])
  })
}

# trackstatus: class=Seurat, feature=test_set_layers, status=done
for (layer_key in names(active_assay@layers)) {
  test_that(paste0("from_Seurat retains layer: ", layer_key), {
    msg <- message_if_known(
      backend = "from_Seurat",
      slot = c("layers"),
      dtype = layer_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(layer_key %in% names(ad$layers))
    expect_true(
      all.equal(
        as.matrix(t(active_assay@layers[[layer_key]])),
        as.matrix(ad$layers[[layer_key]]),
        check.attributes = FALSE
      ),
      info = paste0("layer_key: ", layer_key)
    )
  })
}

# trackstatus: class=Seurat, feature=test_set_uns
for (uns_key in names(obj@misc)) {
  test_that(paste0("from_Seurat retains uns key: ", uns_key), {
    msg <- message_if_known(
      backend = "from_Seurat",
      slot = c("uns"),
      dtype = uns_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(uns_key %in% names(ad$uns))
    expect_equal(ad$uns[[uns_key]], obj@misc[[uns_key]], ignore_attr = TRUE)
  })
}

# trackstatus: class=Seurat, feature=test_set_pca, status=done
test_that("from_Seurat retains pca", {
  msg <- message_if_known(
    backend = "from_Seurat",
    slot = c("obsm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  expect_equal(
    ad$obsm[["X_pca"]],
    Embeddings(obj, reduction = "pca"),
    ignore_attr = TRUE
  )
  expect_equal(
    ad$varm[["pca"]],
    Loadings(obj, reduction = "pca"),
    ignore_attr = TRUE
  )
})

test_that("from_Seurat retains umap", {
  msg <- message_if_known(
    backend = "from_Seurat",
    slot = c("obsm"),
    dtype = "umap",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  expect_equal(
    ad$obsm[["X_umap"]],
    Embeddings(obj, reduction = "umap"),
    ignore_attr = TRUE
  )
})

# trackstatus: class=Seurat, feature=test_set_graphs, status=done
test_that("from_Seurat retains connectivities", {
  msg <- message_if_known(
    backend = "from_Seurat",
    slot = c("graphs"),
    dtype = "connectivities",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  expect_equal(
    as.matrix(ad$obsp[["connectivities"]]),
    as.matrix(obj@graphs[["RNA_nn"]]),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(ad$obsp[["snn"]]),
    as.matrix(obj@graphs[["RNA_snn"]]),
    ignore_attr = TRUE
  )
})
