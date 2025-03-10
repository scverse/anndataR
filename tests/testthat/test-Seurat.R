test_that("to_Seurat() fails gracefully", {
  expect_error(to_Seurat(), regexp = "adata.*is missing")
  expect_error(to_Seurat("foo"), regexp = "must be a <AbstractAnnData>")
})

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_var = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

skip_if_not_installed("Seurat")
library(Seurat)

##################
# TEST TO_SEURAT #
##################

seu <- ad$to_Seurat()

test_that("to_Seurat retains number of observations and features", {
  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)

  # trackstatus: class=Seurat, feature=test_get_var_names, status=done
  expect_equal(rownames(seu), rownames(ad$var))
  # tracksstatus: class=Seurat, feature=test_get_obs_names, status=done
  expect_equal(colnames(seu), rownames(ad$obs))
})

# trackstatus: class=Seurat, feature=test_get_obs, status=done
for (obs_key in colnames(ad$obs)) {
  test_that(paste0("to_Seurat retains obs key: ", obs_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("obs"),
      dtype = obs_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obs_key %in% colnames(seu@meta.data))
    expect_equal(
      seu@meta.data[[obs_key]],
      ad$obs[[obs_key]],
      info = paste0("obs_key: ", obs_key)
    )
  })
}

# trackstatus: class=Seurat, feature=test_get_var, status=done
for (var_key in colnames(ad$var)) {
  test_that(paste0("to_Seurat retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    active_assay <- seu[[DefaultAssay(seu)]]
    expect_true(var_key %in% colnames(active_assay@meta.data))
    expect_equal(active_assay@meta.data[[var_key]], ad$var[[var_key]])
  })
}


# trackstatus: class=Seurat, feature=test_get_layers, status=done
for (layer_key in names(ad$layers)) {
  test_that(paste0("to_Seurat retains layer: ", layer_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("layers"),
      dtype = layer_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    active_assay <- seu@assays[[seu@active.assay]]

    expect_true(layer_key %in% names(active_assay@layers))
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

# trackstatus: class=Seurat, feature=test_get_uns, status=done
for (uns_key in names(ad$uns)) {
  test_that(paste0("to_Seurat retains uns key: ", uns_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("uns"),
      dtype = uns_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(uns_key %in% names(seu@misc))
    expect_equal(seu@misc[[uns_key]], ad$uns[[uns_key]], ignore_attr = TRUE)
  })
}

# trackstatus: class=Seurat, feature=test_get_pca, status=done
test_that("to_Seurat retains pca dimred", {
  msg <- message_if_known(
    backend = "to_Seurat",
    slot = c("obsm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  expect_true("pca" %in% names(seu@reductions))
  expect_equal(
    Embeddings(seu, reduction = "pca"),
    ad$obsm[["X_pca"]],
    ignore_attr = TRUE
  )
  expect_equal(
    Loadings(seu, reduction = "pca"),
    ad$varm[["PCs"]],
    ignore_attr = TRUE
  )
})


####################
# TEST FROM_SEURAT #
####################

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
