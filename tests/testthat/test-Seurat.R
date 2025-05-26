test_that("to_Seurat() fails gracefully", {
  expect_error(to_Seurat(), regexp = "adata.*is missing")
  expect_error(to_Seurat("foo"), regexp = "must be a <AbstractAnnData>")
})

skip_if_not_installed("Seurat")
library(Seurat)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

seu <- ad$as_Seurat()

test_that("as_Seurat retains number of observations and features", {
  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)

  # trackstatus: class=Seurat, feature=test_get_var_names, status=done
  expect_equal(rownames(seu), rownames(ad$var))
  # trackstatus: class=Seurat, feature=test_get_obs_names, status=done
  expect_equal(colnames(seu), rownames(ad$obs))
})

# trackstatus: class=Seurat, feature=test_get_obs, status=done
for (obs_key in colnames(ad$obs)) {
  test_that(paste0("as_Seurat retains obs key: ", obs_key), {
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
  test_that(paste0("as_Seurat retains var key: ", var_key), {
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
  test_that(paste0("as_Seurat retains layer: ", layer_key), {
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
  test_that(paste0("as_Seurat retains uns key: ", uns_key), {
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

test_that("as_Seurat retains pca dimred", {
  msg <- message_if_known(
    backend = "to_Seurat",
    slot = c("obsm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  # trackstatus: class=Seurat, feature=test_get_obsm, status=wip
  expect_true("X_pca" %in% names(seu@reductions))
  expect_equal(
    Embeddings(seu, reduction = "X_pca"),
    ad$obsm[["X_pca"]],
    ignore_attr = TRUE
  )

  # trackstatus: class=Seurat, feature=test_get_varm, status=wip
  expect_equal(
    Loadings(seu, reduction = "X_pca"),
    ad$varm[["PCs"]],
    ignore_attr = TRUE
  )
})

test_that("as_Seurat works with list mappings", {
  expect_no_error(
    ad$as_Seurat(
      object_metadata_mapping = as.list(.to_Seurat_guess_object_metadata(ad)),
      layers_mapping = as.list(.to_Seurat_guess_layers(ad)),
      assay_metadata_mapping = as.list(.to_Seurat_guess_assay_metadata(ad)),
      reduction_mapping = as.list(.to_Seurat_guess_reductions(ad)),
      graph_mapping = as.list(.to_Seurat_guess_graphs(ad)),
      misc_mapping = as.list(.to_Seurat_guess_misc(ad))
    )
  )

  expect_error(
    ad$as_Seurat(
      reduction_mapping = list(numeric = "numeric_matrix")
    )
  )
})

test_that("as_Seurat works with a vector reduction_mapping", {
  expect_no_error(
    ad$as_Seurat(
      reduction_mapping = c(numeric = "numeric_matrix")
    )
  )
})

test_that("as_Seurat works with unnamed mappings", {
  expect_no_error(
    ad$as_Seurat(
      object_metadata_mapping = unname(.to_Seurat_guess_object_metadata(ad)),
      layers_mapping = c(
        na.omit(unname(.to_Seurat_guess_layers(ad))),
        counts = NA
      ),
      assay_metadata_mapping = unname(.to_Seurat_guess_assay_metadata(ad)),
      graph_mapping = unname(.to_Seurat_guess_graphs(ad)),
      misc_mapping = unname(.to_Seurat_guess_misc(ad))
    )
  )
})

test_that("deprecated to_Seurat() works", {
  expect_warning(seu <- ad$to_Seurat())
  expect_s4_class(seu, "Seurat")
})

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

ad <- as_AnnData(obj)

test_that("as_AnnData retains number of observations and features", {
  expect_equal(ad$n_obs(), 200L)
  expect_equal(ad$n_vars(), 100L)

  # trackstatus: class=Seurat, feature=test_set_var_names, status=done
  expect_equal(rownames(ad$var), rownames(obj))
  # trackstatus: class=Seurat, feature=test_set_obs_names, status=done
  expect_equal(rownames(ad$obs), colnames(obj))
})

# trackstatus: class=Seurat, feature=test_set_obs, status=done
for (obs_key in colnames(obj@meta.data)) {
  test_that(paste0("as_AnnData retains obs key: ", obs_key), {
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
  test_that(paste0("as_AnnData retains var key: ", var_key), {
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
  test_that(paste0("as_AnnData retains layer: ", layer_key), {
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

# trackstatus: class=Seurat, feature=test_set_uns, status=done
for (uns_key in names(obj@misc)) {
  test_that(paste0("as_AnnData retains uns key: ", uns_key), {
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

test_that("as_AnnData retains pca", {
  msg <- message_if_known(
    backend = "from_Seurat",
    slot = c("obsm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )
  skip_if(!is.null(msg), message = msg)

  # trackstatus: class=Seurat, feature=test_set_obsm, status=wip
  expect_equal(
    ad$obsm[["pca"]],
    Embeddings(obj, reduction = "pca"),
    ignore_attr = TRUE
  )

  # trackstatus: class=Seurat, feature=test_set_varm, status=wip
  expect_equal(
    ad$varm[["pca"]],
    Loadings(obj, reduction = "pca"),
    ignore_attr = TRUE
  )
})

test_that("as_AnnData retains umap", {
  msg <- message_if_known(
    backend = "from_Seurat",
    slot = c("obsm"),
    dtype = "umap",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  expect_equal(
    ad$obsm[["umap"]],
    Embeddings(obj, reduction = "umap"),
    ignore_attr = TRUE
  )
})

# trackstatus: class=Seurat, feature=test_set_obsp, status=done
test_that("as_AnnData retains connectivities", {
  msg <- message_if_known(
    backend = "from_Seurat",
    slot = c("graphs"),
    dtype = "connectivities",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  expect_equal(
    as.matrix(ad$obsp[["nn"]]),
    as.matrix(obj@graphs[["RNA_nn"]]),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(ad$obsp[["snn"]]),
    as.matrix(obj@graphs[["RNA_snn"]]),
    ignore_attr = TRUE
  )
})

test_that("as_AnnData works with v3 Assays", {
  obj_v3_assay <- obj
  expect_warning(
    obj_v3_assay[["RNA"]] <- as(Seurat::GetAssay(obj, "RNA"), "Assay")
  )

  adata_v3_assay <- from_Seurat(obj_v3_assay)

  expect_identical(
    to_R_matrix(adata_v3_assay$layers$counts),
    SeuratObject::GetAssayData(obj_v3_assay, layer = "counts")
  )
})

test_that("as_AnnData works with list mappings", {
  active_assay <- SeuratObject::DefaultAssay(obj)
  expect_no_error(
    from_Seurat(
      obj,
      layers_mapping = as.list(.from_Seurat_guess_layers(obj, active_assay)),
      obs_mapping = as.list(.from_Seurat_guess_obs(obj, active_assay)),
      var_mapping = as.list(.from_Seurat_guess_var(obj, active_assay)),
      obsm_mapping = as.list(.from_Seurat_guess_obsms(obj, active_assay)),
      varm_mapping = as.list(.from_Seurat_guess_varms(obj, active_assay)),
      obsp_mapping = as.list(.from_Seurat_guess_obsps(obj, active_assay)),
      varp_mapping = as.list(.from_Seurat_guess_varps(obj)),
      uns_mapping = as.list(.from_Seurat_guess_uns(obj))
    )
  )
})

test_that("as_AnnData works with unnamed mappings", {
  active_assay <- SeuratObject::DefaultAssay(obj)
  expect_no_error(
    from_Seurat(
      obj,
      layers_mapping = unname(.from_Seurat_guess_layers(obj, active_assay)),
      obs_mapping = unname(.from_Seurat_guess_obs(obj, active_assay)),
      var_mapping = unname(.from_Seurat_guess_var(obj, active_assay)),
      obsm_mapping = unname(.from_Seurat_guess_obsms(obj, active_assay)),
      varm_mapping = unname(.from_Seurat_guess_varms(obj, active_assay)),
      obsp_mapping = unname(.from_Seurat_guess_obsps(obj, active_assay)),
      varp_mapping = unname(.from_Seurat_guess_varps(obj)),
      uns_mapping = unname(.from_Seurat_guess_uns(obj))
    )
  )
})
