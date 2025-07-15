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

active_assay <- obj[[DefaultAssay(obj)]]

ad <- as_AnnData(obj)

test_that("as_AnnData (Seurat) retains number of observations and features", {
  expect_equal(ad$n_obs(), 200L)
  expect_equal(ad$n_vars(), 100L)

  # trackstatus: class=Seurat, feature=test_set_var_names, status=done
  expect_equal(rownames(ad$var), rownames(obj))
  # trackstatus: class=Seurat, feature=test_set_obs_names, status=done
  expect_equal(rownames(ad$obs), colnames(obj))
})

# trackstatus: class=Seurat, feature=test_set_obs, status=done
for (obs_key in colnames(Misc(obj))) {
  test_that(paste0("as_AnnData (Seurat) retains obs key: ", obs_key), {
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
      obj[[obs_key]],
      info = paste0("obs_key: ", obs_key)
    )
  })
}

# trackstatus: class=Seurat, feature=test_set_var, status=done
for (var_key in colnames(active_assay[[]])) {
  test_that(paste0("as_AnnData (Seurat) retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "from_Seurat",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(var_key %in% colnames(ad$var))
    expect_equal(ad$var[[var_key]], active_assay[[var_key]][[1]])
  })
}

# trackstatus: class=Seurat, feature=test_set_layers, status=done
for (layer_key in names(active_assay[])) {
  test_that(paste0("as_AnnData (Seurat) retains layer: ", layer_key), {
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
        as.matrix(t(active_assay[layer_key])),
        as.matrix(ad$layers[[layer_key]]),
        check.attributes = FALSE
      ),
      info = paste0("layer_key: ", layer_key)
    )
  })
}

# trackstatus: class=Seurat, feature=test_set_uns, status=done
for (uns_key in names(Misc(obj))) {
  test_that(paste0("as_AnnData (Seurat) retains uns key: ", uns_key), {
    msg <- message_if_known(
      backend = "from_Seurat",
      slot = c("uns"),
      dtype = uns_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(uns_key %in% names(ad$uns))
    expect_equal(ad$uns[[uns_key]], Misc(obj)[[uns_key]], ignore_attr = TRUE)
  })
}

test_that("as_AnnData (Seurat) retains pca", {
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

test_that("as_AnnData (Seurat) retains umap", {
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
test_that("as_AnnData (Seurat) retains connectivities", {
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
    as.matrix(Graphs(obj, "RNA_nn")),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(ad$obsp[["snn"]]),
    as.matrix(Graphs(obj, "RNA_snn")),
    ignore_attr = TRUE
  )
})

test_that("as_AnnData (Seurat) works with v3 Assays", {
  obj_v3_assay <- obj
  expect_warning(
    obj_v3_assay[["RNA"]] <- as(Seurat::GetAssay(obj, "RNA"), "Assay")
  )

  adata_v3_assay <- as_AnnData(obj_v3_assay)

  expect_identical(
    to_R_matrix(adata_v3_assay$layers$counts),
    SeuratObject::GetAssayData(obj_v3_assay, layer = "counts")
  )
})

test_that("as_AnnData (Seurat) works with list mappings", {
  active_assay <- SeuratObject::DefaultAssay(obj)
  expect_no_error(
    as_AnnData(
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

test_that("as_AnnData (Seurat) works with unnamed mappings", {
  active_assay <- SeuratObject::DefaultAssay(obj)
  expect_no_error(
    as_AnnData(
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

test_that("as_AnnData (Seurat) works with empty mappings", {
  expect_warning(as_AnnData(obj, layers_mapping = NULL), "argument is empty")
  expect_warning(as_AnnData(obj, layers_mapping = c()), "argument is empty")
})
