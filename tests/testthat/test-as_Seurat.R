test_that("as_Seurat() fails gracefully", {
  expect_error(as_Seurat(), regexp = "adata.*is missing")
  expect_error(as_Seurat("foo"), regexp = "must be a <AbstractAnnData>")
})

skip_if_not_installed("Seurat")
library(Seurat)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

layers_mapping <- c(NA, names(ad$layers))
names(layers_mapping) <- c("counts", names(ad$layers))
seu <- ad$as_Seurat(layers_mapping = layers_mapping)

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

    expect_true(obs_key %in% colnames(seu[[]]))
    expect_equal(
      seu@meta.data[[obs_key]],
      ad$obs[[obs_key]],
      info = paste0("obs_key: ", obs_key),
      ignore_attr = TRUE
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
    expect_true(var_key %in% colnames(active_assay[[]]))
    expect_equal(active_assay[[var_key]][[1]], ad$var[[var_key]])
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

    active_assay <- seu[[DefaultAssay(seu)]]
    expect_true(layer_key %in% active_assay[])
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

    expect_true(uns_key %in% names(Misc(seu)))
    expect_equal(Misc(seu)[[uns_key]], ad$uns[[uns_key]], ignore_attr = TRUE)
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
  expect_true("X_pca" %in% Reductions(seu))
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
      x_mapping = "counts",
      object_metadata_mapping = as.list(.as_Seurat_guess_object_metadata(ad)),
      layers_mapping = as.list(.as_Seurat_guess_layers(ad)),
      assay_metadata_mapping = as.list(.as_Seurat_guess_assay_metadata(ad)),
      reduction_mapping = as.list(.as_Seurat_guess_reductions(ad)),
      graph_mapping = as.list(.as_Seurat_guess_graphs(ad)),
      misc_mapping = as.list(.as_Seurat_guess_misc(ad))
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
      x_mapping = "counts",
      reduction_mapping = c(numeric = "numeric_matrix")
    )
  )
})

test_that("as_Seurat works with unnamed mappings", {
  expect_no_error(
    ad$as_Seurat(
      object_metadata_mapping = unname(.as_Seurat_guess_object_metadata(ad)),
      layers_mapping = c(
        na.omit(unname(.as_Seurat_guess_layers(ad))),
        counts = NA
      ),
      assay_metadata_mapping = unname(.as_Seurat_guess_assay_metadata(ad)),
      graph_mapping = unname(.as_Seurat_guess_graphs(ad)),
      misc_mapping = unname(.as_Seurat_guess_misc(ad))
    )
  )
})

test_that("deprecated to_Seurat() works", {
  expect_warning(seu <- ad$to_Seurat(x_mapping = "counts"))
  expect_s4_class(seu, "Seurat")
})
