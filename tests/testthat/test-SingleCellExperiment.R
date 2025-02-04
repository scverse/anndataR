known_issues <- read_known_issues()

test_that("to_SingleCellExperiment with inmemoryanndata", { # nolint
  library(SingleCellExperiment)

  ad <- generate_dataset(n_obs = 10L, n_var = 20L, format = "AnnData")

  ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
  ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

  sce <- ad$to_SingleCellExperiment()

  expect_equal(nrow(sce), 20)
  expect_equal(ncol(sce), 10)

  # trackstatus: class=SingleCellExperiment, feature=test_get_var_names, status=done
  expect_equal(rownames(sce), rownames(ad$var))

  # trackstatus: class=SingleCellExperiment, feature=test_get_obs_names, status=done
  expect_equal(colnames(sce), rownames(ad$obs))

  # check whether all obs keys are found in the sce metadata
  # trackstatus: class=SingleCellExperiment, feature=test_get_obs, status=done
  for (obs_key in colnames(ad$obs)) {
    expect_true(obs_key %in% colnames(colData(sce)))
    expect_equal(colData(sce)[[obs_key]], ad$obs[[obs_key]], info = paste0("obs_key: ", obs_key))
  }

  # check whether all var keys are found in the sce assay metadata
  # trackstatus: class=SingleCellExperiment, feature=test_get_var, status=done
  for (var_key in colnames(ad$var)) {
    expect_true(var_key %in% colnames(rowData(sce)))
    expect_equal(rowData(sce)[[var_key]], ad$var[[var_key]], info = paste0("var_key: ", var_key))
  }

  # check whether layers are found in the sce assays
  # trackstatus: class=SingleCellExperiment, feature=test_get_layers, status=done
  for (layer_key in names(ad$layers)) {
    expect_true(layer_key %in% names(assays(sce)))
    expect_true(
      all.equal(assay(sce, layer_key), t(ad$layers[[layer_key]]), check.attributes = FALSE),
      info = paste0("layer_key: ", layer_key)
    )
  }

  # check whether all obsp keys are found in the colPairs
  # trackstatus: class=SingleCellExperiment, feature=test_get_obsp, status=done
  for (obsp_key in names(ad$obsp)) {
    expect_true(obsp_key %in% names(colPairs(sce)))

    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("obsp"),
      dtype = obsp_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    sce_matrix <- as.matrix(colPair(sce, obsp_key, asSparse = TRUE))
    ad_matrix <- as.matrix(ad$obsp[[obsp_key]])

    expect_equal(sce_matrix, ad_matrix, info = paste0("obsp_key: ", obsp_key))
  }

  # check whether all varp keys are found in the rowPairs
  # trackstatus: class=SingleCellExperiment, feature=test_get_varp, status=done
  for (varp_key in names(ad$varp)) {
    expect_true(varp_key %in% names(rowPairs(sce)))

    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("obsp"),
      dtype = varp_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    sce_matrix <- as.matrix(rowPair(sce, varp_key, asSparse = TRUE))
    ad_matrix <- as.matrix(ad$varp[[varp_key]])

    expect_equal(sce_matrix, ad_matrix, info = paste0("varp_key: ", varp_key))

  }

  # check reduction
  # check if ad$obsm[["X_pca"]] is found in the sce SingleCellExperiment::sampleFactors(reducedDims$PCA)
  # and check if ad$varm[["PCs"]] is found in the sce SingleCellExperiment::featureLoadings(reducedDims$PCA)

  expect_true("pca" %in% names(reducedDims(sce)))
  expect_equal(
    reducedDims(sce)$pca,
    ad$obsm[["X_pca"]],
    info = "reducedDims(pca)"
  )
  expect_equal(
    featureLoadings(sce)$pca,
    ad$varm[["PC"]],
    info = "featureLoadings(pca)"
  )

  # TODO: uns keys?

})


# TODO gracefully failing

test_that("from_SingleCellExperiment works", {
  # reverse from to_SCE --> assays to layers, colData to obs, rowData to var, colPairs to obsp, rowPairs to varp

  skip_if_not_installed("SingleCellExperiment")
  library(SingleCellExperiment)

  ad1 <- generate_dataset(n_obs = 10L, n_var = 20L, format = "AnnData")

  ad1$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
  ad1$varm[["PCs"]] <- matrix(1:100, 20, 5)

  sce <- ad1$to_SingleCellExperiment()

  ad <- from_SingleCellExperiment(sce)

  # trackstatus: class=SingleCellExperiment, feature=test_set_var_names, status=done
  expect_equal(rownames(sce), rownames(ad$var))

  # trackstatus: class=SingleCellExperiment, feature=test_set_obs_names, status=done
  expect_equal(colnames(sce), rownames(ad$obs))

  # check whether all sce coldata keys are found in the obs
  # trackstatus: class=SingleCellExperiment, feature=test_set_obs, status=done
  for (obs_key in colnames(colData(sce))) {
    expect_true(obs_key %in% colnames(ad$obs))
    expect_equal(ad$obs[[obs_key]], colData(sce)[[obs_key]], info = paste0("obs_key: ", obs_key))
  }

  # check whether all sce rowdata keys are found in the var
  # trackstatus: class=SingleCellExperiment, feature=test_set_var, status=done
  for (var_key in colnames(rowData(sce))) {
    expect_true(var_key %in% colnames(ad$var))
    expect_equal(ad$var[[var_key]], rowData(sce)[[var_key]], info = paste0("var_key: ", var_key))
  }

  # check whether assays are found in the layers
  # not all true, sometimes they change matrix type --> is this a known issue or not?
  for (layer_key in names(assays(sce))) {
    expect_true(layer_key %in% names(ad$layers), info = paste0("layer_key: ", layer_key))
    expect_true(
      all.equal(as.matrix(ad$layers[[layer_key]]), as.matrix(t(assay(sce, layer_key))), check.attributes = FALSE),
      info = paste0("layer_key: ", layer_key)
    )
  }

  # check whether all colPairs keys are found in the obsp
  for (obsp_key in names(colPairs(sce))) {
    expect_true(obsp_key %in% names(ad$obsp))
    expect_equal(
      ad$obsp[[obsp_key]], colPairs(sce, asSparse = TRUE)[[obsp_key]], check.attributes = FALSE,
      info = paste0("obsp_key: ", obsp_key)
    )
  }

  # check whether all rowPairs keys are found in the varp
  for (varp_key in names(rowPairs(sce))) {
    expect_true(varp_key %in% names(ad$varp))
    expect_equal(ad$varp[[varp_key]], rowPairs(sce, asSparse = TRUE)[[varp_key]], info = paste0("varp_key: ", varp_key))
  }

  # check whether all metadata keys are found in the uns
  for (uns_key in names(metadata(sce))) {
    expect_true(uns_key %in% names(ad$uns))
    expect_equal(ad$uns[[uns_key]], metadata(sce)[[uns_key]], info = paste0("uns_key: ", uns_key))
  }

  # TODO: obsm keys? varm keys? --> test a reduction

})
