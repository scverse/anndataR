skip_if_not_installed("SingleCellExperiment")
library(SingleCellExperiment)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

sce <- ad$to_SingleCellExperiment()

test_that("to_SCE retains nr of observations and features", {
  expect_equal(nrow(sce), 20)
  expect_equal(ncol(sce), 10)

  # trackstatus: class=SingleCellExperiment, feature=test_get_var_names, status=done
  expect_equal(rownames(sce), rownames(ad$var))
  # trackstatus: class=SingleCellExperiment, feature=test_get_obs_names, status=done
  expect_equal(colnames(sce), rownames(ad$obs))
})

# trackstatus: class=SingleCellExperiment, feature=test_get_obs, status=done
for (obs_key in colnames(ad$obs)) {
  test_that(paste0("to_SCE retains obs key: ", obs_key), {
    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("obs"),
      dtype = obs_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obs_key %in% colnames(colData(sce)))
    expect_equal(
      colData(sce)[[obs_key]],
      ad$obs[[obs_key]],
      info = paste0("obs_key: ", obs_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_get_var, status=done
for (var_key in colnames(ad$var)) {
  test_that(paste0("to_SCE retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(var_key %in% colnames(rowData(sce)))
    expect_equal(
      rowData(sce)[[var_key]],
      ad$var[[var_key]],
      info = paste0("var_key: ", var_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_get_layers, status=done
for (layer_key in names(ad$layers)) {
  test_that(paste0("to_SCE retains layer: ", layer_key), {
    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("layers"),
      dtype = layer_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(layer_key %in% names(assays(sce)))
    expect_true(
      all.equal(
        as.matrix(t(assay(sce, layer_key))),
        as.matrix(ad$layers[[layer_key]]),
        check.attributes = FALSE
      ),
      info = paste0("layer_key: ", layer_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_get_obsp, status=done
for (obsp_key in names(ad$obsp)) {
  test_that(paste0("to_SCE retains obsp key: ", obsp_key), {
    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("obsp"),
      dtype = obsp_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obsp_key %in% names(colPairs(sce)))

    sce_matrix <- as.matrix(colPair(sce, obsp_key, asSparse = TRUE))
    ad_matrix <- as.matrix(ad$obsp[[obsp_key]])

    expect_equal(sce_matrix, ad_matrix, info = paste0("obsp_key: ", obsp_key))
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_get_varp, status=done
for (varp_key in names(ad$varp)) {
  test_that(paste0("to_SCE retains varp key: ", varp_key), {
    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("obsp"),
      dtype = varp_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(varp_key %in% names(rowPairs(sce)))

    sce_matrix <- as.matrix(rowPair(sce, varp_key, asSparse = TRUE))
    ad_matrix <- as.matrix(ad$varp[[varp_key]])

    expect_equal(sce_matrix, ad_matrix, info = paste0("varp_key: ", varp_key))
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_get_uns, status=done
for (uns_key in names(ad$uns)) {
  test_that(paste0("to_SCE retains uns key: ", uns_key), {
    msg <- message_if_known(
      backend = "to_SCE",
      slot = c("uns"),
      dtype = uns_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(uns_key %in% names(metadata(sce)))
    expect_equal(
      metadata(sce)[[uns_key]],
      ad$uns[[uns_key]],
      info = paste0("uns_key: ", uns_key)
    )
  })
}

test_that("to_SCE retains pca dimred", {
  msg <- message_if_known(
    backend = "to_SCE",
    slot = c("obsm", "varm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )
  skip_if(!is.null(msg), message = msg)

  # trackstatus: class=SingleCellExperiment, feature=test_get_obsm, status=wip
  expect_true("pca" %in% names(reducedDims(sce)))
  expect_equal(
    sampleFactors(reducedDims(sce)$pca),
    ad$obsm[["X_pca"]],
    ignore.attributes = TRUE
  )

  # trackstatus: class=SingleCellExperiment, feature=test_get_varm, status=wip
  expect_equal(
    featureLoadings(reducedDims(sce)$pca),
    ad$varm[["PCs"]]
  )
})

test_that("to_SCE works with list mappings", {
  expect_no_error(
    ad$to_SingleCellExperiment(
      assays_mapping = as.list(.to_SCE_guess_assays(ad)),
      colData_mapping = as.list(.to_SCE_guess_all(ad, "obs")),
      rowData_mapping = as.list(.to_SCE_guess_all(ad, "var")),
      reduction_mapping = as.list(.to_SCE_guess_reduction(ad)),
      colPairs_mapping = as.list(.to_SCE_guess_all(ad, "obsp")),
      rowPairs_mapping = as.list(.to_SCE_guess_all(ad, "varp")),
      metadata_mapping = as.list(.to_SCE_guess_all(ad, "uns"))
    )
  )

  expect_error(
    ad$to_SingleCellExperiment(
      reduction_mapping = list(numeric = "numeric_matrix")
    )
  )
})

test_that("to_SCE works with a vector reduction_mapping", {
  expect_no_error(
    ad$to_SingleCellExperiment(
      reduction_mapping = c(numeric = "numeric_matrix")
    )
  )
})

test_that("to_SCE works with unnamed mappings", {
  expect_no_error(
    ad$to_SingleCellExperiment(
      assays_mapping = unname(.to_SCE_guess_assays(ad)),
      colData_mapping = unname(.to_SCE_guess_all(ad, "obs")),
      rowData_mapping = unname(.to_SCE_guess_all(ad, "var")),
      colPairs_mapping = unname(.to_SCE_guess_all(ad, "obsp")),
      rowPairs_mapping = unname(.to_SCE_guess_all(ad, "varp")),
      metadata_mapping = unname(.to_SCE_guess_all(ad, "uns"))
    )
  )
})

# TODO gracefully failing

skip_if_not_installed("SingleCellExperiment")
library(SingleCellExperiment)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

# TODO: Build an SCE rather than converting
sce <- ad$to_SingleCellExperiment()
ad <- as_AnnData(sce)

test_that("from_SCE retains observations and features", {
  expect_equal(nrow(sce), 20)
  expect_equal(ncol(sce), 10)

  # trackstatus: class=SingleCellExperiment, feature=test_set_obs_names, status=done
  expect_equal(colnames(sce), rownames(ad$obs))
  # trackstatus: class=SingleCellExperiment, feature=test_set_var_names, status=done
  expect_equal(rownames(sce), rownames(ad$var))
})

# trackstatus: class=SingleCellExperiment, feature=test_set_obs, status=done
for (obs_key in colnames(colData(sce))) {
  test_that(paste0("from_SCE retains obs key: ", obs_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("obs"),
      dtype = obs_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obs_key %in% colnames(ad$obs))
    expect_equal(
      ad$obs[[obs_key]],
      colData(sce)[[obs_key]],
      info = paste0("obs_key: ", obs_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_set_var, status=done
for (var_key in colnames(rowData(sce))) {
  test_that(paste0("from_SCE retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(var_key %in% colnames(ad$var))
    expect_equal(
      ad$var[[var_key]],
      rowData(sce)[[var_key]],
      info = paste0("var_key: ", var_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_set_layers, status=done
for (layer_key in names(assays(sce))) {
  test_that(paste0("from_SCE retains layer: ", layer_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("layers"),
      dtype = layer_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(layer_key %in% names(ad$layers))
    expect_true(
      all.equal(
        as.matrix(ad$layers[[layer_key]]),
        as.matrix(t(assay(sce, layer_key))),
        ignore_attr = TRUE
      ),
      info = paste0("layer_key: ", layer_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_set_obsp, status=done
for (obsp_key in names(colPairs(sce))) {
  test_that(paste0("from_SCE retains obsp key: ", obsp_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("obsp"),
      dtype = obsp_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obsp_key %in% names(ad$obsp))
    expect_equal(
      ad$obsp[[obsp_key]],
      colPairs(sce, asSparse = TRUE)[[obsp_key]],
      ignore_attr = TRUE,
      info = paste0("obsp_key: ", obsp_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_set_varp, status=done
for (varp_key in names(rowPairs(sce))) {
  test_that(paste0("from_SCE retains varp key: ", varp_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("varp"),
      dtype = varp_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(varp_key %in% names(ad$varp))
    expect_equal(
      ad$varp[[varp_key]],
      rowPairs(sce, asSparse = TRUE)[[varp_key]],
      info = paste0("varp_key: ", varp_key)
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=test_set_uns, status=done
for (uns_key in names(metadata(sce))) {
  test_that(paste0("from_SCE retains uns key: ", uns_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("uns"),
      dtype = uns_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(uns_key %in% names(ad$uns))
    expect_equal(
      ad$uns[[uns_key]],
      metadata(sce)[[uns_key]],
      info = paste0("uns_key: ", uns_key)
    )
  })
}

test_that("from_SCE retains pca dimred", {
  msg <- message_if_known(
    backend = "from_SCE",
    slot = c("obsm", "varm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )
  skip_if(!is.null(msg), message = msg)

  # trackstatus: class=SingleCellExperiment, feature=test_set_obsm, status=wip
  expect_true("X_pca" %in% names(ad$obsm))
  expect_equal(
    sampleFactors(reducedDims(sce)$X_pca),
    ad$obsm[["X_pca"]]
  )

  # trackstatus: class=SingleCellExperiment, feature=test_set_varm, status=wip
  expect_true("X_pca" %in% names(ad$varm))
  expect_equal(
    featureLoadings(reducedDims(sce)$X_pca),
    ad$varm[["X_pca"]]
  )
})

test_that("from_SCE works with list mappings", {
  expect_no_error(
    as_AnnData(
      sce,
      layers_mapping = as.list(.from_SCE_guess_layers(sce, NULL)),
      obs_mapping = as.list(.from_SCE_guess_all(sce, colData)),
      var_mapping = as.list(.from_SCE_guess_all(sce, rowData)),
      obsm_mapping = as.list(.from_SCE_guess_obsm(sce)),
      varm_mapping = as.list(.from_SCE_guess_varm(sce)),
      obsp_mapping = as.list(.from_SCE_guess_obspvarp(sce, colPairs)),
      varp_mapping = as.list(.from_SCE_guess_obspvarp(sce, rowPairs)),
      uns_mapping = as.list(.from_SCE_guess_all(sce, S4Vectors::metadata))
    )
  )
})

test_that("from_SCE works with unnamed mappings", {
  expect_no_error(
    as_AnnData(
      sce,
      layers_mapping = unname(.from_SCE_guess_layers(sce, NULL)),
      obs_mapping = unname(.from_SCE_guess_all(sce, colData)),
      var_mapping = unname(.from_SCE_guess_all(sce, rowData)),
      obsm_mapping = unname(.from_SCE_guess_obsm(sce)),
      varm_mapping = unname(.from_SCE_guess_varm(sce)),
      obsp_mapping = unname(.from_SCE_guess_obspvarp(sce, colPairs)),
      varp_mapping = unname(.from_SCE_guess_obspvarp(sce, rowPairs)),
      uns_mapping = unname(.from_SCE_guess_all(sce, S4Vectors::metadata))
    )
  )
})
