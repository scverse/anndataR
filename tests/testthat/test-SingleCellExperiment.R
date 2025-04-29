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

# TODO gracefully failing

skip_if_not_installed("SingleCellExperiment")
library(SingleCellExperiment)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

# TODO: Build an SCE rather than converting
sce <- ad$to_SingleCellExperiment()
ad <- from_SingleCellExperiment(sce)

test_that("from_SCE retains observatoins and features", {
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

  print(ad)
  print(sce)

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

test_that("from_SingleCellExperiment() works with Zarr", {
  ## 0-dimensioned
  sce0 <- SingleCellExperiment::SingleCellExperiment()
  dimnames(sce0) <- list(character(0), character(0))

  ## complete
  x <- matrix(1:15, 3, 5)
  layers <- list(A = matrix(15:1, 3, 5), B = matrix(LETTERS[1:15], 3, 5))
  obs <- data.frame(cell = 1:3, row.names = LETTERS[1:3])
  var <- data.frame(gene = 1:5, row.names = letters[1:5])
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = lapply(c(list(x), layers), t),
    colData = obs,
    rowData = var
  )
  dimnames <- dimnames(sce)

  store0 <- pizzarr::MemoryStore$new()
  store <- pizzarr::MemoryStore$new()

  ad0 <- from_SingleCellExperiment(sce0, "ZarrAnnData", store = store0)
  ad <- from_SingleCellExperiment(sce, "ZarrAnnData", store = store)

  # trackstatus: class=SingleCellExperiment, feature=test_set_X, status=done
  expect_identical(ad0$X, NULL)
  expect_identical(ad$X, x)
  # trackstatus: class=SingleCellExperiment, feature=test_set_obs, status=done
  expect_identical(ad0$obs, data.frame(row.names = character(0)))
  expect_identical(ad$obs, obs)
  # trackstatus: class=SingleCellExperiment, feature=test_set_var, status=done
  expect_identical(ad0$var, data.frame(row.names = character(0)))
  expect_identical(ad$var, var)
  # trackstatus: class=SingleCellExperiment, feature=test_set_obs_names, status=done
  expect_identical(rownames(ad0$obs), character(0))
  expect_identical(rownames(ad$obs), dimnames[[2]])
  # trackstatus: class=SingleCellExperiment, feature=test_set_var_names, status=done
  expect_identical(rownames(ad0$var), character(0))
  expect_identical(rownames(ad$var), dimnames[[1]])
  # trackstatus: class=SingleCellExperiment, feature=test_set_layers, status=done
  layers0 <- list()
  expect_identical(ad0$layers, layers0)
  expect_identical(ad$layers, layers)
})
