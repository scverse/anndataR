skip_if_not_installed("SingleCellExperiment")
library(SingleCellExperiment)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)

# TODO: Build an SCE rather than converting
sce <- ad$as_SingleCellExperiment()
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
for (layer_key in assayNames(sce)) {
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

test_that("from_SCE works with empty mappings", {
  expect_warning(as_AnnData(sce, layers_mapping = NULL), "argument is empty")
  expect_warning(as_AnnData(sce, layers_mapping = c()), "argument is empty")
})

obs_mapping <- c(obs1 = "character", obs2 = "numeric")
var_mapping <- c(var1 = "character", var2 = "numeric")
layers_mapping <- c(counts = "integer_matrix", data = "numeric_matrix")

ad_partial <- as_AnnData(
  sce,
  layers_mapping = layers_mapping,
  obs_mapping = obs_mapping,
  var_mapping = var_mapping,
  obsm_mapping = list(),
  varm_mapping = list(),
  obsp_mapping = list(),
  varp_mapping = list(),
  uns_mapping = list()
)

for (obs_key in names(obs_mapping)) {
  obs_from <- obs_mapping[[obs_key]]
  test_that(paste0("from_SCE retains obs key: ", obs_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("obs"),
      dtype = obs_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obs_key %in% colnames(ad_partial$obs))
    expect_equal(
      ad_partial$obs[[obs_key]],
      colData(sce)[[obs_from]],
      info = paste0("obs_key: ", obs_key)
    )
  })
}

for (var_key in names(var_mapping)) {
  var_from <- var_mapping[[var_key]]
  test_that(paste0("from_SCE retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "from_SCE",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(var_key %in% colnames(ad_partial$var))
    expect_equal(
      ad_partial$var[[var_key]],
      rowData(sce)[[var_from]],
      info = paste0("var_key: ", var_key)
    )
  })
}

test_that(paste0("from_SCE does not copy unmapped structs"), {
  msg <- message_if_known(
    backend = "from_SCE",
    slot = c("obsm", "varm", "obsp", "varp", "uns"),
    dtype = "unmapped",
    process = "convert",
    known_issues = known_issues
  )
  skip_if(!is.null(msg), message = msg)

  expect_true(is.null(ad_partial$obsm))
  expect_true(is.null(ad_partial$varm))
  expect_true(is.null(ad_partial$obsp))
  expect_true(is.null(ad_partial$varp))
  expect_true(is.null(ad_partial$uns))
})
