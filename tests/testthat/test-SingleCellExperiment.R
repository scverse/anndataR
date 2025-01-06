dummy <- generate_dataset(10L, 20L)

test_that("to_SingleCellExperiment with inmemoryanndata", {
  library(SingleCellExperiment)

  ad <- generate_dataset(n_obs = 10L, n_var = 20L, format = "AnnData")

  sce <- ad$to_SingleCellExperiment()

  expect_equal(nrow(sce), 20)
  expect_equal(ncol(sce), 10)

  # trackstatus: class=SingleCellExperiment, feature=test_get_var_names, status=done
  expect_equal(rownames(sce), rownames(dummy$var))

  # trackstatus: class=SingleCellExperiment, feature=test_get_obs_names, status=done
  expect_equal(colnames(sce), rownames(dummy$obs))

  # check whether all obs keys are found in the sce metadata
  # trackstatus: class=SingleCellExperiment, feature=test_get_obs, status=done
  for (obs_key in colnames(dummy$obs)) {
    expect_true(obs_key %in% colnames(colData(sce)))
    expect_equal(colData(sce)[[obs_key]], dummy$obs[[obs_key]], info = paste0("obs_key: ", obs_key))
  }

  # check whether all var keys are found in the sce assay metadata
  # trackstatus: class=SingleCellExperiment, feature=test_get_var, status=done
  for (var_key in colnames(dummy$var)) {
    expect_true(var_key %in% colnames(rowData(sce)))
    expect_equal(rowData(sce)[[var_key]], dummy$var[[var_key]], info = paste0("var_key: ", var_key))
  }

  # check whether layers are found in the sce assays
  for (layer_key in names(dummy$layers)) {
    expect_true(layer_key %in% names(assays(sce)))
    expect_true(all.equal(assay(sce, layer_key), t(dummy$layers[[layer_key]]), check.attributes = FALSE), info = paste0("layer_key: ", layer_key)) ## matrix dimensions and dimnames are different
  }

  # check whether all obsp keys are found in the colPairs
  for (obsp_key in names(dummy$obsp)) {
    expect_true(obsp_key %in% names(colPairs(sce)))
    expect_equal(colData(sce)[[obsp_key]], dummy$obsp[[obsp_key]], info = paste0("obsm_key: ", obsm_key))
  }

  # check whether all varp keys are found in the rowPairs
  for (varp_key in names(dummy$varp)) {
    expect_true(varp_key %in% names(rowPairs(sce)))
    expect_equal(rowData(sce)[[varp_key]], dummy$varp[[varp_key]], info = paste0("varm_key: ", varm_key))
  }

  # TODO: obsm keys? varm keys? --> test a reduction

})

# test_that("to_SingleCellExperiment() works", {
#   ad <- AnnData(
#     X = matrix(1:5, 3L, 5L),
#     obs = data.frame(row.names = letters[1:3], cell = 1:3),
#     var = data.frame(row.names = LETTERS[1:5], gene = 1:5)
#   )
#   ad0 <- AnnData(
#     obs = data.frame(row.names = letters[1:5]),
#     var = data.frame(row.names = LETTERS[1:10])
#   )

#   # conversion works
#   expect_no_error(sce <- to_SingleCellExperiment(ad))
#   expect_no_error(sce0 <- to_SingleCellExperiment(ad0))
#   expect_true(validObject(sce))
#   expect_true(validObject(sce0))

#   expect_identical(dim(sce), rev(ad$shape()))
#   expect_identical(dim(sce0), rev(ad0$shape()))

#   # trackstatus: class=SingleCellExperiment, feature=test_get_obs_names, status=done
#   # trackstatus: class=SingleCellExperiment, feature=test_get_var_names, status=done
#   expect_identical(dimnames(sce), list(LETTERS[1:5], letters[1:3]))
#   expect_identical(dimnames(sce0), list(LETTERS[1:10], letters[1:5]))

#   # trackstatus: class=SingleCellExperiment, feature=test_get_var, status=done
#   var_ <- as.data.frame(SummarizedExperiment::rowData(sce))
#   expect_identical(var_, ad$var)
#   var0_ <- as.data.frame(SummarizedExperiment::rowData(sce0))
#   expect_identical(var0_, ad0$var)

#   # trackstatus: class=SingleCellExperiment, feature=test_get_obs, status=done
#   obs_ <- as.data.frame(SummarizedExperiment::colData(sce))
#   expect_identical(obs_, ad$obs)
#   obs0_ <- as.data.frame(SummarizedExperiment::colData(sce0))
#   expect_identical(obs0_, ad0$obs)

#   # trackstatus: class=SingleCellExperiment, feature=test_get_X, status=done
#   expect_identical(
#     SummarizedExperiment::assay(sce, withDimnames = FALSE),
#     t(ad$X)
#   )
#   expect_error(
#     SummarizedExperiment::assay(sce0, withDimnames = FALSE)
#   )
# })

# test_that("from_SingleCellExperiment() works", {
#   ## 0-dimensioned
#   sce0 <- SingleCellExperiment::SingleCellExperiment()
#   dimnames(sce0) <- list(character(0), character(0))

#   ## complete
#   x <- matrix(1:15, 3, 5)
#   layers <- list(A = matrix(15:1, 3, 5), B = matrix(LETTERS[1:15], 3, 5))
#   obs <- data.frame(cell = 1:3, row.names = LETTERS[1:3])
#   var <- data.frame(gene = 1:5, row.names = letters[1:5])
#   sce <- SingleCellExperiment::SingleCellExperiment(
#     assays = lapply(c(list(x), layers), t),
#     colData = obs,
#     rowData = var
#   )
#   dimnames <- dimnames(sce)

#   ad0 <- from_SingleCellExperiment(sce0, "InMemory")
#   ad <- from_SingleCellExperiment(sce, "InMemory")

#   # trackstatus: class=SingleCellExperiment, feature=test_set_X, status=done
#   expect_identical(ad0$X, NULL)
#   expect_identical(ad$X, x)
#   # trackstatus: class=SingleCellExperiment, feature=test_set_obs, status=done
#   expect_identical(ad0$obs, data.frame(row.names = character(0)))
#   expect_identical(ad$obs, obs)
#   # trackstatus: class=SingleCellExperiment, feature=test_set_var, status=done
#   expect_identical(ad0$var, data.frame(row.names = character(0)))
#   expect_identical(ad$var, var)
#   # trackstatus: class=SingleCellExperiment, feature=test_set_obs_names, status=done
#   expect_identical(ad0$obs_names, character(0))
#   expect_identical(ad$obs_names, dimnames[[2]])
#   # trackstatus: class=SingleCellExperiment, feature=test_set_var_names, status=done
#   expect_identical(ad0$var_names, character(0))
#   expect_identical(ad$var_names, dimnames[[1]])
#   # trackstatus: class=SingleCellExperiment, feature=test_set_layers, status=done
#   layers0 <- list()
#   names(layers0) <- character()
#   expect_identical(ad0$layers, layers0)
#   expect_identical(ad$layers, layers)
# })
