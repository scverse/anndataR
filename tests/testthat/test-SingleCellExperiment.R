test_that("to_SingleCellExperiment() works", {
  ad <- InMemoryAnnData$new(
    X = matrix(1:5, 3L, 5L),
    obs = data.frame(cell = 1:3),
    var = data.frame(gene = 1:5),
    obs_names = LETTERS[1:3],
    var_names = letters[1:5]
  )
  ad0 <- InMemoryAnnData$new(
    obs_names = LETTERS[1:6],
    var_names = letters[1:8]
  )

  # conversion works
  expect_no_error(sce <- to_SingleCellExperiment(ad))
  expect_no_error(sce0 <- to_SingleCellExperiment(ad0))
  expect_true(validObject(sce))
  expect_true(validObject(sce0))

  expect_identical(dim(sce), rev(ad$shape()))
  expect_identical(dim(sce0), rev(ad0$shape()))

  # trackstatus: class=SingleCellExperiment, feature=test_get_obs_names, status=done
  # trackstatus: class=SingleCellExperiment, feature=test_get_var_names, status=done
  expect_identical(dimnames(sce), list(letters[1:5], LETTERS[1:3]))
  expect_identical(dimnames(sce0), list(letters[1:8], LETTERS[1:6]))

  # trackstatus: class=SingleCellExperiment, feature=test_get_var, status=done
  var_ <- as.data.frame(SummarizedExperiment::rowData(sce))
  rownames(var_) <- NULL
  expect_identical(var_, ad$var)
  var0_ <- as.data.frame(SummarizedExperiment::rowData(sce0))
  rownames(var0_) <- NULL
  expect_identical(var0_, ad0$var)

  # trackstatus: class=SingleCellExperiment, feature=test_get_obs, status=done
  obs_ <- as.data.frame(SummarizedExperiment::colData(sce))
  rownames(obs_) <- NULL
  expect_identical(obs_, ad$obs)
  obs0_ <- as.data.frame(SummarizedExperiment::colData(sce0))
  rownames(obs0_) <- NULL
  expect_identical(obs0_, ad0$obs)

  # trackstatus: class=SingleCellExperiment, feature=test_get_X, status=done
  expect_identical(
    SummarizedExperiment::assay(sce, withDimnames = FALSE),
    t(ad$X)
  )
  expect_error(
    SummarizedExperiment::assay(sce0, withDimnames = FALSE)
  )
})

test_that("from_SingleCellExperiment() works", {
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

  rownames(obs) <- NULL
  rownames(var) <- NULL

  ad0 <- from_SingleCellExperiment(sce0, "InMemory")
  ad <- from_SingleCellExperiment(sce, "InMemory")

  # trackstatus: class=SingleCellExperiment, feature=test_set_X, status=done
  expect_identical(ad0$X, NULL)
  expect_identical(ad$X, x)
  # trackstatus: class=SingleCellExperiment, feature=test_set_obs, status=done
  expect_identical(ad0$obs, data.frame())
  expect_identical(ad$obs, obs)
  # trackstatus: class=SingleCellExperiment, feature=test_set_var, status=done
  expect_identical(ad0$var, data.frame())
  expect_identical(ad$var, var)
  # trackstatus: class=SingleCellExperiment, feature=test_set_obs_names, status=done
  expect_identical(ad0$obs_names, character(0))
  expect_identical(ad$obs_names, dimnames[[2]])
  # trackstatus: class=SingleCellExperiment, feature=test_set_var_names, status=done
  expect_identical(ad0$var_names, character(0))
  expect_identical(ad$var_names, dimnames[[1]])
  # trackstatus: class=SingleCellExperiment, feature=test_set_layers, status=done
  layers0 <- list()
  names(layers0) <- character()
  expect_identical(ad0$layers, layers0)
  expect_identical(ad$layers, layers)
})
