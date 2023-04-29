test_that("to_SingleCellExperiment() works", {
  ad <- InMemoryAnnData$new(
    X = matrix(1:5, 3L, 5L),
    obs = data.frame(cell = 1:3),
    var = data.frame(gene = 1:5)
  )

  expect_no_error(sce <- to_SingleCellExperiment(ad))
  expect_true(validObject(sce))
  expect_identical(dim(sce), rev(ad$shape()))

  expect_identical(dimnames(sce), NULL)
  df <- as.data.frame(SummarizedExperiment::rowData(sce))
  rownames(df) <- NULL
  expect_identical(df, ad$var)
  df <- as.data.frame(SummarizedExperiment::colData(sce))
  rownames(df) <- NULL
  expect_identical(df, ad$obs)

  expect_identical(
    SummarizedExperiment::assay(sce, withDimnames = FALSE),
    t(ad$X)
  )

  ## dimnames() from *_names
  ad <- InMemoryAnnData$new(
    X = matrix(1:5, 3L, 5L),
    obs = data.frame(cell = 1:3),
    var = data.frame(gene = 1:5),
    obs_names = LETTERS[1:3],
    var_names = letters[1:5]
  )
  expect_no_condition(sce <- to_SingleCellExperiment(ad))
  expect_identical(dimnames(sce), list(letters[1:5], LETTERS[1:3]))
})
