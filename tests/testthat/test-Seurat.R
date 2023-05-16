library(Matrix)

dummy <- dummy_data(10L, 20L)

test_that("to_Seurat with inmemoryanndata", {
  ad <- InMemoryAnnData$new(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
  )

  seu <- suppressWarnings(to_Seurat(ad))

  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)
  expect_equal(rownames(seu), dummy$var_names)
  expect_equal(colnames(seu), dummy$obs_names)

  # todo: check whether X can be retrieved
})

# todo: check substitution from _ to -

test_that("to_Seurat() fails gracefully", {
  expect_error(to_Seurat(), regexp = "obj.*is missing")
  expect_error(to_Seurat("foo"), regexp = "AbstractAnnData.*not TRUE")
})
