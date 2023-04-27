library(Matrix)

test_that("create inmemory anndata", {
  # construct dummy X
  X <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = .1)

  # construct dummy obs
  obs <- data.frame(
    row.names = paste0("cell_", seq_len(10)),
    cell_type = sample(c("tcell", "bcell"), 10, replace = TRUE),
    cluster = sample.int(3, 10, replace = TRUE)
  )

  # construct dummy obs
  var <- data.frame(
    row.names = paste0("gene", seq_len(20)),
    geneinfo = sample(c("a", "b", "c"), 20, replace = TRUE)
  )

  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var
  )

  expect_equal(ad$X, X)
  expect_equal(ad$obs, obs)
  expect_equal(ad$var, var)
})
