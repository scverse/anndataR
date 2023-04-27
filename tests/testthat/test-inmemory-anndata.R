library(Matrix)

test_that("create inmemory anndata", {
  # construct dummy X
  X <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = .1)

  # construct dummy obs
  obs <- data.frame(
    cell_type = sample(c("tcell", "bcell"), 10, replace = TRUE),
    cluster = sample.int(3, 10, replace = TRUE)
  )

  # construct dummy obs names
  obs_names <- paste0("cell_", seq_len(10))

  # construct dummy var
  var <- data.frame(
    geneinfo = sample(c("a", "b", "c"), 20, replace = TRUE)
  )

  # construct dummy var names
  var_names <- paste0("gene", seq_len(20))

  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  expect_equal(ad$X, X)
  expect_equal(ad$obs, obs)
  expect_equal(ad$var, var)
  expect_identical(ad$shape(), c(10L, 20L))
})

test_that("InMemoryAnnData$new() fails gracefully", {
  expect_error(InMemoryAnnData$new())
  expect_error(InMemoryAnnData$new(obs = data.frame(x = 1:3)))
  expect_error(InMemoryAnnData$new(var = data.frame(x = 1:3)))
})


test_that("InMemoryAnnData$new() validates dimensions", {
  X <- matrix(0L, 3L, 5L)
  expect_error(InMemoryAnnData$new(X = X, obs = data.frame(x = 1:5)))
  expect_error(InMemoryAnnData$new(X = X, var = data.frame(x = 1:3)))
})