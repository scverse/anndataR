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
  expect_identical(ad$dim(), dim(X))
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

test_that("dim() works", {
  ad <- InMemoryAnnData$new(shape = c(1L, 1L))
  expect_identical(ad$dim(), c(1L, 1L))

  ## from obs and var, or if that fails X
  ad <- InMemoryAnnData$new(obs = data.frame(x = 1:3), shape = c(3L, 0L))
  expect_identical(ad$dim(), c(3L, 0L))
  ad <- InMemoryAnnData$new(var = data.frame(x = 1:5), shape = c(0L, 5L))
  expect_identical(ad$dim(), c(0L, 5L))
  ad <- InMemoryAnnData$new(X = matrix(0L, 3L, 5L))
  expect_identical(ad$dim(), c(3L, 5L))
})
