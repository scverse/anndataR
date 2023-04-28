library(Matrix)

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

test_that("create inmemory anndata", {
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

test_that("with empty obs", {
  ad <- InMemoryAnnData$new(
    obs = data.frame(),
    var = var,
    obs_names = c(),
    var_names = var_names
  )
  expect_identical(ad$shape(), c(0L, 20L))
})

test_that("with empty var", {
  ad <- InMemoryAnnData$new(
    obs = obs,
    var = data.frame(),
    obs_names = obs_names,
    var_names = c()
  )
  expect_identical(ad$shape(), c(10L, 0L))
})

test_that("InMemoryAnnData$new() fails gracefully", {
  expect_error(InMemoryAnnData$new())
  expect_error(InMemoryAnnData$new(obs = data.frame(x = 1:3)))
  expect_error(InMemoryAnnData$new(var = data.frame(x = 1:3)))
})

test_that("InMemoryAnnData$new produces a warning if rownames are found", {
  # check X with rownames
  X_with_rownames <- X
  rownames(X_with_rownames) <- obs_names

  expect_warning({
    InMemoryAnnData$new(
      X = X_with_rownames,
      obs = obs,
      var = var,
      obs_names = obs_names,
      var_names = var_names
    )
  })

  # check X with obsnames
  X_with_colnames <- X
  colnames(X_with_colnames) <- var_names

  expect_warning({
    InMemoryAnnData$new(
      X = X_with_colnames,
      obs = obs,
      var = var,
      obs_names = obs_names,
      var_names = var_names
    )
  })

  # check obs with rownames
  obs_with_rownames <- obs
  rownames(obs_with_rownames) <- obs_names

  expect_warning({
    InMemoryAnnData$new(
      X = X,
      obs = obs_with_rownames,
      var = var,
      obs_names = obs_names,
      var_names = var_names
    )
  })

  # check var with rownames
  var_with_rownames <- var
  rownames(var_with_rownames) <- var_names

  expect_warning({
    InMemoryAnnData$new(
      X = X,
      obs = obs,
      var = var_with_rownames,
      obs_names = obs_names,
      var_names = var_names
    )
  })

})

test_that("'layers' works", {
  ## layers test helper function
  layers_test <- function(obs, var, layers) {
    expect_no_condition({
      ad <- InMemoryAnnData$new(obs = obs, var = var, layers = layers)
    })
    expect_identical(ad$layers, layers)
  }

  obs <- var <- data.frame()
  layers_test(obs, var, NULL)
  layers_test(obs, var, setNames(list(), character()))
  layers_test(obs, var, list(A = matrix(0, 0, 0)))
  ## must be a named list
  expect_error(InMemoryAnnData$new(obs = obs, var = var, list()))

  obs <- data.frame(x = 1:3)
  var <- data.frame(y = 1:5)
  layers_test(obs, var, NULL)
  layers_test(obs, var, list(A = matrix(0, 3, 5)))
  layers_test(obs, var, list(A = matrix(0, 3, 5), B = matrix(0, 3, 5)))

  ## must be a named list
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = list()))
  layers <- list(matrix(0, 3, 5))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
  ## non-trivial names
  layers <- list(A = matrix(0, 3, 5), matrix(0, 3, 5))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
  ## matching dimensions
  layers <- list(A = matrix(0, 0, 0))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
  layers <- list(A = matrix(0, 3, 5), B = matrix(0, 5, 3))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
})

test_that("*_keys() works", {
  obs <- var <- data.frame()
  ad <- InMemoryAnnData$new(obs = obs, var = var)
  expect_identical(ad$obs_keys(), character())
  expect_identical(ad$var_keys(), character())

  obs <- data.frame(x = 1:3)
  var <- data.frame(y = 1:5)
  ad <- InMemoryAnnData$new(obs = obs, var = var)
  expect_identical(ad$obs_keys(), "x")
  expect_identical(ad$var_keys(), "y")
})
