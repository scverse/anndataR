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

layers = list(
  X2 = X * 2,
  X3 = X * 3
)

# GETTERS ----------------------------------------------------------------
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
  x_with_rownames <- X
  rownames(x_with_rownames) <- obs_names

  expect_warning({
    InMemoryAnnData$new(
      X = x_with_rownames,
      obs = obs,
      var = var,
      obs_names = obs_names,
      var_names = var_names
    )
  })

  # check X with obsnames
  x_with_colnames <- X
  colnames(x_with_colnames) <- var_names

  expect_warning({
    InMemoryAnnData$new(
      X = x_with_colnames,
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
  layers_test(obs, var, list(A = matrix(0, 3, 5), B = matrix(1, 3, 5)))

  ## must be a named list
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = list()))
  layers <- list(matrix(0, 3, 5))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
  ## non-trivial names
  layers <- list(A = matrix(0, 3, 5), matrix(1, 3, 5))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
  ## matching dimensions
  layers <- list(A = matrix(0, 0, 0))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
  layers <- list(A = matrix(0, 3, 5), B = matrix(1, 5, 3))
  expect_error(InMemoryAnnData$new(obs = obs, var = var, layers = layers))
})

test_that("*_keys() works", {
  obs <- var <- data.frame()
  ad <- InMemoryAnnData$new(obs = obs, var = var)
  expect_identical(ad$layers_keys(), NULL)
  expect_identical(ad$obs_keys(), character())
  expect_identical(ad$var_keys(), character())

  layers <- setNames(list(), character())
  ad <- InMemoryAnnData$new(obs = obs, var = var, layers = layers)
  expect_identical(ad$layers_keys(), character())

  layers <- list(A = matrix(0, 3, 5), B = matrix(1, 3, 5))
  obs <- data.frame(x = 1:3)
  var <- data.frame(y = 1:5)
  ad <- InMemoryAnnData$new(obs = obs, var = var, layers = layers)
  expect_identical(ad$layers_keys(), c("A", "B"))
  expect_identical(ad$obs_keys(), "x")
  expect_identical(ad$var_keys(), "y")
})

# SETTERS -----------------------------------------------------------------


test_that("write to X", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  X2 <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = .1)
  ad$X <- X2

  expect_equal(ad$X, X2)

  # change row in X
  ad$X[2, ] <- 10
  expect_equal(ad$X[2, ], rep(10, 20L))

  # change column in X
  ad$X[, 3] <- 5
  expect_equal(ad$X[, 3], rep(5, 10L))
})


test_that("write to obs", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  obs2 <- data.frame(
    foo = sample(letters, 10, replace = TRUE),
    bar = sample.int(4, 10, replace = TRUE),
    zing = sample(c(TRUE, FALSE), 10, replace = TRUE)
  )
  ad$obs <- obs2

  expect_equal(ncol(ad$obs), 3)
  expect_equal(nrow(ad$obs), 10)
  expect_equal(ad$obs, obs2)

  # change row in obs
  obs2row <- data.frame(foo = "a", bar = 3, zing = FALSE)
  ad$obs[2, ] <- obs2row
  expect_equal(ad$obs[2, ], obs2row, ignore_attr = TRUE)

  # change column in obs
  ad$obs[, 3] <- FALSE
  expect_equal(ad$obs[, 3], rep(FALSE, 10L))
})

test_that("write to var", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  var2 <- data.frame(
    foo = sample(letters, 20L, replace = TRUE),
    bar = sample.int(4, 20L, replace = TRUE),
    zing = sample(c(TRUE, FALSE), 20L, replace = TRUE)
  )
  ad$var <- var2

  expect_equal(ncol(ad$var), 3)
  expect_equal(nrow(ad$var), 20)
  expect_equal(ad$var, var2)

  # change row in var
  var2row <- data.frame(foo = "a", bar = 3, zing = FALSE)
  ad$var[2, ] <- var2row
  expect_equal(ad$var[2, ], var2row, ignore_attr = TRUE)

  # change column in var
  ad$var[, 3] <- FALSE
  expect_equal(ad$var[, 3], rep(FALSE, 20L))
})


test_that("write to obs_names", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  obs_names2 <- letters[1:10]
  ad$obs_names <- obs_names2

  expect_equal(ad$obs_names, obs_names2)

  # change value in obs_names
  ad$obs_names[2:3] <- c("foo", "bar")
  expect_equal(ad$obs_names[1:4], c("a", "foo", "bar", "d"))
})

test_that("write to var_names", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  var_names2 <- LETTERS[1:20]
  ad$var_names <- var_names2

  expect_equal(ad$var_names, var_names2)

  # change value in var_names
  ad$var_names[2:3] <- c("foo", "bar")
  expect_equal(ad$var_names[1:4], c("A", "foo", "bar", "D"))
})

test_that("write to layers", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names,
    layers = layers
  )

  ## element assignment
  ad$layers[["X2"]] <- X + 2
  expect_equal(ad$layers[["X2"]], X + 2)

  ad$layers[["X4"]] <- X + 4
  expect_equal(ad$layers[["X4"]], X + 4)

  ## list assignment
  ad$layers <- rev(layers)
  expect_equal(ad$layers, rev(layers))

  ## remove
  ad$layers <- NULL
  expect_true(is.null(ad$layers))
  ## add to NULL
  ad$layers <- layers
  expect_identical(ad$layers, layers)
})
