dummy <- dummy_data(10L, 20L)

file <- system.file("extdata", "example.h5ad", package = "anndataR")
adata <- read_h5ad(file, to = "InMemoryAnnData")

# GETTERS ----------------------------------------------------------------
test_that("create inmemory anndata", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
  )

  # trackstatus: class=InMemoryAnnData, feature=test_get_X, status=done
  expect_equal(ad$X, dummy$X)
  # trackstatus: class=InMemoryAnnData, feature=test_get_obs, status=done
  expect_equal(ad$obs, dummy$obs)
  # trackstatus: class=InMemoryAnnData, feature=test_get_var, status=done
  expect_equal(ad$var, dummy$var)
  # trackstatus: class=InMemoryAnnData, feature=test_get_obs_names, status=done
  expect_equal(ad$obs_names, dummy$obs_names)
  # trackstatus: class=InMemoryAnnData, feature=test_get_var_names, status=done
  expect_equal(ad$var_names, dummy$var_names)
  expect_identical(ad$shape(), c(10L, 20L))
})

test_that("with empty obs", {
  ad <- AnnData(
    obs = data.frame(),
    var = dummy$var,
    obs_names = character(0),
    var_names = dummy$var_names
  )
  expect_identical(ad$shape(), c(0L, 20L))
})

test_that("with empty var", {
  ad <- AnnData(
    obs = dummy$obs,
    var = data.frame(),
    obs_names = dummy$obs_names,
    var_names = character(0)
  )
  expect_identical(ad$shape(), c(10L, 0L))
})

test_that("AnnData() fails gracefully", {
  expect_error(AnnData())
  expect_error(AnnData(obs = data.frame(x = 1:3)))
  expect_error(AnnData(var = data.frame(x = 1:3)))
})

test_that("InMemoryAnnData$new produces a warning if rownames are found", {
  # check X with rownames
  x_with_rownames <- dummy$X
  rownames(x_with_rownames) <- dummy$obs_names

  expect_warning({
    AnnData(
      X = x_with_rownames,
      obs = dummy$obs,
      var = dummy$var,
      obs_names = dummy$obs_names,
      var_names = dummy$var_names
    )
  })

  # check X with obsnames
  x_with_colnames <- dummy$X
  colnames(x_with_colnames) <- dummy$var_names

  expect_warning({
    AnnData(
      X = x_with_colnames,
      obs = dummy$obs,
      var = dummy$var,
      obs_names = dummy$obs_names,
      var_names = dummy$var_names
    )
  })

  # check obs with rownames
  obs_with_rownames <- dummy$obs
  rownames(obs_with_rownames) <- dummy$obs_names

  expect_warning({
    AnnData(
      X = dummy$X,
      obs = obs_with_rownames,
      var = dummy$var,
      obs_names = dummy$obs_names,
      var_names = dummy$var_names
    )
  })

  # check var with rownames
  var_with_rownames <- dummy$var
  rownames(var_with_rownames) <- dummy$var_names

  expect_warning({
    AnnData(
      X = dummy$X,
      obs = dummy$obs,
      var = var_with_rownames,
      obs_names = dummy$obs_names,
      var_names = dummy$var_names
    )
  })
})


# trackstatus: class=InMemoryAnnData, feature=test_get_layers, status=done
test_that("'layers' works", {
  ## layers test helper function
  layers_test <- function(obs_names, var_names, layers) {
    expect_no_condition({
      ad <- AnnData(obs_names = obs_names, var_names = var_names, layers = layers)
    })
    expect_identical(ad$layers, layers)
  }

  obs_names <- var_names <- character(0)
  layers_test(obs_names, var_names, NULL)
  layers_test(obs_names, var_names, setNames(list(), character()))
  layers_test(obs_names, var_names, list(A = matrix(0, 0, 0)))
  ## must be a named list
  expect_error(AnnData(obs_names = obs_names, var_names = var_names, list()))

  obs_names <- letters[1:3]
  var_names <- LETTERS[1:5]
  layers_test(obs_names, var_names, NULL)
  layers_test(obs_names, var_names, list(A = matrix(0, 3, 5)))
  layers_test(obs_names, var_names, list(A = matrix(0, 3, 5), B = matrix(1, 3, 5)))

  ## must be a named list
  layers <- list(matrix(0, 3, 5))
  expect_error(AnnData(obs_names = obs_names, var_names = var_names, layers = layers))
  ## non-trivial names
  layers <- list(A = matrix(0, 3, 5), matrix(1, 3, 5))
  expect_error(AnnData(obs_names = obs_names, var_names = var_names, layers = layers))
  ## matching dimensions
  layers <- list(A = matrix(0, 0, 0))
  expect_error(AnnData(obs_names = obs_names, var_names = var_names, layers = layers))
  layers <- list(A = matrix(0, 3, 5), B = matrix(1, 5, 3))
  expect_error(AnnData(obs_names = obs_names, var_names = var_names, layers = layers))
})

test_that("reading obsm works", {
  obsm <- adata$obsm
  expect_true(is.list(obsm), "list")
  expect_equal(
    names(obsm),
    c("X_pca", "X_umap")
  )
})

test_that("reading varm works", {
  varm <- adata$varm
  expect_true(is.list(varm), "list")
  expect_equal(
    names(varm),
    c("PCs")
  )
})

test_that("*_keys() works", {
  obs_names <- var_names <- character(0)
  ad <- AnnData(obs_names = obs_names, var_names = var_names)
  expect_identical(ad$layers_keys(), NULL)
  expect_identical(ad$obs_keys(), character())
  expect_identical(ad$var_keys(), character())

  layers <- setNames(list(), character())
  ad <- AnnData(obs_names = obs_names, var_names = var_names, layers = layers)
  expect_identical(ad$layers_keys(), character())

  obs_names <- letters[1:3]
  var_names <- letters[1:5]
  layers <- list(A = matrix(0, 3, 5), B = matrix(1, 3, 5))
  obs <- data.frame(x = 1:3)
  var <- data.frame(y = 1:5)
  ad <- AnnData(obs_names = obs_names, var_names = var_names, obs = obs, var = var, layers = layers)
  expect_identical(ad$layers_keys(), c("A", "B"))
  expect_identical(ad$obs_keys(), "x")
  expect_identical(ad$var_keys(), "y")
})

# SETTERS -----------------------------------------------------------------

# trackstatus: class=InMemoryAnnData, feature=test_set_X, status=done
test_that("write to X", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
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

# trackstatus: class=InMemoryAnnData, feature=test_set_obs, status=done
test_that("write to obs", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
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

# trackstatus: class=InMemoryAnnData, feature=test_set_var, status=done
test_that("write to var", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
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

# trackstatus: class=InMemoryAnnData, feature=test_set_obs_names, status=done
test_that("write to obs_names", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
  )

  obs_names2 <- letters[1:10]
  ad$obs_names <- obs_names2

  expect_equal(ad$obs_names, obs_names2)

  # change value in obs_names
  ad$obs_names[2:3] <- c("foo", "bar")
  expect_equal(ad$obs_names[1:4], c("a", "foo", "bar", "d"))
})

# trackstatus: class=InMemoryAnnData, feature=test_set_var_names, status=done
test_that("write to var_names", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names
  )

  var_names2 <- LETTERS[1:20]
  ad$var_names <- var_names2

  expect_equal(ad$var_names, var_names2)

  # change value in var_names
  ad$var_names[2:3] <- c("foo", "bar")
  expect_equal(ad$var_names[1:4], c("A", "foo", "bar", "D"))
})

# trackstatus: class=InMemoryAnnData, feature=test_set_layers, status=done
test_that("write to layers", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names,
    layers = dummy$layers
  )

  ## element assignment
  ad$layers[["X2"]] <- dummy$X + 2
  expect_equal(ad$layers[["X2"]], dummy$X + 2)

  ad$layers[["X4"]] <- dummy$X + 4
  expect_equal(ad$layers[["X4"]], dummy$X + 4)

  ## list assignment
  ad$layers <- rev(dummy$layers)
  expect_equal(ad$layers, rev(dummy$layers))

  ## remove
  ad$layers <- NULL
  expect_true(is.null(ad$layers))
  ## add to NULL
  ad$layers <- dummy$layers
  expect_identical(ad$layers, dummy$layers)
})
