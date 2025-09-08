dummy <- generate_dataset(10L, 20L)

# GETTERS ----------------------------------------------------------------
test_that("Creating InMemoryAnnData works", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var
  )

  # trackstatus: class=InMemoryAnnData, feature=test_get_X, status=done
  # X should now have dimnames added on-the-fly
  expected_x <- dummy$X
  dimnames(expected_x) <- list(rownames(dummy$obs), rownames(dummy$var))
  expect_equal(ad$X, expected_x)
  # trackstatus: class=InMemoryAnnData, feature=test_get_obs, status=done
  expect_equal(ad$obs, dummy$obs)
  # trackstatus: class=InMemoryAnnData, feature=test_get_var, status=done
  expect_equal(ad$var, dummy$var)
  # trackstatus: class=InMemoryAnnData, feature=test_get_obs_names, status=done
  expect_equal(ad$obs_names, rownames(dummy$obs))
  # trackstatus: class=InMemoryAnnData, feature=test_get_var_names, status=done
  expect_equal(ad$var_names, rownames(dummy$var))
  expect_identical(ad$shape(), c(10L, 20L))
})

test_that("Creating InMemoryAnnData works with empty obs", {
  ad <- AnnData(
    obs = data.frame(),
    var = dummy$var
  )
  expect_identical(ad$shape(), c(0L, 20L))
})

test_that("Creating InMemoryAnnData works with empty var", {
  ad <- AnnData(
    obs = dummy$obs,
    var = data.frame()
  )
  expect_identical(ad$shape(), c(10L, 0L))
})

test_that("Creating AnnData works with only X, no obs or var", {
  X <- dummy$X
  dimnames(X) <- list(
    rownames(dummy$obs),
    rownames(dummy$var)
  )
  ad <- AnnData(X = X)
  expect_identical(ad$shape(), c(10L, 20L))
  expect_identical(ad$obs_names, rownames(dummy$obs))
  expect_identical(ad$var_names, rownames(dummy$var))
})

# trackstatus: class=InMemoryAnnData, feature=test_get_layers, status=done
test_that("'layers' works", {
  ## layers test helper function
  layers_test <- function(obs, var, layers) {
    expect_no_condition({
      ad <- AnnData(obs = obs, var = var, layers = layers)
    })
    # Expected layers should have dimnames added
    expected_layers <- layers
    if (!is.null(expected_layers)) {
      for (i in seq_along(expected_layers)) {
        if (!is.null(expected_layers[[i]])) {
          dimnames(expected_layers[[i]]) <- list(rownames(obs), rownames(var))
        }
      }
    }
    expect_identical(ad$layers, expected_layers)
  }

  obs <- var <- data.frame()
  layers_test(obs, var, NULL)
  layers_test(obs, var, setNames(list(), character()))
  layers_test(obs, var, list(A = matrix(0, 0, 0)))
  ## must be a named list
  expect_error(AnnData(obs, var, list()))

  obs <- data.frame(row.names = letters[1:3])
  var <- data.frame(row.names = LETTERS[1:5])
  layers_test(obs, var, NULL)
  layers_test(obs, var, list(A = matrix(0, 3, 5)))
  layers_test(obs, var, list(A = matrix(0, 3, 5), B = matrix(1, 3, 5)))

  ## must be a named list
  layers <- list(matrix(0, 3, 5))
  expect_error(AnnData(obs = obs, var = var, layers = layers))
  ## matching dimensions
  layers <- list(A = matrix(0, 0, 0))
  expect_error(AnnData(obs = obs, var = var, layers = layers))
  layers <- list(A = matrix(0, 3, 5), B = matrix(1, 5, 3))
  expect_error(AnnData(obs = obs, var = var, layers = layers))
})

test_that("*_keys() works", {
  obs <- var <- data.frame()
  ad <- AnnData(obs = obs, var = var)
  expect_identical(ad$layers_keys(), NULL)
  expect_identical(ad$obs_keys(), character())
  expect_identical(ad$var_keys(), character())

  layers <- setNames(list(), character())
  ad <- AnnData(obs = obs, var = var, layers = layers)
  expect_identical(ad$layers_keys(), character())

  layers <- list(A = matrix(0, 3, 5), B = matrix(1, 3, 5))
  obs <- data.frame(row.names = letters[1:3], x = 1:3)
  var <- data.frame(row.names = letters[1:5], y = 1:5)
  ad <- AnnData(obs = obs, var = var, layers = layers)
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
    var = dummy$var
  )

  X2 <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = .1)
  ad$X <- X2

  # X should have dimnames added on-the-fly
  expected_x2 <- X2
  dimnames(expected_x2) <- list(ad$obs_names, ad$var_names)
  expect_equal(ad$X, expected_x2)

  # change row in X
  ad$X[2, ] <- 10
  expected_row <- rep(10, 20L)
  names(expected_row) <- ad$var_names
  expect_equal(ad$X[2, ], expected_row)

  # change column in X
  ad$X[, 3] <- 5
  expected_col <- rep(5, 10L)
  names(expected_col) <- ad$obs_names
  expect_equal(ad$X[, 3], expected_col)
})

# trackstatus: class=InMemoryAnnData, feature=test_set_obs, status=done
test_that("write to obs", {
  ad <- AnnData(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var
  )

  obs2 <- data.frame(
    foo = sample(letters, 10, replace = TRUE),
    bar = sample.int(4, 10, replace = TRUE),
    zing = sample(c(TRUE, FALSE), 10, replace = TRUE)
  )
  ad$obs <- obs2

  expect_equal(ncol(ad$obs), 3)
  expect_equal(nrow(ad$obs), 10)
  # obs should now have rownames added on-the-fly
  expected_obs2 <- obs2
  rownames(expected_obs2) <- ad$obs_names
  expect_equal(ad$obs, expected_obs2)

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
    var = dummy$var
  )

  var2 <- data.frame(
    foo = sample(letters, 20L, replace = TRUE),
    bar = sample.int(4, 20L, replace = TRUE),
    zing = sample(c(TRUE, FALSE), 20L, replace = TRUE)
  )
  ad$var <- var2

  expect_equal(ncol(ad$var), 3)
  expect_equal(nrow(ad$var), 20)
  # var should now have rownames added on-the-fly
  expected_var2 <- var2
  rownames(expected_var2) <- ad$var_names
  expect_equal(ad$var, expected_var2)

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
    var = dummy$var
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
    var = dummy$var
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
    layers = dummy$layers
  )

  ## element assignment
  ad$layers[["X2"]] <- dummy$X + 2
  expected_x2 <- dummy$X + 2
  dimnames(expected_x2) <- list(ad$obs_names, ad$var_names)
  expect_equal(ad$layers[["X2"]], expected_x2)

  ad$layers[["X4"]] <- dummy$X + 4
  expected_x4 <- dummy$X + 4
  dimnames(expected_x4) <- list(ad$obs_names, ad$var_names)
  expect_equal(ad$layers[["X4"]], expected_x4)

  ## list assignment
  ad$layers <- rev(dummy$layers)
  expected_layers <- rev(dummy$layers)
  for (i in seq_along(expected_layers)) {
    if (!is.null(expected_layers[[i]])) {
      dimnames(expected_layers[[i]]) <- list(ad$obs_names, ad$var_names)
    }
  }
  expect_equal(ad$layers, expected_layers)

  ## remove
  ad$layers <- NULL
  expect_true(is.null(ad$layers))
  ## add to NULL
  ad$layers <- dummy$layers
  expected_layers <- dummy$layers
  for (i in seq_along(expected_layers)) {
    if (!is.null(expected_layers[[i]])) {
      dimnames(expected_layers[[i]]) <- list(ad$obs_names, ad$var_names)
    }
  }
  expect_identical(ad$layers, expected_layers)
})
