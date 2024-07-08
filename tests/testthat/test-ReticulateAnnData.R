library(reticulate)

dummy <- generate_dataset(10L, 20L)

# GETTERS ----------------------------------------------------------------
test_that("test getters", {
  adata <- ReticulateAnnData$new(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    layers = dummy$layers
  )

  # trackstatus: class=ReticulateAnnData, feature=test_get_X, status=done
  expect_equal(adata$X, dummy$X)
  # trackstatus: class=ReticulateAnnData, feature=test_get_obs, status=done
  expect_equal(adata$obs, dummy$obs)
  # trackstatus: class=ReticulateAnnData, feature=test_get_var, status=done
  expect_equal(adata$var, dummy$var)
  expect_identical(adata$shape(), c(10L, 20L))
  # trackstatus: class=ReticulateAnnData, feature=test_get_obs_names, status=done
  expect_equal(adata$obs_names, dummy$obs_names)
  # trackstatus: class=ReticulateAnnData, feature=test_get_var_names, status=done
  expect_equal(adata$var_names, dummy$var_names)
  # trackstatus: class=ReticulateAnnData, feature=test_get_layers, status=done
  expect_equal(adata$layers, dummy$layers)
})


# SETTERS ----------------------------------------------------------------

# trackstatus: class=ReticulateAnnData, feature=test_set_X, status=done
test_that("write to X", {
  ad <- ReticulateAnnData$new(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var
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

# trackstatus: class=ReticulateAnnData, feature=test_set_obs, status=done
test_that("write to obs", {
  ad <- ReticulateAnnData$new(
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
  expect_equal(ad$obs, obs2)

  # change row in obs
  obs2row <- data.frame(foo = "a", bar = 3, zing = FALSE)
  ad$obs[2, ] <- obs2row
  expect_equal(ad$obs[2, ], obs2row, ignore_attr = TRUE)

  # change column in obs
  ad$obs[, 3] <- FALSE
  expect_equal(ad$obs[, 3], rep(FALSE, 10L))
})

# trackstatus: class=ReticulateAnnData, feature=test_set_var, status=done
test_that("write to var", {
  ad <- ReticulateAnnData$new(
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
  expect_equal(ad$var, var2)

  # change row in var
  var2row <- data.frame(foo = "a", bar = 3, zing = FALSE)
  ad$var[2, ] <- var2row
  expect_equal(ad$var[2, ], var2row, ignore_attr = TRUE)

  # change column in var
  ad$var[, 3] <- FALSE
  expect_equal(ad$var[, 3], rep(FALSE, 20L))
})

# trackstatus: class=ReticulateAnnData, feature=test_set_obs_names, status=done
test_that("write to obs_names", {
  ad <- ReticulateAnnData$new(
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

# trackstatus: class=ReticulateAnnData, feature=test_set_var_names, status=done
test_that("write to var_names", {
  ad <- ReticulateAnnData$new(
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

# trackstatus: class=ReticulateAnnData, feature=test_set_layers, status=done
test_that("write to layers", {
  ad <- ReticulateAnnData$new(
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    layers = dummy$layers
  )

  ## element assignment
  ad$layers[["X2"]] <- dummy$X + 2
  expect_equal(ad$layers[["X2"]], as.matrix(dummy$X + 2))

  ad$layers[["X4"]] <- dummy$X + 4
  expect_equal(ad$layers[["X4"]], as.matrix(dummy$X + 4))

  ## list assignment
  ad$layers <- rev(dummy$layers)
  expect_equal(ad$layers, rev(dummy$layers))

  ## remove
  ad$layers <- NULL
  expect_true(length(ad$layers) == 0)

  ## add to NULL
  ad$layers <- dummy$layers
  expect_identical(ad$layers, dummy$layers)
})
