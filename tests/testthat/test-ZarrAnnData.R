skip_if_not_installed("Rarr")

file <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(file, exdir = td)
store <- file.path(td, "example_v2.zarr")

test_that("opening Zarr works", {
  adata <- ZarrAnnData$new(store, mode = "r")
  expect_true(inherits(adata, "ZarrAnnData"))
})

adata <- ZarrAnnData$new(store, mode = "r")

# GETTERS ----------------------------------------------------------------
# trackstatus: class=ZarrAnnData, feature=test_get_X, status=done
test_that("reading X works", {
  X <- adata$X
  expect_s4_class(X, "dgRMatrix")
  expect_equal(dim(X), c(50, 100))
})

# trackstatus: class=ZarrAnnData, feature=test_get_layers, status=done
test_that("reading layers works", {
  layers <- adata$layers
  expect_true(is.list(layers), "list")
  expect_equal(
    names(layers),
    c("counts", "csc_counts", "dense_X", "dense_counts")
  )
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

# trackstatus: class=ZarrAnnData, feature=test_get_obsp, status=done
test_that("reading obsp works", {
  obsp <- adata$obsp
  expect_true(is.list(obsp), "list")
  expect_equal(
    names(obsp),
    c("connectivities", "distances")
  )
})

# trackstatus: class=ZarrAnnData, feature=test_get_varp, status=done
test_that("reading varp works", {
  varp <- adata$varp
  expect_true(is.list(varp), "list")
  expect_equal(
    names(varp),
    c("test_varp")
  )
})

# trackstatus: class=ZarrAnnData, feature=test_get_obs, status=done
test_that("reading obs works", {
  obs <- adata$obs
  expect_s3_class(obs, "data.frame")
  expect_equal(
    colnames(obs),
    c(
      "Float",
      "FloatNA",
      "Int",
      "IntNA",
      "Bool",
      "BoolNA",
      "n_genes_by_counts",
      "log1p_n_genes_by_counts",
      "total_counts",
      "log1p_total_counts",
      "leiden"
    )
  )
})

# trackstatus: class=ZarrAnnData, feature=test_get_var, status=done
test_that("reading var works", {
  var <- adata$var
  expect_s3_class(var, "data.frame")
  expect_equal(
    colnames(var),
    c(
      "String",
      "n_cells_by_counts",
      "mean_counts",
      "log1p_mean_counts",
      "pct_dropout_by_counts",
      "total_counts",
      "log1p_total_counts",
      "highly_variable",
      "means",
      "dispersions",
      "dispersions_norm"
    )
  )
})

# trackstatus: class=ZarrAnnData, feature=test_get_obs_names, status=done
test_that("reading obs names works", {
  obs_names <- adata$obs_names
  expect_vector(obs_names, ptype = character(), size = 50)
})

# trackstatus: class=ZarrAnnData, feature=test_get_var_names, status=done
test_that("reading var names works", {
  var_names <- adata$var_names
  expect_vector(var_names, ptype = character(), size = 100)
})

# SETTERS ----------------------------------------------------------------
test_that("creating empty Zarr works", {
  empty_store <- tempfile(fileext = ".zarr")
  expect_silent(ZarrAnnData$new(empty_store))
  unlink(empty_store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_X, status=done
test_that("writing X works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(zarr$X <- X)
  unlink(store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_layers, status=done
test_that("writing layers works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(zarr$layers <- list(layer1 = X, layer2 = X))
  unlink(store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_obs, status=done
test_that("writing obs works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  obs <- data.frame(
    Letters = LETTERS[1:10],
    Numbers = 1:10,
    row.names = paste0("Row", 1:10)
  )
  zarr$obs <- obs
  expect_identical(zarr$obs_names, paste0("Row", 1:10))
  unlink(store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_var, status=done
test_that("writing var works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  var <- data.frame(
    Letters = LETTERS[1:20],
    Numbers = 1:20,
    row.names = paste0("Row", 1:20)
  )
  zarr$var <- var
  expect_identical(zarr$var_names, paste0("Row", 1:20))
  unlink(store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_obs_names, status=done
test_that("writing obs names works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  zarr$obs_names <- LETTERS[1:10]
  expect_identical(zarr$obs_names, LETTERS[1:10])
  unlink(store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_var_names, status=done
test_that("writing var names works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  zarr$var_names <- LETTERS[1:20]
  expect_identical(zarr$var_names, LETTERS[1:20])
  unlink(store, recursive = TRUE)
})

# trackstatus: class=ZarrAnnData, feature=test_set_obsm, status=done
test_that("writing obsm works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)
  obsm_x <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  zarr$obsm <- list(X = obsm_x)
  # obsm should now have rownames added on-the-fly
  expected_obsm_x <- obsm_x
  rownames(expected_obsm_x) <- zarr$obs_names
  expect_identical(zarr$obsm$X, expected_obsm_x)
})

# trackstatus: class=ZarrAnnData, feature=test_set_varm, status=done
test_that("writing varm works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  varm_x <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  zarr$varm <- list(PCs = varm_x)
  # varm should now have rownames added on-the-fly
  expected_varm_x <- varm_x
  rownames(expected_varm_x) <- zarr$var_names
  expect_identical(zarr$varm$PCs, expected_varm_x)
})

# trackstatus: class=ZarrAnnData, feature=test_set_obsp, status=done
test_that("writing obsp works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)

  obsp_x <- matrix(rnorm(10 * 10), nrow = 10, ncol = 10)
  zarr$obsp <- list(connectivities = obsp_x)
  # obsp should now have dimnames added on-the-fly
  expected_obsp_x <- obsp_x
  dimnames(expected_obsp_x) <- list(zarr$obs_names, zarr$obs_names)
  expect_identical(zarr$obsp$connectivities, expected_obsp_x)
})

# trackstatus: class=ZarrAnnData, feature=test_set_varp, status=done
test_that("writing varp works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)
  varp_x <- matrix(rnorm(20 * 20), nrow = 20, ncol = 20)
  zarr$varp <- list(connectivities = varp_x)
  # varp should now have dimnames added on-the-fly
  expected_varp_x <- varp_x
  dimnames(expected_varp_x) <- list(zarr$var_names, zarr$var_names)
  expect_identical(zarr$varp$connectivities, expected_varp_x)
})

# trackstatus: class=ZarrAnnData, feature=test_set_uns, status=done
test_that("writing uns works", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store = store)
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  zarr <- ZarrAnnData$new(store, obs = obs, var = var)
  zarr$uns <- list(
    foo = "bar",
    baz = c(1, 2, 3),
    nested = list(
      nested_foo = "nested_bar",
      nested_baz = c(4L, 5L, 6L)
    )
  )
  expect_identical(zarr$uns$foo, "bar")
  expect_equal(zarr$uns$baz, c(1, 2, 3), ignore_attr = TRUE)
  expect_identical(zarr$uns$nested$nested_foo, "nested_bar")
  expect_equal(zarr$uns$nested$nested_baz, c(4L, 5L, 6L), ignore_attr = TRUE)
})

# ERROR HANDLING ---------------------------------------------------------
test_that("opening a non-existent path in read mode errors", {
  expect_error(
    ZarrAnnData$new(tempfile(fileext = ".zarr"), mode = "r"),
    "does not exist"
  )
})

test_that("opening a non-existent path in r+ mode errors", {
  expect_error(
    ZarrAnnData$new(tempfile(fileext = ".zarr"), mode = "r+"),
    "does not exist"
  )
})

test_that("opening an existing file in exclusive-create mode errors", {
  store <- tempfile(fileext = ".zarr")
  create_zarr(store)
  on.exit(unlink(store, recursive = TRUE))
  expect_error(
    ZarrAnnData$new(store, mode = "w-"),
    "already exists"
  )
})

test_that("writing to a read-only store errors", {
  store <- tempfile(fileext = ".zarr")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  ZarrAnnData$new(store, obs = obs, var = var)
  on.exit(unlink(store, recursive = TRUE))

  zarr_ro <- ZarrAnnData$new(store, mode = "r")
  expect_error(
    zarr_ro$X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20),
    "read-only"
  )
})

# CONVERSION -------------------------------------------------------------
test_that("as_ZarrAnnData() round-trip from InMemoryAnnData works", {
  mem <- AnnData(
    X = matrix(1:20, nrow = 4, ncol = 5),
    obs = data.frame(a = 1:4, row.names = paste0("obs", 1:4)),
    var = data.frame(b = 1:5, row.names = paste0("var", 1:5)),
    layers = list(counts = matrix(21:40, nrow = 4, ncol = 5)),
    uns = list(foo = "bar")
  )

  store <- tempfile(fileext = ".zarr")
  on.exit(unlink(store, recursive = TRUE))

  zarr <- as_ZarrAnnData(mem, file = store)
  expect_true(inherits(zarr, "ZarrAnnData"))
  expect_equal(zarr$obs_names, mem$obs_names)
  expect_equal(zarr$var_names, mem$var_names)
  expect_equal(as.matrix(zarr$X), mem$X, ignore_attr = TRUE)
  expect_equal(
    as.matrix(zarr$layers$counts),
    mem$layers$counts,
    ignore_attr = TRUE
  )
  expect_equal(zarr$uns$foo, mem$uns$foo)
})
