skip_if_not_installed("rhdf5")

requireNamespace("vctrs")

file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("opening H5AD works", {
  adata <- HDF5AnnData$new(file, mode = "r")
  expect_true(inherits(adata, "HDF5AnnData"))
  adata$close()
})

adata <- HDF5AnnData$new(file, mode = "r")

# GETTERS ----------------------------------------------------------------
# trackstatus: class=HDF5AnnData, feature=test_get_X, status=done
test_that("reading X works", {
  X <- adata$X
  expect_s4_class(X, "dgRMatrix")
  expect_equal(dim(X), c(50, 100))
})

# trackstatus: class=HDF5AnnData, feature=test_get_layers, status=done
test_that("reading layers works", {
  layers <- adata$layers
  expect_true(is.list(layers), "list")
  expect_equal(
    names(layers),
    c("counts", "csc_counts", "dense_X", "dense_counts")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obsm, status=done
test_that("reading obsm works", {
  obsm <- adata$obsm
  expect_true(is.list(obsm), "list")
  expect_equal(
    names(obsm),
    c("X_pca", "X_umap")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_varm, status=done
test_that("reading varm works", {
  varm <- adata$varm
  expect_true(is.list(varm), "list")
  expect_equal(
    names(varm),
    c("PCs")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obsp, status=done
test_that("reading obsp works", {
  obsp <- adata$obsp
  expect_true(is.list(obsp), "list")
  expect_equal(
    names(obsp),
    c("connectivities", "distances")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_varp, status=done
test_that("reading varp works", {
  varp <- adata$varp
  expect_true(is.list(varp), "list")
  expect_equal(
    names(varp),
    c("test_varp")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obs, status=done
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

# trackstatus: class=HDF5AnnData, feature=test_get_var, status=done
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

# trackstatus: class=HDF5AnnData, feature=test_get_obs_names, status=done
test_that("reading obs names works", {
  obs_names <- adata$obs_names
  expect_vector(obs_names, ptype = character(), size = 50)
})

# trackstatus: class=HDF5AnnData, feature=test_get_var_names, status=done
test_that("reading var names works", {
  var_names <- adata$var_names
  expect_vector(var_names, ptype = character(), size = 100)
})

# SETTERS ----------------------------------------------------------------
test_that("creating empty H5AD works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  expect_silent(HDF5AnnData$new(file = h5ad_file))
})

# trackstatus: class=HDF5AnnData, feature=test_set_X, status=done
test_that("writing X works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(h5ad$X <- X)
})

# trackstatus: class=HDF5AnnData, feature=test_set_layers, status=done
test_that("writing layers works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(h5ad$layers <- list(layer1 = X, layer2 = X))
})

# trackstatus: class=HDF5AnnData, feature=test_set_obs, status=done
test_that("writing obs works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  obs <- data.frame(
    Letters = LETTERS[1:10],
    Numbers = 1:10,
    row.names = paste0("Row", 1:10)
  )
  h5ad$obs <- obs
  expect_identical(h5ad$obs_names, paste0("Row", 1:10))
})

# trackstatus: class=HDF5AnnData, feature=test_set_var, status=done
test_that("writing var works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  var <- data.frame(
    Letters = LETTERS[1:20],
    Numbers = 1:20,
    row.names = paste0("Row", 1:20)
  )
  h5ad$var <- var
  expect_identical(h5ad$var_names, paste0("Row", 1:20))
})

# trackstatus: class=HDF5AnnData, feature=test_set_obs_names, status=done
test_that("writing obs names works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  h5ad$obs_names <- LETTERS[1:10]
  expect_identical(h5ad$obs_names, LETTERS[1:10])
})

# trackstatus: class=HDF5AnnData, feature=test_set_var_names, status=done
test_that("writing var names works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  h5ad$var_names <- LETTERS[1:20]
  expect_identical(h5ad$var_names, LETTERS[1:20])
})

# trackstatus: class=HDF5AnnData, feature=test_set_obsm, status=done
test_that("writing obsm works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  obsm_x <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  h5ad$obsm <- list(X = obsm_x)
  expect_identical(h5ad$obsm$X, obsm_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_varm, status=done
test_that("writing varm works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)
  varm_x <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  h5ad$varm <- list(PCs = varm_x)
  expect_identical(h5ad$varm$PCs, varm_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_obsp, status=done
test_that("writing obsp works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  obsp_x <- matrix(rnorm(10 * 10), nrow = 10, ncol = 10)
  h5ad$obsp <- list(connectivities = obsp_x)
  expect_identical(h5ad$obsp$connectivities, obsp_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_varp, status=done
test_that("writing varp works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  varp_x <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  h5ad$varp <- list(connectivities = varp_x)
  expect_identical(h5ad$varp$connectivities, varp_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_uns, status=done
test_that("writing uns works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  h5ad$uns <- list(
    foo = "bar",
    baz = c(1, 2, 3),
    nested = list(
      nested_foo = "nested_bar",
      nested_baz = c(4L, 5L, 6L)
    )
  )
  expect_identical(h5ad$uns$foo, "bar")
  expect_identical(h5ad$uns$baz, c(1, 2, 3))
  expect_identical(h5ad$uns$nested$nested_foo, "nested_bar")
  expect_identical(h5ad$uns$nested$nested_baz, c(4L, 5L, 6L))
})
