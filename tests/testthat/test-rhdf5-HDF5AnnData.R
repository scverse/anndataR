skip_if_not_installed("rhdf5")

requireNamespace("vctrs")

file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("opening H5AD works", {
  adata <- HDF5AnnData$new(file, rhdf5 = TRUE)
  expect_true(inherits(adata, "HDF5AnnData"))
  adata$close()
})

adata <- HDF5AnnData$new(file, rhdf5 = TRUE)

# GETTERS ----------------------------------------------------------------
# trackstatus: class=HDF5AnnData, feature=test_get_rhdf5_X, status=done
test_that("reading X works", {
  X <- adata$X
  expect_s4_class(X, "dgRMatrix")
  expect_equal(dim(X), c(50, 100))
})

# trackstatus: class=HDF5AnnData, feature=test_get_rhdf5_layers, status=done
test_that("reading layers works", {
  layers <- adata$layers
  expect_true(is.list(layers), "list")
  expect_equal(
    names(layers),
    c("counts", "csc_counts", "dense_X", "dense_counts")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_obsm, status=done
test_that("reading obsm works", {
  obsm <- adata$obsm
  expect_true(is.list(obsm), "list")
  expect_equal(
    names(obsm),
    c("X_pca", "X_umap")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_varm, status=done
test_that("reading varm works", {
  varm <- adata$varm
  expect_true(is.list(varm), "list")
  expect_equal(
    names(varm),
    c("PCs")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_obs, status=done
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

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_var, status=done
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

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_obs_names, status=done
test_that("reading obs names works", {
  obs_names <- adata$obs_names
  expect_vector(obs_names, ptype = character(), size = 50)
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_var_names, status=done
test_that("reading var names works", {
  var_names <- adata$var_names
  expect_vector(var_names, ptype = character(), size = 100)
})

# SETTERS ----------------------------------------------------------------
test_that("creating empty H5AD works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  expect_silent(HDF5AnnData$new(file = h5ad_file, rhdf5 = TRUE))
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_X, status=done
test_that("writing X works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var, rhdf5 = TRUE)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(h5ad$X <- X)
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_layers, status=done
test_that("writing layers works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var, rhdf5 = TRUE)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(h5ad$layers <- list(layer1 = X, layer2 = X))
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_obs, status=done
test_that("writing obs works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var, rhdf5 = TRUE)

  obs <- data.frame(
    Letters = LETTERS[1:10],
    Numbers = 1:10,
    row.names = paste0("Row", 1:10)
  )
  h5ad$obs <- obs
  expect_identical(h5ad$obs_names, paste0("Row", 1:10))
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_var, status=done
test_that("writing var works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var, rhdf5 = TRUE)

  var <- data.frame(
    Letters = LETTERS[1:20],
    Numbers = 1:20,
    row.names = paste0("Row", 1:20)
  )
  h5ad$var <- var
  expect_identical(h5ad$var_names, paste0("Row", 1:20))
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_obs_names, status=done
test_that("writing obs names works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var, rhdf5 = TRUE)

  h5ad$obs_names <- LETTERS[1:10]
  expect_identical(h5ad$obs_names, LETTERS[1:10])
})

# trackstatus: class=HDF5AnnData, feature=test_set_rhdf5_var_names, status=done
test_that("writing var names works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var, rhdf5 = TRUE)

  h5ad$var_names <- LETTERS[1:20]
  expect_identical(h5ad$var_names, LETTERS[1:20])
})

test_that("deprecated to_InMemoryAnnData() works", {
  expect_warning(mem_ad <- adata$to_InMemoryAnnData())
  expect_true(inherits(mem_ad, "InMemoryAnnData"))
})
