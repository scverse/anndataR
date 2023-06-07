skip_if_not_installed("rhdf5")

file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("opening H5AD works", {
  adata <- HDF5AnnData$new(file)
  expect_true(is(adata, "HDF5AnnData"))
})

adata <- HDF5AnnData$new(file)

# GETTERS ----------------------------------------------------------------
# trackstatus: class=HDF5AnnData, feature=test_get_X, status=wip
test_that("reading X works", {
  X <- adata$X
  expect_s4_class(X, "dgRMatrix")
  expect_equal(dim(X), c(50, 100))
})

# trackstatus: class=HDF5AnnData, feature=test_get_layers, status=wip
test_that("reading layers works", {
  layers <- adata$layers
  expect_true(is.list(layers), "list")
  expect_equal(
    names(layers),
    c("counts", "csc_counts", "dense_X", "dense_counts")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obs, status=wip
test_that("reading obs works", {
  obs <- adata$obs
  expect_s3_class(obs, "data.frame")
  expect_equal(
    colnames(obs),
    c(
      "Float", "FloatNA", "Int", "IntNA", "Bool", "BoolNA", "n_genes_by_counts",
      "log1p_n_genes_by_counts", "total_counts", "log1p_total_counts", "leiden"
    )
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_var, status=wip
test_that("reading var works", {
  var <- adata$var
  expect_s3_class(var, "data.frame")
  expect_equal(
    colnames(var),
    c(
      "String", "n_cells_by_counts", "mean_counts", "log1p_mean_counts",
      "pct_dropout_by_counts", "total_counts", "log1p_total_counts",
      "highly_variable", "means", "dispersions", "dispersions_norm"
    )
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obs_names, status=wip
test_that("reading obs names works", {
  obs_names <- adata$obs_names
  expect_vector(obs_names, ptype = character(), size = 50)
})

# trackstatus: class=HDF5AnnData, feature=test_get_var_names, status=wip
test_that("reading var names works", {
  var_names <- adata$var_names
  expect_vector(var_names, ptype = character(), size = 100)
})

# SETTERS ----------------------------------------------------------------
