library(reticulate)

X <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = .1)
obs <- data.frame(
  cell_type = sample(c("tcell", "bcell"), 10, replace = TRUE),
  cluster = sample.int(3, 10, replace = TRUE)
)
obs_names <- paste0("cell_", seq_len(10))
var <- data.frame(
  geneinfo = sample(c("a", "b", "c"), 20, replace = TRUE)
)
var_names <- paste0("gene", seq_len(20))
layers <- list(
  X2 = matrix(2L, nrow = 10, ncol = 20),
  X3 = matrix(2.5, nrow = 10, ncol = 20),
  X4 = Matrix::rsparsematrix(10, 20, .2)
)
admem <- InMemoryAnnData$new(
  X = X,
  obs = obs,
  var = var,
  obs_names = obs_names,
  var_names = var_names,
  layers = layers
)

# GETTERS ----------------------------------------------------------------
test_that("read from file, test getters", {
  withr::local_file("foo.h5ad")

  to_reticulate(admem)$write_h5ad("foo.h5ad")

  adata <- ReticulateAnnData$read_h5ad("foo.h5ad")

  expect_equal(adata$X, X)
  expect_equal(adata$obs, obs, ignore_attr = TRUE)
  expect_equal(adata$var, var, ignore_attr = TRUE)
  expect_equal(adata$obs_names, obs_names, ignore_attr = TRUE)
  expect_equal(adata$var_names, var_names, ignore_attr = TRUE)
  expect_equal(adata$layers, layers, ignore_attr = TRUE)
})


test_that("convert from inmem, test getters", {
  adata <- to_reticulate(admem)

  expect_equal(adata$X, X)
  expect_equal(adata$obs, obs)
  expect_equal(adata$var, var)
  expect_identical(adata$shape(), c(10L, 20L))
  expect_equal(adata$obs_names, obs_names)
  expect_equal(adata$var_names, var_names)
  expect_equal(adata$layers, layers)
})


# SETTERS ----------------------------------------------------------------