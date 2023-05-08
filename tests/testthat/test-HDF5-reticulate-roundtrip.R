# helper function to skip tests if we don't have the 'foo' module
skip_if_no_anndata <- function() {
  requireNamespace("reticulate")
  have_anndata <- reticulate::py_module_available("anndata")
  if (!have_anndata) {
    skip("anndata not available for testing")
  }
}

skip_if_no_anndata()

# construct dummy objects
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
  X2 = X * 2,
  X3 = X * 3
)

test_that("test reticulate+Python -> h5ad -> R+HDF5AnnData", {
  filename <- withr::local_file("file.h5ad")

  # create anndata in python
  obs_ <- obs
  rownames(obs_) <- obs_names
  var_ <- var
  rownames(var_) <- var_names
  ad <- anndata::AnnData(
    X = X,
    obs = obs_,
    var = var_
  )

  # write to file
  ad$write_h5ad(filename)

  # read from file
  ad_new <- HDF5AnnData$new(filename)

  # check whether ad_new$obs == obs and so on
  expect_equal(ad_new$n_obs(), nrow(obs))
  expect_equal(ad_new$n_vars(), nrow(var))
})

# It's not possible to do this roundtrip yet
# nolint start
# test_that("test R+HDF5AnnData -> h5ad -> reticulate+Python", {
#   # create anndata in R
#   ad <- InMemoryAnnData$new(
#     X = X,
#     obs = obs,
#     var = var,
#     obs_names = obs_names,
#     var_names = var_names
#   )

#   # convert and write to file
#   ad$to_HDF5AnnData("file.h5ad")

#   # read from file, hopefully Python anndata is able to read this
#   ad_new <- anndata::read_h5ad("file.h5ad")

#   # check whether ad_new$obs == obs and so on
#   expect_equal(ad_new$n_obs(), nrow(obs))
#   expect_equal(ad_new$n_vars(), nrow(var))
# })
# nolint end
