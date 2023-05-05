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

layers <- list(
  X2 = X * 2,
  X3 = X * 3
)

test_that("test Python -> h5ad -> R", {
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
  ad$write_h5ad("file.h5ad")

  # read from file
  ad_new <- HDF5AnnData$new("file.h5ad")

  # check whether ad_new$obs == obs and so on
})

test_that("test R -> h5ad -> Python", {
  # create anndata in R
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  # write to file
  ad$write_h5ad("file.h5ad")

  # read from file, hopefully Python anndata is able to read this
  ad_new <- anndata::read_h5ad("file.h5ad")

  # check whether ad_new$obs == obs and so on
})