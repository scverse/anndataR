library(Matrix)

# construct dummy X
X <- Matrix::rsparsematrix(nrow = 10, ncol = 20, density = .1)

# construct dummy obs
obs <- data.frame(
  cell_type = sample(c("tcell", "bcell"), 10, replace = TRUE),
  cluster = sample.int(3, 10, replace = TRUE)
)

# construct dummy obs names
obs_names <- paste0("cell", seq_len(10))

# construct dummy var
var <- data.frame(
  geneinfo = sample(c("a", "b", "c"), 20, replace = TRUE)
)

# construct dummy var names
var_names <- paste0("gene", seq_len(20))


test_that("to_Seurat with inmemoryanndata", {
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  seu <- to_Seurat(ad)

  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)
  expect_equal(rownames(seu), var_names)
  expect_equal(colnames(seu), obs_names)

  # todo: check whether X can be retrieved
})

# todo: check substitution from _ to -

test_that("to_Seurat() fails gracefully", {
  expect_error(to_Seurat(), regexp = "obj.*is missing")
  expect_error(to_Seurat("foo"), regexp = "AbstractAnnData.*not TRUE")
})
