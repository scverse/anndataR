test_that("obsm/varm validation works", {
  N_OBS <- 5
  N_VAR <- 3

  mtx <- matrix(
    0,
    N_OBS,
    N_VAR
  )

  adata <- AnnData(
    X = mtx,
    obs = data.frame(row.names = as.character(1:N_OBS)),
    var = data.frame(row.names = as.character(1:N_VAR))
  )

  adata$obsm <- list(PCA = matrix(0, N_OBS, 4))
  adata$varm <- list(PCs = matrix(0, N_VAR, 4))

  expect_error(adata$obsm <- list(PCA = matrix(0, 4, 4)))
  expect_error(adata$varm <- list(PCs = matrix(0, 4, 4)))
})

test_that("obsp/varp validation works", {
  N_OBS <- 5
  N_VAR <- 3

  adata <- AnnData(
    obs = data.frame(row.names = as.character(1:N_OBS)),
    var = data.frame(row.names = as.character(1:N_VAR))
  )

  adata$obsp <- list(graph1 = matrix(0, N_OBS, N_OBS))
  adata$varp <- list(graph1 = matrix(0, N_VAR, N_VAR))

  expect_error(adata$obsp <- list(graph1 = matrix(0, 4, 4)))
  expect_error(adata$varp <- list(graph1 = matrix(0, 4, 4)))
})
