skip_if_not_installed("rhdf5")

obs_names <- paste0("cell", 1:3)
var_names <- paste0("gene", 1:4)

df_adata <- function() {
  AnnData(
    X = matrix(1:12, nrow = 3),
    obs = data.frame(row.names = obs_names),
    var = data.frame(row.names = var_names),
    obsm = list(df = data.frame(a = 1:3, b = c("x", NA, "z"))),
    varm = list(df = data.frame(a = c(0.5, NA, 1.5, 2), b = letters[1:4]))
  )
}

test_that("obsm and varm data frames are written with the obs and var names as index", {
  file <- withr::local_file(tempfile(fileext = ".h5ad"))
  write_h5ad(df_adata(), file, mode = "w")

  expect_identical(
    as.vector(rhdf5::h5read(file, "obsm/df/_index")),
    obs_names
  )
  expect_identical(
    as.vector(rhdf5::h5read(file, "varm/df/_index")),
    var_names
  )
})

test_that("obsm and varm data frames round-trip with their values", {
  file <- withr::local_file(tempfile(fileext = ".h5ad"))
  write_h5ad(df_adata(), file, mode = "w")
  adata <- read_h5ad(file, as = "InMemoryAnnData")

  expect_identical(rownames(adata$varm$df), var_names)
  expect_identical(adata$varm$df$a, c(0.5, NA, 1.5, 2))
  expect_identical(adata$varm$df$b, letters[1:4])
  expect_identical(rownames(adata$obsm$df), obs_names)
  expect_identical(adata$obsm$df$b, c("x", NA, "z"))
})

test_that("Python anndata reads obsm and varm data frames written by R", {
  skip_if_no_anndata_py()
  file <- withr::local_file(tempfile(fileext = ".h5ad"))
  write_h5ad(df_adata(), file, mode = "w")

  ad <- reticulate::import("anndata", convert = FALSE)
  adata_py <- ad$read_h5ad(file)
  varm_index <- reticulate::py_to_r(adata_py$varm[["df"]]$index$tolist())
  obsm_index <- reticulate::py_to_r(adata_py$obsm[["df"]]$index$tolist())

  expect_identical(unlist(varm_index), var_names)
  expect_identical(unlist(obsm_index), obs_names)
})
