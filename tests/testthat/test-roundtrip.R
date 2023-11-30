library(testthat)
library(Matrix)

h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")

base_file <- system.file("extdata", "example.h5ad", package = "anndataR")

gen_adata <- function(type) {
  N_OBS <- 10
  N_VAR <- 15
  obs_names <- paste0("obs_", seq_len(N_OBS))
  var_names <- paste0("var_", seq_len(N_VAR))
  adata <- AnnData(
    X = rsparsematrix(N_OBS, N_VAR, 0.1),
    obs_names = obs_names,
    var_names = var_names,
    layers = list(
      dense = matrix(seq_len(N_OBS * N_VAR), N_OBS, N_VAR),
      sparse = rsparsematrix(N_OBS, N_VAR, 0.1)
    ),
    obsm = list(
      dense = matrix(seq_len(N_OBS * 5), N_OBS, 5),
      sparse = rsparsematrix(N_OBS, 5, 0.1)
    ),
    varm = list(
      dense = matrix(seq_len(N_VAR * 5), N_VAR, 5),
      sparse = rsparsematrix(N_VAR, 5, 0.1)
    ),
    obsp = list(
      dense = matrix(seq_len(N_OBS * N_OBS), N_OBS, N_OBS),
      sparse = rsparsematrix(N_OBS, N_OBS, 0.1)
    ),
    varp = list(
      dense = matrix(seq_len(N_VAR * N_VAR), N_VAR, N_VAR),
      sparse = rsparsematrix(N_VAR, N_VAR, 0.1)
    )
  )
  if (type == "HDF5AnnData") {
    tempfile(pattern = "hdf5_write_", fileext = ".h5ad")
    write_h5ad(adata, h5ad_file)
    read_h5ad(h5ad_file, to = type)
  } else if (type == "InMemoryAnnData") {
    adata
  } else {
    stop(paste0("Unknown type: ", type))
  }
}

check_round_trip <- function(expected, type) {
  h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")
  write_h5ad(expected, h5ad_file)
  actual <- read_h5ad(h5ad_file, to = type)

  expect_equal(actual, expected)
}

for (typ in c("HDF5AnnData", "InMemoryAnnData")) {
  test_that(paste("round trip w/ example data for", typ), {
    adata <- read_h5ad(base_file, to = typ)
    check_round_trip(adata, typ)
  })
  test_that(paste("round trip w/ generated data for", typ), {
    adata <- gen_adata(typ)
    check_round_trip(adata, typ)
  })
}
