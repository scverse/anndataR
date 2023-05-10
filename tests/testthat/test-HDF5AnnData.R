library(rhdf5)

file_path <- system.file("extdata", "krumsiek11_augmented_sparse_v0-8.h5ad", package = "anndataR")

# >>> ad.read_h5ad("inst/extdata/krumsiek11_augmented_v0-8.h5ad")
# AnnData object with n_obs × n_vars = 640 × 11
#     obs: 'cell_type', 'dummy_num', 'dummy_num2', 'dummy_int', 'dummy_int2', 'dummy_bool', 'dummy_bool2'
#     var: 'dummy_str'
#     uns: 'dummy_bool', 'dummy_bool2', 'dummy_category', 'dummy_int', 'dummy_int2', 'highlights', 'iroot'


# GETTERS ----------------------------------------------------------------
# trackstatus: class=HDF5AnnData, feature=test_get_X, status=wip
test_that("read X", {
  h5_file <- H5Fopen(file_path)

  tryCatch({
    adata <- HDF5AnnData$new(h5_file)

    X <- adata$X
    expect_equal(nrow(X), 640L)
    expect_equal(ncol(X), 11L)

    # todo: check content of X
  }, finally = {
    H5Fclose(h5_file)
  })
})

# trackstatus: class=HDF5AnnData, feature=test_get_obs, status=wip
test_that("read obs", {
  h5_file <- H5Fopen(file_path)

  tryCatch({
    adata <- HDF5AnnData$new(h5_file)

    obs <- adata$obs
    expect_equal(nrow(obs), 640L)

    # todo: check content of obs
  }, finally = {
    H5Fclose(h5_file)
  })
})

# trackstatus: class=HDF5AnnData, feature=test_get_var, status=wip
test_that("read var", {
  h5_file <- H5Fopen(file_path)

  tryCatch({
    adata <- HDF5AnnData$new(h5_file)

    var <- adata$var
    expect_equal(nrow(var), 11L)

    # todo: check content of var
  }, finally = {
    H5Fclose(h5_file)
  })
})

# SETTERS ----------------------------------------------------------------
