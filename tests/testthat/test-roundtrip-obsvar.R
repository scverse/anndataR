skip_if_no_anndata()
skip_if_not_installed("rhdf5")

data <- generate_dataset_as_list(10L, 20L)

for (name in names(data$obs)) {
  test_that(paste("roundtrip with obs and var '", name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$X,
      obs = data$obs[, name, drop = FALSE],
      var = data$var[, name, drop = FALSE],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # write to file
    filename <- withr::local_file(paste0("roundtrip_obsvar_", name, ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$obs[[name]],
      data$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    expect_equal(
      ad_new$var[[name]],
      data$var[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obs)) {
  test_that(paste0("reticulate->hdf5 with obs and var '", name, "'"), {
    ad <- anndata::AnnData(
      obs = data$obs[, name, drop = FALSE],
      var = data$var[, name, drop = FALSE]
    )
    ad$obs_names <- data$obs_names
    ad$var_names <- data$var_names

    # write to file
    filename <- withr::local_file(paste0("reticulate_to_hdf5_obsvar_", name, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obs[[name]],
      data$obs[[name]],
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obs)) {
  test_that(paste0("hdf5->reticulate with obs and var '", name, "'"), {
    # write to file
    filename <- withr::local_file(paste0("hdf5_to_reticulate_obsvar_", name, ".h5ad"))

    # strip rownames
    obs <- data$obs[, name, drop = FALSE]
    var <- data$var[, name, drop = FALSE]
    rownames(obs[[name]]) <- NULL
    rownames(var[[name]]) <- NULL

    # create anndata
    ad <- HDF5AnnData$new(
      file = filename,
      obs = obs,
      var = var,
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    obs_ <- ad_new$obs[[name]]
    dimnames(obs_) <- list(NULL, NULL)
    expect_equal(
      obs_,
      data$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}
