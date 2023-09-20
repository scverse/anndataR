skip_if_no_anndata()
skip_if_not_installed("rhdf5")

data <- generate_dataset_as_list(10L, 20L)

for (name in names(data$obsm)) {
  test_that(paste("roundtrip with obsm and varm '", name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$X,
      obsm = data$obsm[name],
      varm = data$varm[name],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # write to file
    filename <- withr::local_file(paste0("roundtrip_obsmvarm_", name, ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$obsm[[name]],
      data$obsm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    expect_equal(
      ad_new$varm[[name]],
      data$varm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obsm)) {
  test_that(paste0("reticulate->hdf5 with obsm and varm '", name, "'"), {
    ad <- anndata::AnnData(
      X = data$X,
      obsm = data$obsm[name],
      varm = data$varm[name]
    )
    ad$obs_names <- data$obs_names
    ad$var_names <- data$var_names

    # write to file
    filename <- withr::local_file(paste0("reticulate_to_hdf5_obsmvarm_", name, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obsm[[name]],
      data$obsm[[name]],
      tolerance = 1e-10
    )
    expect_equal(
      ad_new$varm[[name]],
      data$varm[[name]],
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obsm)) {
  test_that(paste0("hdf5->reticulate with obsm and varm '", name, "'"), {
    # write to file
    filename <- withr::local_file(paste0("hdf5_to_reticulate_obsmvarm_", name, ".h5ad"))
    ad <- HDF5AnnData$new(
      file = filename,
      X = data$X,
      obsm = data$obsm[name],
      varm = data$varm[name],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    obsm_ <- ad_new$obsm[[name]]
    dimnames(obsm_) <- list(NULL, NULL)
    expect_equal(
      obsm_,
      data$obsm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    varm_ <- ad_new$varm[[name]]
    dimnames(varm_) <- list(NULL, NULL)
    expect_equal(
      varm_,
      data$varm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}
