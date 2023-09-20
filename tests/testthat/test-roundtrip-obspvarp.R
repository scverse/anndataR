skip_if_no_anndata()
skip_if_not_installed("rhdf5")

data <- generate_dataset_as_list(10L, 20L)

for (name in names(data$obsp)) {
  test_that(paste0("roundtrip with obsp and varp '", name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$X,
      obsp = data$obsp[name],
      varp = data$varp[name],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # write to file
    filename <- withr::local_file(paste0("roundtrip_obspvarp_", name, ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$obsp[[name]],
      data$obsp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    expect_equal(
      ad_new$varp[[name]],
      data$varp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obsp)) {
  test_that(paste0("reticulate->hdf5 with obsp and varp '", name, "'"), {
    ad <- anndata::AnnData(
      X = data$X,
      obsp = data$obsp[name],
      varp = data$varp[name]
    )
    ad$obs_names <- data$obs_names
    ad$var_names <- data$var_names

    # write to file
    filename <- withr::local_file(paste0("reticulate_to_hdf5_obspvarp_", name, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obsp[[name]],
      data$obsp[[name]],
      tolerance = 1e-10
    )
    expect_equal(
      ad_new$varp[[name]],
      data$varp[[name]],
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obsp)) {
  test_that(paste0("hdf5->reticulate with obsp and varp '", name, "'"), {
    # write to file
    filename <- withr::local_file(paste0("hdf5_to_reticulate_obspvarp_", name, ".h5ad"))
    ad <- HDF5AnnData$new(
      file = filename,
      X = data$X,
      obsp = data$obsp[name],
      varp = data$varp[name],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    obsp_ <- ad_new$obsp[[name]]
    dimnames(obsp_) <- list(NULL, NULL)
    expect_equal(
      obsp_,
      data$obsp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    varp_ <- ad_new$varp[[name]]
    dimnames(varp_) <- list(NULL, NULL)
    expect_equal(
      varp_,
      data$varp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}
