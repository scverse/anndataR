skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

test_names <- names(data$obs)

# TODO: re-enable tests
test_names <- test_names[!grepl("_with_nas", test_names)]

for (name in test_names) {
  test_that(paste0("roundtrip with obs and var '", name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$X,
      obs = data$obs[, name, drop = FALSE],
      var = data$var[, name, drop = FALSE]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$obs[[name]],
      data$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$var[[name]],
      data$var[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in test_names) {
  test_that(paste0("reticulate->hdf5 with obs and var '", name, "'"), {
    ad <- anndata::AnnData(
      obs = data$obs[, name, drop = FALSE],
      var = data$var[, name, drop = FALSE]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obs[[name]],
      data$obs[[name]],
      tolerance = 1e-6
    )
  })
}

for (name in test_names) {
  test_that(paste0("hdf5->reticulate with obs and var '", name, "'"), {
    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))

    # strip rownames
    obs <- data$obs[, name, drop = FALSE]
    var <- data$var[, name, drop = FALSE]

    # create anndata
    ad <- AnnData(
      obs = obs,
      var = var
    )
    write_h5ad(ad, filename)

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    obs_ <- ad_new$obs[[name]]
    expect_equal(
      obs_,
      data$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
