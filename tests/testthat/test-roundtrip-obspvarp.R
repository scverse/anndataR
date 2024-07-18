skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

test_names <- names(data$obsp)

# TODO: re-enable test
test_names <- c()

for (name in test_names) {
  test_that(paste0("roundtrip with obsp and varp '", name, "'"), {
    # create anndata
    ad <- AnnData(
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE],
      obsp = data$obsp[name],
      varp = data$varp[name]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$obsp[[name]],
      data$obsp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$varp[[name]],
      data$varp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in test_names) {
  test_that(paste0("reticulate->hdf5 with obsp and varp '", name, "'"), {
    # create anndata
    ad <- anndata::AnnData(
      obs = data.frame(row.names = data$obs_names),
      var = data.frame(row.names = data$var_names),
      obsp = data$obsp[name],
      varp = data$varp[name]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obsp[[name]],
      data$obsp[[name]],
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$varp[[name]],
      data$varp[[name]],
      tolerance = 1e-6
    )
  })
}

for (name in test_names) {
  test_that(paste0("hdf5->reticulate with obsp and varp '", name, "'"), {
    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))

    # create anndata
    ad <- AnnData(
      obsp = data$obsp[name],
      varp = data$varp[name],
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE]
    )
    write_h5ad(ad, filename)

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obsp[[name]],
      data$obsp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$varp[[name]],
      data$varp[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
