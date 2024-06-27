skip_if_no_anndata()
skip_if_not_installed("rhdf5")

data <- generate_dataset(10L, 20L)

uns_names <- names(data$uns)
# TODO: re-enable these tests
uns_names <- uns_names[!grepl("_with_nas", uns_names)]
# TODO: re-enable these tests
uns_names <- uns_names[!grepl("_na$", uns_names)]
# TODO: re-enable these tests
uns_names <- uns_names[!grepl("mat_", uns_names)]
# # TODO: re-enable these tests
# uns_names <- uns_names[!uns_names %in% c("vec_factor", "vec_factor_ordered", "vec_logical")] #nolint
# TODO: re-enable these tests
uns_names <- uns_names[uns_names != "list"]
# TODO: re-enable these tests
# uns_names <- uns_names[!uns_names %in% c("scalar_factor", "scalar_factor_ordered", "scalar_logical")] # nolint
# TODO: re-enable these tests
# uns_names <- uns_names[!uns_names %in% c("df_factor", "df_factor_ordered", "df_logical")] # nolint

for (name in uns_names) {
  test_that(paste0("roundtrip with uns '", name, "'"), {
    # create anndata
    ad <- AnnData(
      obs_names = data$obs_names,
      var_names = data$var_names,
      uns = data$uns[name]
    )

    # write to file
    filename <- withr::local_file(paste0("roundtrip_uns_", name, ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$uns,
      data$uns[name],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in uns_names) {
  test_that(paste0("reticulate->hdf5 with uns '", name, "'"), {
    # create anndata
    ad <- anndata::AnnData(
      obs = data.frame(row.names = data$obs_names),
      var = data.frame(row.names = data$var_names),
      uns = data$uns[name]
    )

    # write to file
    filename <- withr::local_file(paste0("reticulate_to_hdf5_uns_", name, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$uns,
      data$uns[name],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in uns_names) {
  test_that(paste0("hdf5->reticulate with uns '", name, "'"), {
    # write to file
    filename <- withr::local_file(paste0("hdf5_to_reticulate_uns_", name, ".h5ad"))

    # make anndata
    ad <- HDF5AnnData$new(
      file = filename,
      obs_names = data$obs_names,
      var_names = data$var_names,
      uns = data$uns[name]
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    expect_equal(
      ad_new$uns,
      data$uns[name],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
