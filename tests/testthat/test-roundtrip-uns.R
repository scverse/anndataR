skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

uns_names <- names(data$uns)
# TODO: re-enable these tests
uns_names <- uns_names[!grepl("_with_nas", uns_names)]
# TODO: re-enable these tests
uns_names <- uns_names[!grepl("_na$", uns_names)]
# TODO: re-enable these tests
uns_names <- uns_names[!grepl("mat_", uns_names)]
# TODO: re-enable these tests
uns_names <- uns_names[!uns_names %in% c("vec_factor", "vec_factor_ordered", "vec_logical")]
# TODO: re-enable these tests
uns_names <- uns_names[uns_names != "list"]
# TODO: re-enable these tests
uns_names <- uns_names[!uns_names %in% c("scalar_factor", "scalar_factor_ordered", "scalar_logical")]
# TODO: re-enable these tests
uns_names <- uns_names[!uns_names %in% c("df_factor", "df_factor_ordered", "df_logical")]

for (name in uns_names) {
  test_that(paste0("roundtrip with uns '", name, "'"), {
    # create anndata
    ad <- AnnData(
      obs = data.frame(row.names = rownames(data$obs)),
      var = data.frame(row.names = rownames(data$var)),
      uns = data$uns[name]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
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
      obs = data.frame(row.names = rownames(data$obs)),
      var = data.frame(row.names = rownames(data$var)),
      uns = data$uns[name]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
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
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))

    # make anndata
    ad <- AnnData(
      obs = data.frame(row.names = rownames(data$obs)),
      var = data.frame(row.names = rownames(data$var)),
      uns = data$uns[name]
    )
    write_h5ad(ad, filename)

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
