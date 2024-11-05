skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

test_names <- names(data$obsm)

# TODO: re-enable this
test_names <- test_names[test_names != "character_with_nas"]

# TODO: Add denseMatrix support to anndata and anndataR
test_names <- test_names[!grepl("_dense", test_names)]

for (name in test_names) {
  test_that(paste0("roundtrip with obsm and varm '", name, "'"), {
    # create anndata
    ad <- AnnData(
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE],
      obsm = data$obsm[name],
      varm = data$varm[name]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$obsm[[name]],
      data$obsm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$varm[[name]],
      data$varm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

# TODO: re-enable these tests
# it seemed like there is a difference in the dimnames of the
# obsm and varm between anndata and anndataR
test_names <- c()

for (name in test_names) {
  test_that(paste0("reticulate->hdf5 with obsm and varm '", name, "'"), {
    # create anndata
    ad <- anndata::AnnData(
      obs = data.frame(row.names = data$obs_names),
      var = data.frame(row.names = data$var_names),
      obsm = data$obsm[name],
      varm = data$varm[name]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obsm[[name]],
      data$obsm[[name]],
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$varm[[name]],
      data$varm[[name]],
      tolerance = 1e-6
    )
  })
}

for (name in test_names) {
  test_that(paste0("hdf5->reticulate with obsm and varm '", name, "'"), {
    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))

    # create anndata
    ad <- AnnData(
      obsm = data$obsm[name],
      varm = data$varm[name],
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE]
    )
    write_h5ad(ad, filename)

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obsm[[name]],
      data$obsm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
    expect_equal(
      ad_new$varm[[name]],
      data$varm[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
