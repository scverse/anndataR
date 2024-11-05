skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

test_names <- names(data$layers)

# TODO: Add denseMatrix support to anndata and anndataR
test_names <- test_names[!grepl("_dense", test_names)]

for (name in test_names) {
  test_that(paste0("roundtrip with layer '", name, "'"), {
    # create anndata
    ad <- AnnData(
      layers = data$layers[name],
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE]
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$layers[[name]],
      data$layers[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in test_names) {
  test_that(paste0("reticulate->hdf5 with layer '", name, "'"), {
    # add rownames
    layers <- data$layers[name]
    obs <- data.frame(row.names = rownames(data$obs))
    var <- data.frame(row.names = rownames(data$var))

    # create anndata
    ad <- anndata::AnnData(
      layers = layers,
      shape = dim(data$X),
      obs = obs,
      var = var
    )

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$layers[[name]],
      data$layers[[name]],
      tolerance = 1e-6
    )
  })
}

r2py_names <- test_names
# TODO: rsparse gets converted to csparse by anndata
r2py_names <- r2py_names[!grepl("rsparse", r2py_names)]

for (name in r2py_names) {
  test_that(paste0("hdf5->reticulate with layer '", name, "'"), {
    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))

    # make anndata
    ad <- AnnData(
      layers = data$layers[name],
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE]
    )
    write_h5ad(ad, filename)

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    layer_ <- ad_new$layers[[name]]
    # anndata returns these layers as CsparseMatrix
    if (grepl("rsparse", name)) {
      layer_ <- as(layer_, "RsparseMatrix")
    }
    # strip rownames
    dimnames(layer_) <- list(NULL, NULL)
    expect_equal(
      layer_,
      data$layers[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
