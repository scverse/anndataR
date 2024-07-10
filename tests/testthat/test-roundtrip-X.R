skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

layer_names <- names(data$layers)
# TODO: Add denseMatrix support to anndata and anndataR
layer_names <- layer_names[!grepl("_dense", layer_names)]

for (name in layer_names) {
  test_that(paste0("roundtrip with X '", name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$layers[[name]],
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
      ad_new$X,
      data$layers[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in layer_names) {
  test_that(paste0("reticulate->hdf5 with X '", name, "'"), {
    # add rownames
    X <- data$layers[[name]]
    obs <- data.frame(row.names = rownames(data$obs))
    var <- data.frame(row.names = rownames(data$var))

    # TODO: remove this?
    if (is.matrix(X) && any(is.na(X))) {
      na_indices <- is.na(X)
      X[na_indices] <- NaN
    }

    # create anndata
    ad <- anndata::AnnData(
      X = X,
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
      ad_new$X,
      data$layers[[name]],
      tolerance = 1e-6
    )
  })
}

r2py_names <- layer_names
# TODO: re-enable -- rsparse gets converted to csparse by anndata
r2py_names <- r2py_names[!grepl("rsparse", r2py_names)]

for (name in r2py_names) {
  test_that(paste0("hdf5->reticulate with X '", name, "'"), {
    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))


    # make anndata
    ad <- AnnData(
      X = data$layers[[name]],
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE]
    )
    write_h5ad(ad, filename)

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    layer_ <- ad_new$X
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
