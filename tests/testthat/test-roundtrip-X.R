skip_if_no_anndata()
skip_if_not_installed("hdf5r")

data <- generate_dataset(10L, 20L)

layer_names <- names(data$layers)
# TODO: re-enable these tests
layer_names <- layer_names[!grepl("_dense", layer_names)]

for (layer_name in layer_names) {
  test_that(paste0("roundtrip with layer '", layer_name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$layers[[layer_name]],
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
      data$layers[[layer_name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}

for (name in layer_names) {
  test_that(paste0("reticulate->hdf5 with layer '", name, "'"), {
    # add rownames
    X <- data$layers[[name]]
    obs <- data.frame(row.names = rownames(data$obs))
    var <- data.frame(row.names = rownames(data$var))

    if (is.matrix(X) && any(is.na(X))) {
      # is.na(value) <- NaN gets ignored
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

for (layer_name in r2py_names) {
  test_that(paste0("hdf5->reticulate with layer '", layer_name, "'"), {
    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))

    # strip rownames
    X <- data$layers[[layer_name]]
    if (!is.null(X)) {
      dimnames(X) <- list(NULL, NULL)
    }

    # make anndata
    ad <- HDF5AnnData$new(
      file = filename,
      X = X,
      obs = data$obs[, c(), drop = FALSE],
      var = data$var[, c(), drop = FALSE]
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    layer_ <- ad_new$X
    if (!is.null(layer_)) {
      dimnames(layer_) <- list(NULL, NULL)
    }
    # anndata returns these layers as CsparseMatrix
    if (grepl("rsparse", layer_name)) {
      layer_ <- as(layer_, "RsparseMatrix")
    }
    expect_equal(
      layer_,
      data$layers[[layer_name]],
      ignore_attr = TRUE,
      tolerance = 1e-6
    )
  })
}
