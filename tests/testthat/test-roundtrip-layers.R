skip_if_no_anndata()
skip_if_not_installed("rhdf5")

data <- generate_dataset_as_list(10L, 20L)

for (layer_name in names(data$layers)) {
  test_that(paste0("roundtrip with layer '", layer_name, "'"), {
    # create anndata
    ad <- AnnData(
      X = data$layers[[layer_name]],
      layers = data$layers[layer_name],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # write to file
    filename <- withr::local_file(paste0("roundtrip_layer_", layer_name, ".h5ad"))
    write_h5ad(ad, filename)

    # read from file
    ad_new <- read_h5ad(filename, to = "HDF5AnnData")

    # expect slots are unchanged
    expect_equal(
      ad_new$layers[[layer_name]],
      data$layers[[layer_name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}

for (layer_name in names(data$layers)) {
  test_that(paste0("reticulate->hdf5 with layer '", layer_name, "'"), {
    ad <- anndata::AnnData(
      X = data$layers[[layer_name]],
      layers = data$layers[layer_name]
    )
    ad$obs_names <- data$obs_names
    ad$var_names <- data$var_names

    # write to file
    filename <- withr::local_file(paste0("reticulate_to_hdf5_layer_", layer_name, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$layers[[layer_name]],
      data$layers[[layer_name]],
      tolerance = 1e-10
    )
  })
}

r2py_names <- names(data$layers)
# rsparse gets converted to csparse by anndata
r2py_names <- r2py_names[!grepl("rsparse", r2py_names)]

for (layer_name in r2py_names) {
  test_that(paste0("hdf5->reticulate with layer '", layer_name, "'"), {
    # write to file
    filename <- withr::local_file(paste0("hdf5_to_reticulate_layer_", layer_name, ".h5ad"))
    ad <- HDF5AnnData$new(
      file = filename,
      X = data$layers[[layer_name]],
      layers = data$layers[layer_name],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    layer_ <- ad_new$layers[[layer_name]]
    dimnames(layer_) <- list(NULL, NULL)
    # anndata returns these layers as CsparseMatrix
    if (grepl("rsparse", layer_name)) {
      layer_ <- as(layer_, "RsparseMatrix")
    }
    expect_equal(
      layer_,
      data$layers[[layer_name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}
