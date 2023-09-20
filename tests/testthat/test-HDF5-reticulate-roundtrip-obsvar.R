skip_if_no_anndata()
skip_if_not_installed("rhdf5")

data <- generate_dataset_as_list(10L, 20L)

for (name in names(data$obs)) {
  test_that(paste0("Python -> R with obs and var '", name, "'"), {
    ad <- anndata::AnnData(
      X = data$X,
      obs = data$obs[, name, drop = FALSE],
      var = data$var[, name, drop = FALSE]
    )
    ad$obs_names <- data$obs_names
    ad$var_names <- data$var_names

    # write to file
    filename <- withr::local_file(paste0("python_to_r_obsvar_", name, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    # expect slots are unchanged
    expect_equal(
      ad_new$obs[[name]],
      data$obs[[name]],
      tolerance = 1e-10
    )
  })
}

for (name in names(data$obs)) {
  test_that(paste0("R -> Python with obs and var '", name, "'"), {
    # write to file
    filename <- withr::local_file(paste0("r_to_python_obsvar_", name, ".h5ad"))
    ad <- HDF5AnnData$new(
      file = filename,
      X = data$X,
      obs = data$obs[, name, drop = FALSE],
      var = data$var[, name, drop = FALSE],
      obs_names = data$obs_names,
      var_names = data$var_names
    )

    # read from file
    ad_new <- anndata::read_h5ad(filename)

    # expect slots are unchanged
    obs_ <- ad_new$obs[[name]]
    dimnames(obs_) <- list(NULL, NULL)
    expect_equal(
      obs_,
      data$obs[[name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  })
}