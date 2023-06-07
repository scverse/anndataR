skip_if_no_anndata()
skip_if_not_installed("rhdf5")

# construct dummy objects
dummy <- dummy_data(10L, 20L)

test_that("test Python -> R", {
  # create anndata in python
  obs_ <- dummy$obs
  var_ <- dummy$var
  ad <- anndata::AnnData(
    X = dummy$X,
    layers = dummy$layers,
    obs = obs_,
    var = var_,
  )
  ad$obs_names <- dummy$obs_names
  ad$var_names <- dummy$var_names

  # write to file
  filename <- withr::local_file("python_to_r.h5ad")
  ad$write_h5ad(filename)

  # read from file
  ad_new <- HDF5AnnData$new(filename)

  # Python writer coerces strings to categoricals (in most cases)
  obs_$cell_type <- factor(obs_$cell_type)
  var_$geneinfo <- factor(var_$geneinfo)

  # expect slots are unchanged
  expect_equal(ad_new$X, dummy$X, tolerance = 1e-7)
  expect_equal(ad_new$obs, obs_, tolerance = 1e-10)
  expect_equal(ad_new$var, var_, tolerance = 1e-10)
  expect_equal(ad_new$obs_names, dummy$obs_names, tolerance = 1e-10)
  expect_equal(ad_new$var_names, dummy$var_names, tolerance = 1e-10)
  expect_equal(ad_new$layers, dummy$layers, tolerance = 1e-10)
})

test_that("test R -> Python", {
  # write to file
  filename <- withr::local_file("r_to_python.h5ad")
  ad <- HDF5AnnData$new(
    file = filename,
    X = dummy$X,
    obs = dummy$obs,
    var = dummy$var,
    obs_names = dummy$obs_names,
    var_names = dummy$var_names,
    layers = dummy$layers
  )

  # read from file
  ad_new <- anndata::read_h5ad(filename)

  # expect slots are unchanged
  X2 <- ad_new$X
  dimnames(X2) <- list(NULL, NULL)
  expect_equal(X2, dummy$X, tolerance = 1e-10)

  obs_ <- ad_new$obs
  rownames(obs_) <- NULL
  expect_equal(obs_, dummy$obs, ignore_attr = TRUE, tolerance = 1e-10)

  var_ <- ad_new$var
  rownames(var_) <- NULL
  expect_equal(var_, dummy$var, ignore_attr = TRUE, tolerance = 1e-10)

  expect_equal(ad_new$obs_names, dummy$obs_names, tolerance = 1e-10) #nolint

  expect_equal(ad_new$var_names, dummy$var_names, tolerance = 1e-10) #nolint

  expect_equal(names(ad_new$layers), names(dummy$layers))
  for (layer_name in names(dummy$layers)) {
    layer_ <- ad_new$layers[[layer_name]]
    dimnames(layer_) <- list(NULL, NULL)
    expect_equal(
      layer_,
      dummy$layers[[layer_name]],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  }
})
