skip_if_no_anndata()

# construct dummy objects
dummy <- dummy_data(10L, 20L)

test_that("test Python -> R", {
  # create anndata in python
  obs_ <- dummy$obs
  rownames(obs_) <- dummy$obs_names
  var_ <- dummy$var
  rownames(var_) <- dummy$var_names
  ad <- anndata::AnnData(
    X = dummy$X,
    layers = dummy$layers,
    obs = obs_,
    var = var_,
  )

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

  obs_ <- obs
  rownames(obs_) <- NULL
  expect_equal(obs_, dummy$obs, tolerance = 1e-10)

  var_ <- var
  rownames(var_) <- NULL
  expect_equal(var_, dummy$var_, tolerance = 1e-10)

  expect_equal(ad_new$obs_names, dummy$obs_names, tolerance = 1e-10)
  expect_equal(ad_new$var_names, dummy$var_names, tolerance = 1e-10)

  expect_equal(names(ad_new$layers), names(dummy$layers))
  for (layer_name in names(dummy$layers)) {
    expect_equal(
      ad_new$layers[[layer_name]],
      dummy$layers[[layer_name]],
      ignore_attributes = TRUE,
      tolerance = 1e-10
    )
  }
})