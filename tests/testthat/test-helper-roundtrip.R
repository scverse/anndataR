test_that("normalise_py_anndata_print() drops the X entry from the layers line", {
  expect_equal(
    normalise_py_anndata_print(c(
      "AnnData object with n_obs × n_vars = 10 × 20",
      "    layers: None (.X)"
    )),
    "AnnData object with n_obs × n_vars = 10 × 20"
  )
  expect_equal(
    normalise_py_anndata_print("    layers: 'counts', None (.X)"),
    "    layers: 'counts'"
  )
  expect_equal(
    normalise_py_anndata_print("    layers: 'counts', 'data', None (.X)"),
    "    layers: 'counts', 'data'"
  )
  expect_equal(
    normalise_py_anndata_print("    layers: 'counts'"),
    "    layers: 'counts'"
  )
  expect_equal(
    normalise_py_anndata_print("    obsm: 'X_pca'"),
    "    obsm: 'X_pca'"
  )
})
