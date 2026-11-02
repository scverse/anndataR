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
    normalise_py_anndata_print("    layers: None (.X), 'counts'"),
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

test_that("normalise_py_anndata_print() leaves quoted layer names alone", {
  expect_equal(
    normalise_py_anndata_print("    layers: 'counts, None (.X)'"),
    "    layers: 'counts, None (.X)'"
  )
  expect_equal(
    normalise_py_anndata_print("    layers: 'None (.X)'"),
    "    layers: 'None (.X)'"
  )
  # the real entry is still removed when it follows such a name
  expect_equal(
    normalise_py_anndata_print("    layers: 'a, None (.X), b', None (.X)"),
    "    layers: 'a, None (.X), b'"
  )
  expect_equal(
    normalise_py_anndata_print("    layers: None (.X), 'counts, None (.X)'"),
    "    layers: 'counts, None (.X)'"
  )
})
