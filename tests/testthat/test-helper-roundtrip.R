for (x_key in c("NA (.X)", "None (.X)")) {
  test_that(
    paste("normalise_anndata_print() drops", x_key, "from the layers line"),
    {
      # only the x layer
      expect_equal(
        normalise_anndata_print(c(
          "AnnData object with n_obs × n_vars = 10 × 20",
          paste0("    layers: ", x_key)
        )),
        "AnnData object with n_obs × n_vars = 10 × 20"
      )
      # trailing
      expect_equal(
        normalise_anndata_print(paste0(
          "    layers: 'counts', 'data', ",
          x_key
        )),
        "    layers: 'counts', 'data'"
      )
      # leading
      expect_equal(
        normalise_anndata_print(paste0("    layers: ", x_key, ", 'counts'")),
        "    layers: 'counts'"
      )
      # middle
      expect_equal(
        normalise_anndata_print(paste0(
          "    layers: 'counts', ",
          x_key,
          ", 'data'"
        )),
        "    layers: 'counts', 'data'"
      )
    }
  )
}

test_that("normalise_anndata_print() leaves other lines alone", {
  expect_equal(
    normalise_anndata_print("    layers: 'counts'"),
    "    layers: 'counts'"
  )
  expect_equal(
    normalise_anndata_print("    obsm: 'X_pca'"),
    "    obsm: 'X_pca'"
  )
  expect_equal(
    normalise_anndata_print("    uns: 'NA (.X)'"),
    "    uns: 'NA (.X)'"
  )
})
