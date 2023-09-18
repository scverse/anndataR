dummy <- dummy_data(10L, 20L)

test_that("to_Seurat with inmemoryanndata", {
  ad <- example_data(10L, 20L)
  # running to_seurat when ad0$X is null probably doesn't make any sense
  ad0 <- InMemoryAnnData$new(
    obs_names = letters[1:5],
    var_names = LETTERS[1:10]
  )

  seu <- suppressWarnings(ad$to_Seurat())
  seu0 <- ad0$to_Seurat()

  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)

  # trackstatus: class=Seurat, feature=test_get_var_names, status=done
  expect_equal(rownames(seu), dummy$var_names)
  expect_equal(rownames(seu0), LETTERS[1:10])

  # trackstatus: class=Seurat, feature=test_get_obs_names, status=done
  expect_equal(colnames(seu), dummy$obs_names)
  expect_equal(colnames(seu0), letters[1:5])

  # check whether all obs keys are found in the seu metadata
  # seurat will have extra metadata columns
  # trackstatus: class=Seurat, feature=test_get_obs, status=done
  for (obs_key in colnames(dummy$obs)) {
    expect_true(obs_key %in% colnames(seu@meta.data))
    expect_equal(seu@meta.data[[obs_key]], dummy$obs[[obs_key]])
  }
  seu0@meta.data

  # check whether all var keys are found in the seu assay metadata
  # trackstatus: class=Seurat, feature=test_get_var, status=done
  active_assay <- seu@assays[[seu@active.assay]]
  for (var_key in colnames(dummy$var)) {
    expect_true(var_key %in% colnames(active_assay@meta.features))
    expect_equal(active_assay@meta.features[[var_key]], dummy$var[[var_key]])
  }
})

test_that("to_Seurat() fails gracefully", {
  expect_error(to_Seurat(), regexp = "obj.*is missing")
  expect_error(to_Seurat("foo"), regexp = "AbstractAnnData.*not TRUE")
})

# TODO: test from_Seurat
