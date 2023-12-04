test_that("generating dummy data works", {
  dataset <- generate_dataset()
  expect_type(dataset, "list")
  expect_setequal(
    names(dataset),
    c("X", "obs", "obsp", "obsm", "obs_names", "var", "varp", "varm", "var_names", "layers", "uns")
  )
  expect_identical(dim(dataset$X), c(10L, 20L))
})

test_that("generating dummy SingleCellExperiment works", {
  dummy <- generate_dataset(format = "SingleCellExperiment")
  expect_s4_class(dummy, "SingleCellExperiment")
})

suppressPackageStartupMessages(library(SeuratObject))

test_that("generating dummy Seurat works", {
  dummy <- generate_dataset(format = "Seurat")
  expect_s4_class(dummy, "Seurat")
})

test_that("arguments of `generate_dataset()` are correct", {
  args <- formals(generate_dataset)

  expect_setequal(eval(args$x_type), names(matrix_generators)[[1]])
  expect_setequal(eval(args$layer_types), names(matrix_generators))
  expect_setequal(eval(args$obs_types), names(vector_generators))
  expect_setequal(eval(args$var_types), names(vector_generators))
  expect_setequal(eval(args$obsm_types), c(names(vector_generators), names(matrix_generators)))
  expect_setequal(eval(args$varm_types), c(names(vector_generators), names(matrix_generators)))
  expect_setequal(eval(args$obsp_types), names(matrix_generators))
  expect_setequal(eval(args$varp_types), names(matrix_generators))

  # If any of these tests are failing, the arguments of generate_dataset()
  # need to be updated to include the new generator types.
})
