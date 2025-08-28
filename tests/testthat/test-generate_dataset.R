test_that("generating dummy data works", {
  dataset <- generate_dataset()
  expect_type(dataset, "list")
  expect_setequal(
    names(dataset),
    .anndata_slots
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


args1 <- formals(generate_dataset)
args2 <- formals(.generate_dataset_as_list)

# If any of these tests are failing, the arguments of generate_dataset()
# need to be updated to include the new generator types.
for (arg in intersect(names(args1), names(args2))) {
  test_that(
    paste0(
      "generate_dataset(): argument '",
      arg,
      "' has correct default value"
    ),
    {
      expect_equal(eval(args1[[arg]]), eval(args2[[arg]]))
    }
  )
}

test_that("Generated obsm/varm have expected dimensions", {
  n_obs <- 10L
  n_vars <- 20L
  dummy <- generate_dataset(
    n_obs = n_obs,
    n_vars = n_vars,
    x_type = NULL,
    obs_types = list(),
    var_types = list(),
    layer_types = list(),
    obsm_types = list("integer_matrix"),
    varm_types = list("integer_matrix"),
    obsp_types = list(),
    varp_types = list(),
    uns_types = list(),
  )

  expect_identical(dim(dummy$obsm$integer_matrix), c(n_obs, n_obs))
  expect_identical(dim(dummy$varm$integer_matrix), c(n_vars, n_vars))
})

test_that("generate_dataset() works with empty types", {
  expect_no_error(
    generate_dataset(
      x_type = NULL,
      obs_types = list(),
      var_types = list(),
      layer_types = list(),
      obsm_types = list(),
      varm_types = list(),
      obsp_types = list(),
      varp_types = list(),
      uns_types = list()
    )
  )
})

test_that("generated logical vectors are as expected", {
  dummy <- generate_dataset(
    n_obs = 10L,
    x_type = NULL,
    obs_types = list("logical", "logical_with_nas"),
    var_types = list(),
    layer_types = list(),
    obsm_types = list(),
    varm_types = list(),
    obsp_types = list(),
    varp_types = list(),
    uns_types = list()
  )

  expect_true(dummy$obs$logical[1])
  expect_false(dummy$obs$logical[2])
  expect_true(sum(dummy$obs$logical == TRUE) == 5)
  expect_true(sum(dummy$obs$logical == FALSE) == 5)

  expect_true(is.na(dummy$obs$logical_with_nas[1]))
  expect_false(dummy$obs$logical_with_nas[2])
  expect_true(dummy$obs$logical_with_nas[3])
  expect_true(sum(dummy$obs$logical_with_nas == TRUE, na.rm = TRUE) == 4)
  expect_true(sum(dummy$obs$logical_with_nas == FALSE, na.rm = TRUE) == 5)
  expect_true(sum(is.na(dummy$obs$logical_with_nas)) == 1)
})
