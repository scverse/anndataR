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


args1 <- formals(generate_dataset)
args2 <- formals(.generate_dataset_as_list)

# If any of these tests are failing, the arguments of generate_dataset()
# need to be updated to include the new generator types.
for (arg in intersect(names(args1), names(args2))) {
  test_that(paste0("generate_dataset(): argument '", arg, "' has correct default value"), {
    expect_equal(eval(args1[[arg]]), eval(args2[[arg]]))
  })
}
