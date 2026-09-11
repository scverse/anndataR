known_issues <- read_known_issues()

dummy <- generate_dataset(
  n_obs = 10L,
  n_vars = 20L,
  format = "AnnData",
)

file <- tempfile(pattern = "h5ad_write_", fileext = ".h5ad")

# TEMP(rhdf5 2.57.12): the dummy dataset contains numeric elements with NAs
# which trigger a spurious rhdf5 warning, see inst/known_issues.yaml
known_issue <- message_if_known(
  backend = "HDF5AnnData",
  slot = "layers",
  dtype = "numeric_matrix_with_nas",
  process = "write",
  known_issues = known_issues
)

test_that("writing H5AD works", {
  skip_if(!is.null(known_issue), message = known_issue)
  expect_no_condition({
    write_h5ad(dummy, file, mode = "w")
  })
})

test_that("reading H5AD to InMemoryAnnData closes file", {
  skip_if(!is.null(known_issue), message = known_issue)
  expect_no_condition({
    read_h5ad(file)
    write_h5ad(dummy, file, mode = "w")
  })
})

test_that("closing HDF5AnnData file works", {
  skip_if(!is.null(known_issue), message = known_issue)
  expect_no_condition({
    adata <- read_h5ad(file, as = "HDF5AnnData")
    write_h5ad(dummy, file, mode = "w")
  })
})
