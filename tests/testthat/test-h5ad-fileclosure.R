dummy <- generate_dataset(
  n_obs = 10L,
  n_vars = 20L,
  format = "AnnData",
)

file <- tempfile(pattern = "h5ad_write_", fileext = ".h5ad")

test_that("writing H5AD works", {
  expect_no_condition({
    write_h5ad(dummy, file, mode = "w")
  })
})

test_that("reading H5AD to InMemoryAnnData closes file", {
  expect_no_condition({
    read_h5ad(file)
    gc()
    write_h5ad(dummy, file, mode = "w")
  })
})

test_that("closing HDF5AnnData file works", {
  expect_no_condition({
    adata <- read_h5ad(file, as = "HDF5AnnData")
    adata$close()
    write_h5ad(dummy, file, mode = "w")
  })
})
