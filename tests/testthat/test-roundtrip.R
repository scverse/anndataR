library(testthat)

h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")
base_file <- system.file("extdata", "example.h5ad", package = "anndataR")

check_round_trip <- function(expected, type) {
  h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")
  write_h5ad(expected, h5ad_file)
  actual <- read_h5ad(h5ad_file, to = type)
  expect_equal(actual, expected)
}

check_round_trip_example <- function(type) {
  check_round_trip(expected = read_h5ad(base_file, to = type),
                   type = type)
}

for (typ in c("HDF5AnnData", "InMemoryAnnData")) {
  # test1
  test_that(paste("round trip w/ example data for", typ), {
    suppressWarnings(
      check_round_trip_example(type = typ)
    )
  })
  # test2
  test_that(paste("round trip w/ generated data for", typ), {
    adata <- dummy_data(output = typ)
    check_round_trip(expected = adata, type = typ)
  })
}

