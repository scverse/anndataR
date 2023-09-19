h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")

base_file <- system.file("extdata", "example.h5ad", package = "anndataR")


check_round_trip <- function(type) {
  h5ad_file <- tempfile(pattern = "hdf5_write_", fileext = ".h5ad")
  orig <- read_h5ad(base_file, to = type)
  write_h5ad(orig, h5ad_file)
  new <- read_h5ad(h5ad_file, to = type)

  expect_equal(new, orig)
}

for (typ in c("HDF5AnnData", "InMemoryAnnData")) {
  test_that(paste("round trip for", typ), {
    check_round_trip(typ)
  })
}
