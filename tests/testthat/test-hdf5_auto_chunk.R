test_that("hdf5_auto_chunk returns integer vector", {
  result <- hdf5_auto_chunk(c(100L, 50L), "double")
  expect_type(result, "integer")
})

test_that("hdf5_auto_chunk preserves dimensionality", {
  result_1d <- hdf5_auto_chunk(1000L, "double")
  expect_length(result_1d, 1)

  result_2d <- hdf5_auto_chunk(c(100L, 50L), "double")
  expect_length(result_2d, 2)

  result_3d <- hdf5_auto_chunk(c(100L, 50L, 30L), "double")
  expect_length(result_3d, 3)
})

test_that("hdf5_auto_chunk doesn't exceed original dims", {
  dims <- c(2000L, 1000L)
  result <- hdf5_auto_chunk(dims, "double")
  expect_true(all(result <= dims))
})

test_that("hdf5_auto_chunk keeps small datasets as-is", {
  # 5 doubles = 40 bytes, well under any chunk target
  result <- hdf5_auto_chunk(5L, "double")
  expect_equal(result, 5L)

  # 10x10 doubles = 800 bytes
  result <- hdf5_auto_chunk(c(10L, 10L), "double")
  expect_equal(result, c(10L, 10L))
})

test_that("hdf5_auto_chunk reduces large datasets", {
  # 2000 x 1000 doubles = 16 MB, should be chunked below 1 MB
  dims <- c(2000L, 1000L)
  result <- hdf5_auto_chunk(dims, "double")
  chunk_bytes <- prod(result) * 8L
  expect_true(chunk_bytes < 1024L * 1024L)
  expect_true(any(result < dims))
})

test_that("hdf5_auto_chunk handles very large datasets", {
  # 50000 x 20000 doubles = 8 GB
  dims <- c(50000L, 20000L)
  result <- hdf5_auto_chunk(dims, "double")
  chunk_bytes <- prod(result) * 8L
  expect_true(chunk_bytes < 1024L * 1024L)
  expect_true(all(result >= 1L))
})

test_that("hdf5_auto_chunk handles 1D vectors", {
  # 2 million doubles = 16 MB
  result <- hdf5_auto_chunk(2000000L, "double")
  chunk_bytes <- result * 8L
  expect_true(chunk_bytes < 1024L * 1024L)
  expect_true(result >= 1L)
})

test_that("hdf5_auto_chunk respects storage mode", {
  dims <- c(2000L, 1000L)

  result_double <- hdf5_auto_chunk(dims, "double")
  result_integer <- hdf5_auto_chunk(dims, "integer")

  # Integer elements are 4 bytes vs 8 bytes for double,

  # so integer chunks can be larger in element count
  expect_true(prod(result_integer) >= prod(result_double))
})

test_that("hdf5_auto_chunk handles all storage modes", {
  dims <- c(500L, 500L)

  expect_no_error(hdf5_auto_chunk(dims, "double"))
  expect_no_error(hdf5_auto_chunk(dims, "integer"))
  expect_no_error(hdf5_auto_chunk(dims, "character"))
  # Unknown mode falls back to 8 bytes
  expect_no_error(hdf5_auto_chunk(dims, "raw"))
})

test_that("hdf5_auto_chunk produces positive chunk dimensions", {
  # Edge case: very large dataset with small element size
  result <- hdf5_auto_chunk(c(100000L, 100000L), "integer")
  expect_true(all(result >= 1L))
})

test_that("hdf5_auto_chunk handles dimension of 1", {
  # Tall skinny matrix: 1000000 x 1
  result <- hdf5_auto_chunk(c(1000000L, 1L), "double")
  expect_equal(result[2], 1L)
  expect_true(result[1] <= 1000000L)

  # Wide flat matrix: 1 x 1000000
  result <- hdf5_auto_chunk(c(1L, 1000000L), "double")
  expect_equal(result[1], 1L)
  expect_true(result[2] <= 1000000L)
})

test_that("hdf5_auto_chunk respects explicit target_size", {
  dims <- c(2000L, 1000L)

  # 64 KB target
  result_64k <- hdf5_auto_chunk(dims, "double", target_size = 64L * 1024L)
  chunk_bytes_64k <- prod(result_64k) * 8L
  expect_true(chunk_bytes_64k <= 64L * 1024L * 1.5)

  # 512 KB target should produce larger chunks than 64 KB target
  result_512k <- hdf5_auto_chunk(dims, "double", target_size = 512L * 1024L)
  expect_true(prod(result_512k) >= prod(result_64k))
})

test_that("hdf5_auto_chunk avoids the 4 GB chunk limit (regression test for #415)", {
  # The original report used n_obs x n_vars = 321427 x 17690 doubles.
  # Without auto-chunking rhdf5 set chunk = dims, which at 8 bytes/element
  # gives ~45 GB per chunk - exceeding HDF5's 4 GB hard limit.
  # See https://github.com/scverse/anndataR/issues/415
  dims <- c(321427L, 17690L)
  result <- hdf5_auto_chunk(dims, "double")
  chunk_bytes <- prod(as.numeric(result)) * 8
  limit_4gb <- 4 * 1024^3
  expect_true(chunk_bytes < limit_4gb)
  # In practice it should be well within the 1 MB auto-chunk ceiling
  expect_true(chunk_bytes <= 1024L * 1024L)
})
