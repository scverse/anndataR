test_that("zarr_chunk_dim matches the chunks of Python's zarr", {
  # zarr.core.chunk_grids._guess_chunks(shape, typesize) in zarr 3.1
  expect_equal(zarr_chunk_dim(3e6, "double"), 93750L)
  expect_equal(zarr_chunk_dim(c(3e6, 4), "double"), c(187500L, 1L))
  expect_equal(zarr_chunk_dim(c(1e5, 2e4), "double"), c(1563L, 313L))
  expect_equal(zarr_chunk_dim(3e6, "integer"), 187500L)
  expect_equal(zarr_chunk_dim(1e8, "double"), 195313L)
})

test_that("zarr_chunk_dim keeps small arrays as one chunk", {
  expect_equal(zarr_chunk_dim(c(10L, 10L), "double"), c(10L, 10L))
  expect_equal(zarr_chunk_dim(5L, "double"), 5L)
})

test_that("zarr_chunk_dim with chunk_size = NULL gives a single chunk", {
  expect_equal(
    zarr_chunk_dim(c(3e6, 4), "double", chunk_size = NULL),
    c(3e6, 4)
  )
})

test_that("zarr_chunk_dim uses a number as the target chunk size in bytes", {
  result <- zarr_chunk_dim(c(1e4, 10), "double", chunk_size = 64L * 1024L)
  expect_lte(prod(result) * 8, 64 * 1024 * 1.5)
  expect_true(all(result <= c(1e4, 10)))
})

test_that("zarr_chunk_dim leaves scalars and empty arrays alone", {
  expect_equal(zarr_chunk_dim(integer(0), "double"), integer(0))
  expect_equal(zarr_chunk_dim(c(0L, 5L), "double"), c(0L, 5L))
})

test_that("zarr_chunk_dim sizes elements as Python anndata writes them", {
  # expected values from zarr.core.chunk_grids._guess_chunks() with the item
  # size Python anndata writes: bool 1, int32 4, float64 8, strings 1
  expect_equal(zarr_chunk_dim(1e8, "logical"), 781250L)
  expect_equal(zarr_chunk_dim(3e6, "character"), 375000L)
  expect_equal(zarr_chunk_dim(3e5, "character"), 150000L)
  expect_equal(zarr_chunk_dim(3e6, "integer"), 187500L)
  expect_equal(zarr_chunk_dim(c(3e5, 50), "integer"), c(37500L, 7L))
  expect_equal(zarr_chunk_dim(c(1e6, 2e3), "double"), c(15625L, 32L))
})
