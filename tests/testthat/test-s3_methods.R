skip_if_not_installed("Matrix")

test_that("S3 method dim.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test dim() method
  result <- dim(ad)
  expect_type(result, "integer")
  expect_length(result, 2)
  expect_equal(result, c(10, 5))
  expect_equal(result[1], ad$n_obs())
  expect_equal(result[2], ad$n_vars())
})

test_that("S3 method nrow.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test nrow() method
  result <- nrow(ad)
  expect_type(result, "integer")
  expect_length(result, 1)
  expect_equal(result, 10)
  expect_equal(result, ad$n_obs())
})

test_that("S3 method ncol.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test ncol() method
  result <- ncol(ad)
  expect_type(result, "integer")
  expect_length(result, 1)
  expect_equal(result, 5)
  expect_equal(result, ad$n_vars())
})

test_that("S3 method rownames.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test rownames() method
  result <- rownames(ad)
  expect_type(result, "character")
  expect_length(result, 10)
  expect_equal(result, ad$obs_names)

  # Test with custom obs names
  custom_names <- paste0("cell_", 1:10)
  ad$obs_names <- custom_names
  expect_equal(rownames(ad), custom_names)
})

test_that("S3 method colnames.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test colnames() method
  result <- colnames(ad)
  expect_type(result, "character")
  expect_length(result, 5)
  expect_equal(result, ad$var_names)

  # Test with custom var names
  custom_names <- paste0("gene_", 1:5)
  ad$var_names <- custom_names
  expect_equal(colnames(ad), custom_names)
})

test_that("S3 method [.AbstractAnnData works with numeric indices", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test row subsetting
  result <- ad[1:3, ]
  expect_s3_class(result, "AnnDataView")
  expect_s3_class(result, "AbstractAnnData")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
  expect_equal(rownames(result), rownames(ad)[1:3])
  expect_equal(colnames(result), colnames(ad))

  # Test column subsetting
  result <- ad[, 2:4]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
  expect_equal(rownames(result), rownames(ad))
  expect_equal(colnames(result), colnames(ad)[2:4])

  # Test combined subsetting
  result <- ad[1:3, 2:4]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_equal(rownames(result), rownames(ad)[1:3])
  expect_equal(colnames(result), colnames(ad)[2:4])
})

test_that("S3 method [.AbstractAnnData works with logical indices", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test row subsetting with logical vector
  logical_rows <- c(
    TRUE,
    FALSE,
    TRUE,
    FALSE,
    TRUE,
    FALSE,
    TRUE,
    FALSE,
    TRUE,
    FALSE
  )
  result <- ad[logical_rows, ]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
  expect_equal(rownames(result), rownames(ad)[logical_rows])

  # Test column subsetting with logical vector
  logical_cols <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  result <- ad[, logical_cols]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), colnames(ad)[logical_cols])

  # Test combined logical subsetting
  result <- ad[logical_rows, logical_cols]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 3)
})

test_that("S3 method [.AbstractAnnData works with character indices", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Set custom names for testing
  custom_obs_names <- paste0("cell_", 1:10)
  custom_var_names <- paste0("gene_", 1:5)
  ad$obs_names <- custom_obs_names
  ad$var_names <- custom_var_names

  # Test row subsetting with character vector
  selected_obs <- c("cell_1", "cell_3", "cell_5")
  result <- ad[selected_obs, ]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
  expect_equal(rownames(result), selected_obs)

  # Test column subsetting with character vector
  selected_vars <- c("gene_2", "gene_4")
  result <- ad[, selected_vars]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), selected_vars)

  # Test combined character subsetting
  result <- ad[selected_obs, selected_vars]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  expect_equal(rownames(result), selected_obs)
  expect_equal(colnames(result), selected_vars)
})

test_that("S3 method [.AbstractAnnData error handling works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test invalid logical vector length for observations
  expect_error(
    ad[c(TRUE, FALSE), ],
    "Logical vector for observations must have length 10"
  )

  # Test invalid logical vector length for variables
  expect_error(
    ad[, c(TRUE, FALSE)],
    "Logical vector for variables must have length 5"
  )

  # Test out-of-bounds numeric indices for observations
  expect_error(
    ad[11, ],
    "Observation indices must be between 1 and 10"
  )

  expect_error(
    ad[0, ],
    "Observation indices must be between 1 and 10"
  )

  # Test out-of-bounds numeric indices for variables
  expect_error(
    ad[, 6],
    "Variable indices must be between 1 and 5"
  )

  expect_error(
    ad[, 0],
    "Variable indices must be between 1 and 5"
  )

  # Test invalid character names for observations
  expect_error(
    ad[c("nonexistent_obs"), ],
    "Observation names not found: nonexistent_obs"
  )

  # Test invalid character names for variables
  expect_error(
    ad[, c("nonexistent_var")],
    "Variable names not found: nonexistent_var"
  )
})

test_that("S3 methods work with HDF5AnnData", {
  skip_if_not_installed("rhdf5")

  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Write to HDF5 and read back
  write_h5ad(ad, h5ad_file)
  h5_ad <- read_h5ad(h5ad_file)

  # Test all S3 methods work with HDF5AnnData
  expect_equal(dim(h5_ad), c(10, 5))
  expect_equal(nrow(h5_ad), 10)
  expect_equal(ncol(h5_ad), 5)
  expect_length(rownames(h5_ad), 10)
  expect_length(colnames(h5_ad), 5)

  # Test subsetting works
  result <- h5_ad[1:3, 1:2]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)

  # Test print works
  expect_output(print(h5_ad), "AnnData object")
})

test_that("S3 methods work with InMemoryAnnData", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Test all S3 methods work with InMemoryAnnData
  expect_equal(dim(ad), c(10, 5))
  expect_equal(nrow(ad), 10)
  expect_equal(ncol(ad), 5)
  expect_length(rownames(ad), 10)
  expect_length(colnames(ad), 5)

  # Test subsetting works
  result <- ad[1:3, 1:2]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)

  # Test print works
  expect_output(print(ad), "AnnData object")
})

test_that("S3 methods are consistent across different backends", {
  skip_if_not_installed("rhdf5")

  # Create InMemoryAnnData
  mem_ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Create HDF5AnnData
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(mem_ad, h5ad_file)
  h5_ad <- read_h5ad(h5ad_file)

  # Test that S3 methods give consistent results
  expect_equal(dim(mem_ad), dim(h5_ad))
  expect_equal(nrow(mem_ad), nrow(h5_ad))
  expect_equal(ncol(mem_ad), ncol(h5_ad))
  expect_equal(rownames(mem_ad), rownames(h5_ad))
  expect_equal(colnames(mem_ad), colnames(h5_ad))

  # Test that subsetting gives consistent results
  mem_subset <- mem_ad[1:5, 1:3]
  h5_subset <- h5_ad[1:5, 1:3]

  expect_equal(dim(mem_subset), dim(h5_subset))
  expect_equal(rownames(mem_subset), rownames(h5_subset))
  expect_equal(colnames(mem_subset), colnames(h5_subset))
})
