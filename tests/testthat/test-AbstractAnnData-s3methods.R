skip_if_not_installed("Matrix")

test_that("S3 method dim.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  result <- dim(ad)
  expect_type(result, "integer")
  expect_length(result, 2)
  expect_equal(result, c(10, 5))
  expect_equal(result[1], ad$n_obs())
  expect_equal(result[2], ad$n_vars())
})

test_that("S3 method nrow.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  result <- nrow(ad)
  expect_type(result, "integer")
  expect_length(result, 1)
  expect_equal(result, 10)
  expect_equal(result, ad$n_obs())
})

test_that("S3 method ncol.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  result <- ncol(ad)
  expect_type(result, "integer")
  expect_length(result, 1)
  expect_equal(result, 5)
  expect_equal(result, ad$n_vars())
})

test_that("S3 method rownames.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

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

  result <- colnames(ad)
  expect_type(result, "character")
  expect_length(result, 5)
  expect_equal(result, ad$var_names)

  # Test with custom var names
  custom_names <- paste0("gene_", 1:5)
  ad$var_names <- custom_names
  expect_equal(colnames(ad), custom_names)
})

test_that("S3 method rownames<-.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  new_names <- paste0("cell_", 1:10)
  rownames(ad) <- new_names

  expect_equal(rownames(ad), new_names)
  expect_equal(ad$obs_names, new_names)

  # Test that the function modifies the object correctly
  new_names2 <- paste0("obs_", 1:10)
  rownames(ad) <- new_names2
  expect_equal(rownames(ad), new_names2)

  # Test assignment result (R's standard behavior for assignment operators)
  result <- rownames(ad) <- paste0("final_obs_", 1:10)
  expect_equal(result, paste0("final_obs_", 1:10))
  expect_equal(rownames(ad), paste0("final_obs_", 1:10))

  # Test with NULL (should work if the underlying implementation supports it)
  expect_no_error({
    rownames(ad) <- NULL
  })

  # Check that rownames are reset to default values
  default_names <- as.character(seq_len(10))
  expect_equal(rownames(ad), default_names)
  expect_equal(ad$obs_names, default_names)
})

test_that("S3 method colnames<-.AbstractAnnData works", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  new_names <- paste0("gene_", 1:5)
  colnames(ad) <- new_names

  expect_equal(colnames(ad), new_names)
  expect_equal(ad$var_names, new_names)

  # Test that the function modifies the object correctly
  new_names2 <- paste0("var_", 1:5)
  colnames(ad) <- new_names2
  expect_equal(colnames(ad), new_names2)

  # Test assignment result (R's standard behavior for assignment operators)
  result <- colnames(ad) <- paste0("final_var_", 1:5)
  expect_equal(result, paste0("final_var_", 1:5))
  expect_equal(colnames(ad), paste0("final_var_", 1:5))

  # Test with NULL (should work if the underlying implementation supports it)
  expect_no_error({
    colnames(ad) <- NULL
  })

  # Check that colnames are reset to default values
  default_names <- as.character(seq_len(5))
  expect_equal(colnames(ad), default_names)
  expect_equal(ad$var_names, default_names)
})

test_that("S3 method [.AbstractAnnData works with numeric indices", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  result <- ad[1:3, ]
  expect_s3_class(result, "AnnDataView")
  expect_s3_class(result, "AbstractAnnData")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
  expect_equal(rownames(result), rownames(ad)[1:3])
  expect_equal(colnames(result), colnames(ad))

  result <- ad[, 2:4]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
  expect_equal(rownames(result), rownames(ad))
  expect_equal(colnames(result), colnames(ad)[2:4])

  result <- ad[1:3, 2:4]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_equal(rownames(result), rownames(ad)[1:3])
  expect_equal(colnames(result), colnames(ad)[2:4])
})

test_that("S3 method [.AbstractAnnData works with logical indices", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

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

  logical_cols <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  result <- ad[, logical_cols]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), colnames(ad)[logical_cols])

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

  selected_obs <- c("cell_1", "cell_3", "cell_5")
  result <- ad[selected_obs, ]
  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
  expect_equal(rownames(result), selected_obs)

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
    "Logical subset of observations must have length 10"
  )

  # Test invalid logical vector length for variables
  expect_error(
    ad[, c(TRUE, FALSE)],
    "Logical subset of variables must have length 5"
  )

  # Test out-of-bounds numeric indices for observations
  expect_error(
    ad[11, ],
    "Integer indices for observations must be between 1 and 10"
  )

  expect_error(
    ad[0, ],
    "Integer indices for observations must be between 1 and 10"
  )

  # Test out-of-bounds numeric indices for variables
  expect_error(
    ad[, 6],
    "Integer indices for variables must be between 1 and 5"
  )

  expect_error(
    ad[, 0],
    "Integer indices for variables must be between 1 and 5"
  )

  # Test invalid character names for observations
  expect_error(
    ad[c("nonexistent_obs"), ],
    "Names of observations not found: nonexistent_obs"
  )

  # Test invalid character names for variables
  expect_error(
    ad[, c("nonexistent_var")],
    "Names of variables not found: nonexistent_var"
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
  expect_output(print(h5_ad), "InMemoryAnnData object")
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
  expect_output(print(ad), "InMemoryAnnData object")
})

test_that("S3 methods integration: setting names and subsetting", {
  ad <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")

  # Set custom names using S3 setter methods
  custom_obs_names <- paste0("cell_", 1:10)
  custom_var_names <- paste0("gene_", 1:5)

  rownames(ad) <- custom_obs_names
  colnames(ad) <- custom_var_names

  # Verify names were set correctly
  expect_equal(rownames(ad), custom_obs_names)
  expect_equal(colnames(ad), custom_var_names)
  expect_equal(ad$obs_names, custom_obs_names)
  expect_equal(ad$var_names, custom_var_names)

  # Test subsetting using the custom names
  selected_cells <- c("cell_1", "cell_3", "cell_5")
  selected_genes <- c("gene_1", "gene_2")

  result <- ad[selected_cells, selected_genes]

  expect_s3_class(result, "AnnDataView")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  expect_equal(rownames(result), selected_cells)
  expect_equal(colnames(result), selected_genes)

  # Test that the subsetting actually worked by checking dimensions
  expect_identical(result$X, ad$X[c(1, 3, 5), c(1, 2), drop = FALSE])

  # Test mixed subsetting: numeric for rows, names for columns
  mixed_result <- ad[1:3, selected_genes]
  expect_equal(nrow(mixed_result), 3)
  expect_equal(ncol(mixed_result), 2)
  expect_equal(rownames(mixed_result), custom_obs_names[1:3])
  expect_equal(colnames(mixed_result), selected_genes)
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
