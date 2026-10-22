test_that("concat works with basic functionality", {
  # Create test datasets
  ad1 <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")
  ad2 <- generate_dataset(n_obs = 8, n_vars = 5, format = "AnnData") 
  
  # Test basic row concatenation
  result <- concat(list(ad1, ad2), axis = "obs")
  
  expect_s3_class(result, "InMemoryAnnData")
  expect_s3_class(result, "AbstractAnnData")
  expect_equal(nrow(result), 18)  # 10 + 8
  expect_equal(ncol(result), 5)
  
  # Test rbind
  result_rbind <- rbind(ad1, ad2)
  expect_equal(dim(result_rbind), c(18, 5))
  
  # Test cbind with different number of variables
  ad3 <- generate_dataset(n_obs = 10, n_vars = 3, format = "AnnData")
  result_cbind <- cbind(ad1, ad3)
  expect_equal(dim(result_cbind), c(10, 8))  # 5 + 3
})

test_that("concat handles mismatched dimensions correctly", {
  ad1 <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")
  ad2 <- generate_dataset(n_obs = 10, n_vars = 3, format = "AnnData")
  
  # Make variable names different to test joins properly
  ad1$var_names <- paste0("gene_A_", 1:5)
  ad2$var_names <- paste0("gene_B_", 1:3)
  
  # Inner join should take intersection of variables (none in this case)
  result_inner <- concat(list(ad1, ad2), axis = "obs", join = "inner")
  expect_equal(ncol(result_inner), 0)  # No common variables
  
  # Outer join should take union
  result_outer <- concat(list(ad1, ad2), axis = "obs", join = "outer") 
  expect_equal(ncol(result_outer), 8)  # 5 + 3 unique variables
  expect_equal(nrow(result_outer), 20)  # 10 + 10
})

test_that("concat handles empty inputs gracefully", {
  expect_error(concat(list()), "adatas must be a non-empty list")
  
  ad1 <- generate_dataset(n_obs = 10, n_vars = 5, format = "AnnData")
  
  # Single object should return copy with warning
  expect_warning(result <- concat(list(ad1)), "Only one object provided")
  expect_equal(dim(result), dim(ad1))
})

test_that("concat preserves metadata correctly", {
  ad1 <- generate_dataset(n_obs = 5, n_vars = 3, format = "AnnData")
  ad2 <- generate_dataset(n_obs = 4, n_vars = 3, format = "AnnData")
  
  # Add some test metadata
  ad1$uns$test_meta <- "from_ad1"
  ad2$uns$test_meta <- "from_ad2"
  ad1$uns$unique_to_ad1 <- "only_in_ad1"
  
  result <- concat(list(ad1, ad2), axis = "obs", label = "batch")
  
  # Check batch labels were added
  expect_true("batch" %in% colnames(result$obs))
  expect_equal(levels(result$obs$batch), c("1", "2"))
  
  # Check metadata merging (should use "unique" strategy by default)
  expect_true("unique_to_ad1" %in% names(result$uns))
})
