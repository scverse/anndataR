test_that("AnnDataView basic subsetting methods work", {
  # Create a base AnnData object using the dataset generator
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Test that basic subsetting via [ returns AnnDataView
  # Using evaluated logical conditions instead of expressions
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  view_obs <- ad[obs_condition, ]
  expect_s3_class(view_obs, "AnnDataView")

  var_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5
  view_var <- ad[, var_condition]
  expect_s3_class(view_var, "AnnDataView")
})

test_that("AnnDataView obs subsetting works", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Test logical vector subsetting
  view1 <- ad[c(TRUE, FALSE, TRUE), ]

  expect_equal(nrow(view1$obs), 2L)
  expect_equal(view1$obs_names, c("cell1", "cell3"))
  expect_equal(dim(view1$X), c(2L, 5L))
  expect_identical(view1$X, ad$X[c(1, 3), , drop = FALSE])

  # Test evaluated logical condition subsetting
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  view2 <- ad[obs_condition, ]

  expect_equal(nrow(view2$obs), 2L)
  expect_equal(view2$obs_names, c("cell1", "cell3"))
  expect_equal(as.character(view2$obs$factor), c("Value1", "Value1"))

  # Test complex evaluated condition
  complex_condition <- ad$obs$factor == "Value1" & ad$obs$integer > 0
  view3 <- ad[complex_condition, ]

  expect_equal(nrow(view3$obs), 1L)
  expect_equal(view3$obs_names, "cell3")
})

test_that("AnnDataView var subsetting works", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Test logical vector subsetting
  view1 <- ad[, c(TRUE, FALSE, TRUE, FALSE, FALSE)]

  expect_equal(nrow(view1$var), 2L)
  expect_equal(view1$var_names, c("gene1", "gene3"))
  expect_equal(dim(view1$X), c(3L, 2L))
  expect_identical(view1$X, ad$X[, c(1, 3), drop = FALSE])

  # Test evaluated logical condition subsetting
  gene_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5
  view2 <- ad[, gene_condition]

  expect_equal(nrow(view2$var), 3L)
  expect_equal(view2$var_names, c("gene1", "gene3", "gene5"))
  expect_equal(as.character(view2$var$factor), c("Value1", "Value1", "Value1"))

  # Test integer column subsetting
  view3 <- ad[, ad$var$integer >= 2] # Should select gene3, gene4, gene5

  expect_equal(nrow(view3$var), 3L)
  expect_equal(view3$var_names, c("gene3", "gene4", "gene5"))
})

test_that("AnnDataView combined subsetting works", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Create logical conditions
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  var_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5

  # Test combined obs and var subsetting
  view <- ad[obs_condition, var_condition]

  expect_equal(dim(view$X), c(2L, 3L)) # 2 obs × 3 vars
  expect_equal(view$obs_names, c("cell1", "cell3"))
  expect_equal(view$var_names, c("gene1", "gene3", "gene5"))

  # Check that obsm is subset correctly (obs on rows)
  expect_equal(nrow(view$obsm$numeric_matrix), 2L)
  expect_equal(ncol(view$obsm$numeric_matrix), 3L) # columns unchanged

  # Check that varm is subset correctly (var on rows)
  expect_equal(nrow(view$varm$numeric_matrix), 3L) # 3 selected vars
  expect_equal(ncol(view$varm$numeric_matrix), 5L) # columns unchanged

  # Check that obsp is subset correctly (obs on both dimensions)
  expect_equal(dim(view$obsp$numeric_matrix), c(2L, 2L))

  # Check that varp is subset correctly (var on both dimensions)
  expect_equal(dim(view$varp$numeric_matrix), c(3L, 3L))
})

test_that("AnnDataView layers subsetting works", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Create logical conditions
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  var_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5

  view <- ad[obs_condition, var_condition]

  # Check that all layers are subset correctly
  # Note: Layers should have the same dimensions as the view's X matrix
  expect_equal(dim(view$X), c(2L, 3L))
  expect_equal(dim(view$layers$numeric_matrix), c(2L, 3L))
  expect_equal(dim(view$layers$numeric_dense), c(2L, 3L))
  expect_equal(dim(view$layers$numeric_csparse), c(2L, 3L))

  # Verify that subsetting maintains the correct values
  expect_identical(
    view$layers$numeric_matrix,
    ad$layers$numeric_matrix[c(1, 3), c(1, 3, 5), drop = FALSE]
  )
})

test_that("AnnDataView no nested views", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Create logical conditions
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  var_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5

  # Test that subsetting a view doesn't create nested views
  view1 <- ad[obs_condition, ]
  view2 <- view1[, var_condition]

  expect_s3_class(view2, "AnnDataView")
  expect_equal(dim(view2$X), c(2L, 3L)) # 2 obs (cell1, cell3) × 3 vars (gene1, gene3, gene5)
  expect_equal(view2$obs_names, c("cell1", "cell3"))
  expect_equal(view2$var_names, c("gene1", "gene3", "gene5"))

  # The result should be equivalent to doing both subsets at once
  view_combined <- ad[obs_condition, var_condition]
  expect_identical(view2$X, view_combined$X)
  expect_identical(view2$obs, view_combined$obs)
  expect_identical(view2$var, view_combined$var)
})

test_that("AnnDataView uns field are unaffected by subsetting", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Create logical condition
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  view <- ad[obs_condition, ]

  # uns should be unaffected by subsetting
  expect_identical(view$uns, ad$uns)
  expect_equal(view$uns$scalar_character, "value_0")
  expect_equal(length(view$uns$vec_integer), 10L)
})

test_that("AnnDataView conversion to InMemoryAnnData works", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Create logical conditions
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  var_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5

  # Create a view with subsetting
  view <- ad[obs_condition, var_condition]

  # Convert to InMemoryAnnData
  result <- view$as_InMemoryAnnData()

  expect_s3_class(result, "InMemoryAnnData")
  expect_equal(dim(result$X), c(2L, 3L))
  expect_equal(result$obs_names, c("cell1", "cell3"))
  expect_equal(result$var_names, c("gene1", "gene3", "gene5"))
  expect_equal(as.character(result$obs$factor), c("Value1", "Value1"))
  expect_equal(as.character(result$var$factor), c("Value1", "Value1", "Value1"))
  
  # Check that layers are preserved
  expect_equal(length(result$layers), length(ad$layers))
  expect_equal(dim(result$layers$numeric_matrix), c(2L, 3L))
})

test_that("AnnDataView error handling works", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Test initialization with non-AbstractAnnData
  expect_error(
    AnnDataView$new("not_an_anndata"),
    "must be an AbstractAnnData object"
  )

  # Test setting values on AnnDataView (create view first)
  # Create logical condition
  obs_condition <- ad$obs$factor == "Value1"

  view <- ad[obs_condition, ]
  expect_error(
    view$X <- matrix(1:15, 3L, 5L),
    "Cannot set X on a AnnDataView object"
  )
  expect_error(
    view$obs <- data.frame(),
    "Cannot set obs on a AnnDataView object"
  )
  expect_error(
    view$var <- data.frame(),
    "Cannot set var on a AnnDataView object"
  )
  expect_error(
    view$layers <- list(),
    "Cannot set layers on a AnnDataView object"
  )

  # Test that AnnDataView inherits from AbstractAnnData
  # Create logical condition
  obs_condition <- ad$obs$factor == "Value1"

  view2 <- ad[obs_condition, ]
  expect_s3_class(view2, "AbstractAnnData")
  expect_s3_class(view2, "AnnDataView")
})

test_that("AnnDataView subsetting fails with invalid conditions", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Test invalid logical vector length
  expect_error(
    ad[c(TRUE, FALSE), ], # Wrong length
    "Logical subset of observations must have length 3"
  )

  expect_error(
    ad[, c(TRUE, FALSE, TRUE)], # Wrong length
    "Logical subset of variables must have length 5"
  )

  # Test invalid condition types
  expect_error(
    ad["invalid", ], # String that's not a valid name
    "Names of observations not found: invalid"
  )

  # Test with character subsetting when only auto-generated names exist
  ad_auto_names <- AnnData(X = matrix(1:15, 3L, 5L))
  # Auto-generated names are "1", "2", "3", so use names that don't exist
  expect_error(
    ad_auto_names[c("invalid1", "invalid2"), ],
    "Names of observations not found: invalid1, invalid2"
  )
})
