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

  # Test numeric index subsetting
  view4 <- ad[c(1, 3), ]

  expect_equal(nrow(view4$obs), 2L)
  expect_equal(view4$obs_names, c("cell1", "cell3"))
  expect_equal(dim(view4$X), c(2L, 5L))

  # Test character name subsetting
  view5 <- ad[c("cell1", "cell3"), ]

  expect_equal(nrow(view5$obs), 2L)
  expect_equal(view5$obs_names, c("cell1", "cell3"))
  expect_equal(dim(view5$X), c(2L, 5L))
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

  # Test numeric index subsetting
  view4 <- ad[, c(1, 3, 5)]

  expect_equal(nrow(view4$var), 3L)
  expect_equal(view4$var_names, c("gene1", "gene3", "gene5"))
  expect_equal(dim(view4$X), c(3L, 3L))

  # Test character name subsetting
  view5 <- ad[, c("gene1", "gene3", "gene5")]

  expect_equal(nrow(view5$var), 3L)
  expect_equal(view5$var_names, c("gene1", "gene3", "gene5"))
  expect_equal(dim(view5$X), c(3L, 3L))
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

  # Test combined numeric index subsetting
  view2 <- ad[c(1, 3), c(1, 3, 5)]

  expect_equal(dim(view2$X), c(2L, 3L))
  expect_equal(view2$obs_names, c("cell1", "cell3"))
  expect_equal(view2$var_names, c("gene1", "gene3", "gene5"))

  # Test combined character name subsetting
  view3 <- ad[c("cell1", "cell3"), c("gene1", "gene3", "gene5")]

  expect_equal(dim(view3$X), c(2L, 3L))
  expect_equal(view3$obs_names, c("cell1", "cell3"))
  expect_equal(view3$var_names, c("gene1", "gene3", "gene5"))

  # Test mixed subsetting (logical obs, numeric var)
  view4 <- ad[ad$obs$factor == "Value1", c(1, 3, 5)]

  expect_equal(dim(view4$X), c(2L, 3L))
  expect_equal(view4$obs_names, c("cell1", "cell3"))
  expect_equal(view4$var_names, c("gene1", "gene3", "gene5"))
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

test_that("AnnDataView obsm subsetting works", {
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

  # obsm should be subset on rows (observations) but columns unchanged
  expect_equal(nrow(view$obsm$numeric_matrix), 2L) # 2 selected obs
  expect_equal(ncol(view$obsm$numeric_matrix), 3L) # columns unchanged
  expect_equal(nrow(view$obsm$numeric_dense), 2L)
  expect_equal(nrow(view$obsm$numeric_csparse), 2L)

  # Verify that subsetting maintains the correct values (obs rows only)
  expect_identical(
    view$obsm$numeric_matrix,
    ad$obsm$numeric_matrix[c(1, 3), , drop = FALSE]
  )
})

test_that("AnnDataView varm subsetting works", {
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

  # varm should be subset on rows (variables) but columns unchanged
  expect_equal(nrow(view$varm$numeric_matrix), 3L) # 3 selected vars
  expect_equal(ncol(view$varm$numeric_matrix), 5L) # columns unchanged
  expect_equal(nrow(view$varm$numeric_dense), 3L)
  expect_equal(nrow(view$varm$numeric_csparse), 3L)

  # Verify that subsetting maintains the correct values (var rows only)
  expect_identical(
    view$varm$numeric_matrix,
    ad$varm$numeric_matrix[c(1, 3, 5), , drop = FALSE]
  )
})

test_that("AnnDataView obsp subsetting works", {
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

  # obsp should be subset on both rows and columns (observations)
  expect_equal(dim(view$obsp$numeric_matrix), c(2L, 2L)) # 2x2 for selected obs
  expect_equal(dim(view$obsp$numeric_dense), c(2L, 2L))
  expect_equal(dim(view$obsp$numeric_csparse), c(2L, 2L))

  # Verify that subsetting maintains the correct values (both dimensions)
  expect_identical(
    view$obsp$numeric_matrix,
    ad$obsp$numeric_matrix[c(1, 3), c(1, 3), drop = FALSE]
  )
})

test_that("AnnDataView varp subsetting works", {
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

  # varp should be subset on both rows and columns (variables)
  expect_equal(dim(view$varp$numeric_matrix), c(3L, 3L)) # 3x3 for selected vars
  expect_equal(dim(view$varp$numeric_dense), c(3L, 3L))
  expect_equal(dim(view$varp$numeric_csparse), c(3L, 3L))

  # Verify that subsetting maintains the correct values (both dimensions)
  expect_identical(
    view$varp$numeric_matrix,
    ad$varp$numeric_matrix[c(1, 3, 5), c(1, 3, 5), drop = FALSE]
  )
})

test_that("AnnDataView avoids nested wrappers when chaining subsets", {
  ad <- generate_dataset(
    n_obs = 3L,
    n_vars = 5L,
    example = TRUE,
    format = "AnnData"
  )

  # Create logical conditions
  obs_condition <- ad$obs$factor == "Value1" # Should select cell1 and cell3
  var_condition <- ad$var$factor == "Value1" # Should select gene1, gene3, gene5

  # Test that subsetting a view doesn't create nested AnnDataView(AnnDataView(...))
  # Instead, it should compute combined indices and create a new view from the original data
  view1 <- ad[obs_condition, ]
  view2 <- view1[, var_condition]

  expect_s3_class(view2, "AnnDataView")
  expect_equal(dim(view2$X), c(2L, 3L)) # 2 obs (cell1, cell3) × 3 vars (gene1, gene3, gene5)
  expect_equal(view2$obs_names, c("cell1", "cell3"))
  expect_equal(view2$var_names, c("gene1", "gene3", "gene5"))

  # Test chaining subsets on the same dimension multiple times
  # First subset obs to select cell1 and cell3, then further subset to just cell3
  view_obs1 <- ad[obs_condition, ] # cell1, cell3
  view_obs2 <- view_obs1[c("cell3"), ] # further subset to just cell3

  expect_s3_class(view_obs2, "AnnDataView")
  expect_equal(dim(view_obs2$X), c(1L, 5L)) # 1 obs × 5 vars
  expect_equal(view_obs2$obs_names, "cell3")

  # Test chaining var subsets multiple times
  # First subset vars to gene1, gene3, gene5, then further subset to gene1, gene5
  view_var1 <- ad[, var_condition] # gene1, gene3, gene5
  view_var2 <- view_var1[, c("gene1", "gene5")] # further subset

  expect_s3_class(view_var2, "AnnDataView")
  expect_equal(dim(view_var2$X), c(3L, 2L)) # 3 obs × 2 vars
  expect_equal(view_var2$var_names, c("gene1", "gene5"))

  # Test complex chaining: obs -> var -> obs -> var
  view_complex <- ad[obs_condition, ][, var_condition][c("cell3"), ][, c(
    "gene1",
    "gene5"
  )]

  expect_s3_class(view_complex, "AnnDataView")
  expect_equal(dim(view_complex$X), c(1L, 2L)) # 1 obs × 2 vars
  expect_equal(view_complex$obs_names, "cell3")
  expect_equal(view_complex$var_names, c("gene1", "gene5"))

  # Verify that chained subsetting produces the same result as combined subsetting
  # This indirectly ensures that the implementation avoids nested wrappers
  # by computing combined indices rather than wrapping views within views
  view_combined <- ad[obs_condition, var_condition]
  expect_identical(view2$X, view_combined$X)
  expect_identical(view2$obs, view_combined$obs)
  expect_identical(view2$var, view_combined$var)

  # Additional verification that the chained view behaves correctly
  expect_equal(view2$n_obs(), 2L)
  expect_equal(view2$n_vars(), 3L)
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
    "must be an.*AbstractAnnData.*object"
  )

  # Test setting values on AnnDataView (create view first)
  # Create logical condition
  obs_condition <- ad$obs$factor == "Value1"

  view <- ad[obs_condition, ]
  expect_error(
    view$X <- matrix(1:15, 3L, 5L),
    "Cannot set X on an.*AnnDataView.*object"
  )
  expect_error(
    view$obs <- data.frame(),
    "Cannot set obs on an.*AnnDataView.*object"
  )
  expect_error(
    view$var <- data.frame(),
    "Cannot set var on an.*AnnDataView.*object"
  )
  expect_error(
    view$layers <- list(),
    "Cannot set layers on an.*AnnDataView.*object"
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

  # Test invalid numeric indices
  expect_error(
    ad[c(1, 5), ], # Index 5 is out of bounds (only 3 obs)
    "Integer indices for observations must be between 1 and 3"
  )

  expect_error(
    ad[, c(1, 10)], # Index 10 is out of bounds (only 5 vars)
    "Integer indices for variables must be between 1 and 5"
  )

  # Test negative indices (not supported)
  expect_error(
    ad[-1, ],
    "Integer indices for observations must be between 1 and 3"
  )

  expect_error(
    ad[, -1],
    "Integer indices for variables must be between 1 and 5"
  )

  # Test invalid character names
  expect_error(
    ad[c("cell1", "invalid_cell"), ],
    "Names of observations not found: invalid_cell"
  )

  expect_error(
    ad[, c("gene1", "invalid_gene")],
    "Names of variables not found: invalid_gene"
  )

  # Test with character subsetting when only auto-generated names exist
  ad_auto_names <- AnnData(X = matrix(1:15, 3L, 5L))
  # Auto-generated names are "1", "2", "3", so use names that don't exist
  expect_error(
    ad_auto_names[c("invalid1", "invalid2"), ],
    "Names of observations not found: invalid1, invalid2"
  )
})
