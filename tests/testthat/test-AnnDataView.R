test_that("AnnDataView basic subsetting methods", {
  # Create a base AnnData object
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A"),
      score = c(1.5, 2.3, 1.8)
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z"),
      important = c(TRUE, FALSE, TRUE, FALSE, TRUE)
    ),
    layers = list(
      counts = matrix(10:24, 3L, 5L),
      logcounts = matrix(log(10:24), 3L, 5L)
    ),
    obsm = list(
      X_pca = matrix(rnorm(6), 3L, 2L)
    ),
    varm = list(
      PCs = matrix(rnorm(10), 5L, 2L)
    ),
    uns = list(
      method = "test",
      params = list(a = 1, b = 2)
    )
  )

  # Test that basic subsetting via [ returns AnnDataView
  # Using evaluated logical conditions instead of expressions
  obs_condition <- ad$obs$cell_type == "A"
  view_obs <- ad[obs_condition, ]
  expect_s3_class(view_obs, "AnnDataView")

  var_condition <- ad$var$gene_type == "X"
  view_var <- ad[, var_condition]
  expect_s3_class(view_var, "AnnDataView")
})

test_that("AnnDataView obs subsetting", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A"),
      score = c(1.5, 2.3, 1.8)
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    )
  )

  # Test logical vector subsetting
  view1 <- ad[c(TRUE, FALSE, TRUE), ]

  expect_equal(nrow(view1$obs), 2L)
  expect_equal(view1$obs_names, c("A", "C"))
  expect_equal(dim(view1$X), c(2L, 5L))
  expect_identical(view1$X, ad$X[c(1, 3), , drop = FALSE])

  # Test evaluated logical condition subsetting
  obs_condition <- ad$obs$cell_type == "A"
  view2 <- ad[obs_condition, ]

  expect_equal(nrow(view2$obs), 2L)
  expect_equal(view2$obs_names, c("A", "C"))
  expect_equal(view2$obs$cell_type, c("A", "A"))

  # Test complex evaluated condition
  complex_condition <- ad$obs$cell_type == "A" & ad$obs$score > 1.6
  view3 <- ad[complex_condition, ]

  expect_equal(nrow(view3$obs), 1L)
  expect_equal(view3$obs_names, "C")
})

test_that("AnnDataView var subsetting", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A")
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z"),
      important = c(TRUE, FALSE, TRUE, FALSE, TRUE)
    )
  )

  # Test logical vector subsetting
  view1 <- ad[, c(TRUE, FALSE, TRUE, FALSE, FALSE)]

  expect_equal(nrow(view1$var), 2L)
  expect_equal(view1$var_names, c("a", "c"))
  expect_equal(dim(view1$X), c(3L, 2L))
  expect_identical(view1$X, ad$X[, c(1, 3), drop = FALSE])

  # Test evaluated logical condition subsetting
  gene_condition <- ad$var$gene_type == "X"
  view2 <- ad[, gene_condition]

  expect_equal(nrow(view2$var), 2L)
  expect_equal(view2$var_names, c("a", "c"))
  expect_equal(view2$var$gene_type, c("X", "X"))

  # Test boolean column subsetting
  view3 <- ad[, ad$var$important]

  expect_equal(nrow(view3$var), 3L)
  expect_equal(view3$var_names, c("a", "c", "e"))
})

test_that("AnnDataView combined subsetting", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A")
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    ),
    obsm = list(
      X_pca = matrix(rnorm(6), 3L, 2L)
    ),
    varm = list(
      PCs = matrix(rnorm(10), 5L, 2L)
    ),
    obsp = list(
      connectivities = matrix(rnorm(9), 3L, 3L)
    ),
    varp = list(
      correlations = matrix(rnorm(25), 5L, 5L)
    )
  )

  # Create logical conditions
  obs_condition <- ad$obs$cell_type == "A"
  var_condition <- ad$var$gene_type == "X"

  # Test combined obs and var subsetting
  view <- ad[obs_condition, var_condition]

  expect_equal(dim(view$X), c(2L, 2L))
  expect_equal(view$obs_names, c("A", "C"))
  expect_equal(view$var_names, c("a", "c"))

  # Check that obsm is subset correctly (obs on rows)
  expect_equal(nrow(view$obsm$X_pca), 2L)
  expect_equal(ncol(view$obsm$X_pca), 2L) # columns unchanged

  # Check that varm is subset correctly (var on rows)
  expect_equal(nrow(view$varm$PCs), 2L)
  expect_equal(ncol(view$varm$PCs), 2L) # columns unchanged

  # Check that obsp is subset correctly (obs on both dimensions)
  expect_equal(dim(view$obsp$connectivities), c(2L, 2L))

  # Check that varp is subset correctly (var on both dimensions)
  expect_equal(dim(view$varp$correlations), c(2L, 2L))
})

test_that("AnnDataView layers subsetting", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A")
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    ),
    layers = list(
      counts = matrix(10:24, 3L, 5L),
      logcounts = matrix(log(10:24), 3L, 5L)
    )
  )

  # Create logical conditions
  obs_condition <- ad$obs$cell_type == "A"
  var_condition <- ad$var$gene_type == "X"

  view <- ad[obs_condition, var_condition]

  # Check that all layers are subset correctly
  expect_equal(dim(view$layers$counts), c(2L, 2L))
  expect_equal(dim(view$layers$logcounts), c(2L, 2L))
  expect_identical(
    view$layers$counts,
    ad$layers$counts[c(1, 3), c(1, 3), drop = FALSE]
  )
})

test_that("AnnDataView no nested views", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A")
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    )
  )

  # Create logical conditions
  obs_condition <- ad$obs$cell_type == "A"
  var_condition <- ad$var$gene_type == "X"

  # Test that subsetting a view doesn't create nested views
  view1 <- ad[obs_condition, ]
  view2 <- view1[, var_condition]

  expect_s3_class(view2, "AnnDataView")
  expect_equal(dim(view2$X), c(2L, 2L))
  expect_equal(view2$obs_names, c("A", "C"))
  expect_equal(view2$var_names, c("a", "c"))

  # The result should be equivalent to doing both subsets at once
  view_combined <- ad[obs_condition, var_condition]
  expect_identical(view2$X, view_combined$X)
  expect_identical(view2$obs, view_combined$obs)
  expect_identical(view2$var, view_combined$var)
})

test_that("AnnDataView uns field (unaffected by subsetting)", {
  ad <- AnnData(
    X = matrix(1:6, 2L, 3L),
    obs = data.frame(
      row.names = LETTERS[1:2],
      cell_type = c("A", "B")
    ),
    var = data.frame(
      row.names = letters[1:3],
      gene_type = c("X", "Y", "Z")
    ),
    uns = list(
      method = "test",
      params = list(a = 1, b = 2)
    )
  )

  # Create logical condition
  obs_condition <- ad$obs$cell_type == "A"
  view <- ad[obs_condition, ]

  # uns should be unaffected by subsetting
  expect_identical(view$uns, ad$uns)
  expect_equal(view$uns$method, "test")
  expect_equal(view$uns$params$a, 1)
})

test_that("AnnDataView conversion to InMemoryAnnData", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(
      row.names = LETTERS[1:3],
      cell_type = c("A", "B", "A")
    ),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    )
  )

  # Create logical conditions
  obs_condition <- ad$obs$cell_type == "A"
  var_condition <- ad$var$gene_type == "X"

  # Create a view with subsetting
  view <- ad[obs_condition, var_condition]

  # Convert to InMemoryAnnData
  result <- view$as_InMemoryAnnData()

  expect_s3_class(result, "InMemoryAnnData")
  expect_equal(dim(result$X), c(2L, 2L))
  expect_equal(result$obs_names, c("A", "C"))
  expect_equal(result$var_names, c("a", "c"))
  expect_equal(result$obs$cell_type, c("A", "A"))
  expect_equal(result$var$gene_type, c("X", "X"))
})

test_that("AnnDataView error handling", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(row.names = LETTERS[1:3], cell_type = c("A", "B", "A")),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    )
  )

  # Test initialization with non-AbstractAnnData
  expect_error(
    AnnDataView$new("not_an_anndata"),
    "must be an AbstractAnnData object"
  )

  # Test setting values on AnnDataView (create view first)
  # Create logical condition
  obs_condition <- ad$obs$cell_type == "A"

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
  obs_condition <- ad$obs$cell_type == "A"

  view2 <- ad[obs_condition, ]
  expect_s3_class(view2, "AbstractAnnData")
  expect_s3_class(view2, "AnnDataView")
})

test_that("AnnDataView subsetting with invalid conditions", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(row.names = LETTERS[1:3], cell_type = c("A", "B", "A")),
    var = data.frame(
      row.names = letters[1:5],
      gene_type = c("X", "Y", "X", "Y", "Z")
    )
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
})
