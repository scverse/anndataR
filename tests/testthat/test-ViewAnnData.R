test_that("AbstractAnnData basic subsetting and renaming methods", {
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
  
  # Test that basic subsetting methods return ViewAnnData
  view_obs <- ad$subset_obs(cell_type == "A")
  expect_s3_class(view_obs, "ViewAnnData")
  
  view_var <- ad$subset_var(gene_type == "X")
  expect_s3_class(view_var, "ViewAnnData")
  
  # Test that basic renaming methods return ViewAnnData
  view_rename <- ad$rename_obs(cell_category = "cell_type")
  expect_s3_class(view_rename, "ViewAnnData")
})

test_that("ViewAnnData obs subsetting through AbstractAnnData", {
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
  view1 <- ad$subset_obs(c(TRUE, FALSE, TRUE))
  
  expect_equal(nrow(view1$obs), 2L)
  expect_equal(view1$obs_names, c("A", "C"))
  expect_equal(dim(view1$X), c(2L, 5L))
  expect_identical(view1$X, ad$X[c(1, 3), , drop = FALSE])
  
  # Test expression-based subsetting
  view2 <- ad$subset_obs(cell_type == "A")
  
  expect_equal(nrow(view2$obs), 2L)
  expect_equal(view2$obs_names, c("A", "C"))
  expect_equal(view2$obs$cell_type, c("A", "A"))
  
  # Test complex expression
  view3 <- ad$subset_obs(cell_type == "A" & score > 1.6)
  
  expect_equal(nrow(view3$obs), 1L)
  expect_equal(view3$obs_names, "C")
})

test_that("ViewAnnData var subsetting through AbstractAnnData", {
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
  view1 <- ad$subset_var(c(TRUE, FALSE, TRUE, FALSE, FALSE))
  
  expect_equal(nrow(view1$var), 2L)
  expect_equal(view1$var_names, c("a", "c"))
  expect_equal(dim(view1$X), c(3L, 2L))
  expect_identical(view1$X, ad$X[, c(1, 3), drop = FALSE])
  
  # Test expression-based subsetting
  view2 <- ad$subset_var(gene_type == "X")
  
  expect_equal(nrow(view2$var), 2L)
  expect_equal(view2$var_names, c("a", "c"))
  expect_equal(view2$var$gene_type, c("X", "X"))
  
  # Test boolean column subsetting
  view3 <- ad$subset_var(important)
  
  expect_equal(nrow(view3$var), 3L)
  expect_equal(view3$var_names, c("a", "c", "e"))
})

test_that("ViewAnnData combined subsetting through AbstractAnnData", {
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
  
  # Test combined obs and var subsetting
  view <- ad$
    subset_obs(cell_type == "A")$
    subset_var(gene_type == "X")
  
  expect_equal(dim(view$X), c(2L, 2L))
  expect_equal(view$obs_names, c("A", "C"))
  expect_equal(view$var_names, c("a", "c"))
  
  # Check that obsm is subset correctly
  expect_equal(nrow(view$obsm$X_pca), 2L)
  
  # Check that varm is subset correctly
  expect_equal(nrow(view$varm$PCs), 2L)
  
  # Check that obsp is subset correctly (both dimensions)
  expect_equal(dim(view$obsp$connectivities), c(2L, 2L))
  
  # Check that varp is subset correctly (both dimensions)
  expect_equal(dim(view$varp$correlations), c(2L, 2L))
})

test_that("ViewAnnData renaming through AbstractAnnData", {
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
    ),
    layers = list(
      counts = matrix(10:24, 3L, 5L),
      logcounts = matrix(log(10:24), 3L, 5L)
    ),
    obsm = list(
      X_pca = matrix(rnorm(6), 3L, 2L)
    ),
    uns = list(
      method = "test"
    )
  )
  
  # Test obs column renaming
  view1 <- ad$rename_obs(cell_category = "cell_type")
  
  expect_true("cell_category" %in% colnames(view1$obs))
  expect_false("cell_type" %in% colnames(view1$obs))
  expect_equal(view1$obs$cell_category, c("A", "B", "A"))
  
  # Test var column renaming
  view2 <- ad$rename_var(gene_class = "gene_type")
  
  expect_true("gene_class" %in% colnames(view2$var))
  expect_false("gene_type" %in% colnames(view2$var))
  
  # Test layers renaming
  view3 <- ad$rename_layers(raw = "counts", log_norm = "logcounts")
  
  expect_true("raw" %in% names(view3$layers))
  expect_true("log_norm" %in% names(view3$layers))
  expect_false("counts" %in% names(view3$layers))
  expect_false("logcounts" %in% names(view3$layers))
  
  # Test obsm renaming
  view4 <- ad$rename_obsm(pca = "X_pca")
  
  expect_true("pca" %in% names(view4$obsm))
  expect_false("X_pca" %in% names(view4$obsm))
  
  # Test uns renaming
  view5 <- ad$rename_uns(algorithm = "method")
  
  expect_true("algorithm" %in% names(view5$uns))
  expect_false("method" %in% names(view5$uns))
})

test_that("ViewAnnData method chaining through AbstractAnnData", {
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
  
  # Test method chaining
  view <- ad$
    subset_obs(cell_type == "A")$
    subset_var(gene_type == "X")$
    rename_obs(cell_category = "cell_type")$
    rename_var(gene_class = "gene_type")
  
  expect_equal(dim(view$X), c(2L, 2L))
  expect_equal(view$obs_names, c("A", "C"))
  expect_equal(view$var_names, c("a", "c"))
  expect_true("cell_category" %in% colnames(view$obs))
  expect_true("gene_class" %in% colnames(view$var))
  expect_false("cell_type" %in% colnames(view$obs))
  expect_false("gene_type" %in% colnames(view$var))
})

test_that("ViewAnnData transformation functions through AbstractAnnData", {
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
  
  # Test X transformation
  view1 <- ad$add_transformation("X", function(x) x * 2)
  
  expect_identical(view1$X, ad$X * 2)
  
  # Test obs transformation
  view2 <- ad$add_transformation("obs", function(obs) {
    obs$score_scaled <- scale(obs$score)[, 1]
    obs
  })
  
  expect_true("score_scaled" %in% colnames(view2$obs))
  
  # Test multiple transformations by chaining
  view3 <- ad$
    add_transformation("X", function(x) x * 2)$
    add_transformation("X", function(x) x + 1)
  
  expect_identical(view3$X, (ad$X * 2) + 1)
})

test_that("ViewAnnData conversion to InMemoryAnnData through AbstractAnnData", {
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
  
  # Create a view with transformations
  view <- ad$
    subset_obs(cell_type == "A")$
    subset_var(gene_type == "X")$
    rename_obs(cell_category = "cell_type")
  
  # Convert to InMemoryAnnData
  result <- view$as_InMemoryAnnData()
  
  expect_s3_class(result, "InMemoryAnnData")
  expect_equal(dim(result$X), c(2L, 2L))
  expect_equal(result$obs_names, c("A", "C"))
  expect_equal(result$var_names, c("a", "c"))
  expect_true("cell_category" %in% colnames(result$obs))
  expect_false("cell_type" %in% colnames(result$obs))
  expect_equal(result$obs$cell_category, c("A", "A"))
})

test_that("ViewAnnData error handling", {
  ad <- AnnData(
    X = matrix(1:15, 3L, 5L),
    obs = data.frame(row.names = LETTERS[1:3], cell_type = c("A", "B", "A")),
    var = data.frame(row.names = letters[1:5], gene_type = c("X", "Y", "X", "Y", "Z"))
  )
  
  # Test initialization with non-AbstractAnnData
  expect_error(ViewAnnData$new("not_an_anndata"), "must be an AbstractAnnData object")
  
  # Test setting values on ViewAnnData (create view first)
  view <- ad$subset_obs(cell_type == "A")
  expect_error(view$X <- matrix(1:15, 3L, 5L), "Cannot set X on a ViewAnnData object")
  expect_error(view$obs <- data.frame(), "Cannot set obs on a ViewAnnData object")
  
  # Test invalid transformation function
  expect_error(ad$add_transformation("X", "not_a_function"), "transform_func must be a function")
  
  # Test invalid slot name for transformation
  expect_error(ad$add_transformation("invalid_slot", function(x) x), 
               "slot_name must be one of")
  
  # Test that ViewAnnData inherits from AbstractAnnData
  view2 <- ad$subset_obs(cell_type == "A")
  expect_s3_class(view2, "AbstractAnnData")
  expect_s3_class(view2, "ViewAnnData")
})
