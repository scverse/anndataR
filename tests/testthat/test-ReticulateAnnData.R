test_that("ReticulateAnnData creation works", {
  skip_if_no_anndata_py()

  # Test creating empty ReticulateAnnData
  adata <- ReticulateAnnData$new(shape = c(10, 5))
  expect_s3_class(adata, "ReticulateAnnData")
  expect_equal(adata$n_obs(), 10)
  expect_equal(adata$n_vars(), 5)

  # Test creating with data
  X <- matrix(1:20, nrow = 4, ncol = 5)
  obs <- data.frame(row.names = paste0("cell", 1:4), group = letters[1:4])
  var <- data.frame(row.names = paste0("gene", 1:5), type = rep("gene", 5))

  adata2 <- ReticulateAnnData$new(
    X = X,
    obs = obs,
    var = var
  )

  expect_s3_class(adata2, "ReticulateAnnData")
  expect_equal(adata2$n_obs(), 4)
  expect_equal(adata2$n_vars(), 5)
  expect_equal(dim(adata2$X), c(4, 5))
})

test_that("ReticulateAnnData wrapping Python object works", {
  skip_if_no_anndata_py()

  library(reticulate)
  ad <- import("anndata", convert = FALSE)

  # Create a Python AnnData object
  py_x <- r_to_py(matrix(1:20, nrow = 4, ncol = 5))
  py_adata <- ad$AnnData(X = py_x)

  # Wrap it with ReticulateAnnData
  r_adata <- ReticulateAnnData$new(py_anndata = py_adata)

  expect_s3_class(r_adata, "ReticulateAnnData")
  expect_equal(r_adata$n_obs(), 4)
  expect_equal(r_adata$n_vars(), 5)
  expect_equal(dim(r_adata$X), c(4, 5))
})

test_that("ReticulateAnnData slot access works", {
  skip_if_no_anndata_py()

  X <- matrix(1:20, nrow = 4, ncol = 5)
  obs <- data.frame(row.names = paste0("cell", 1:4), group = letters[1:4])
  var <- data.frame(row.names = paste0("gene", 1:5), type = rep("gene", 5))
  layers <- list(
    raw = matrix(21:40, nrow = 4, ncol = 5),
    scaled = matrix(41:60, nrow = 4, ncol = 5)
  )
  uns <- list(
    metadata = "test",
    params = list(param1 = 1, param2 = "value")
  )

  adata <- ReticulateAnnData$new(
    X = X,
    obs = obs,
    var = var,
    layers = layers,
    uns = uns
  )

  # Test X
  # trackstatus: class=ReticulateAnnData, feature=test_get_X, status=done
  expect_equal(dim(adata$X), c(4, 5))

  # Test obs
  # trackstatus: class=ReticulateAnnData, feature=test_get_obs, status=done
  expect_s3_class(adata$obs, "data.frame")
  expect_equal(nrow(adata$obs), 4)

  # Test var
  # trackstatus: class=ReticulateAnnData, feature=test_get_var, status=done
  expect_s3_class(adata$var, "data.frame")
  expect_equal(nrow(adata$var), 5)

  # Test layers
  # trackstatus: class=ReticulateAnnData, feature=test_get_layers, status=done
  expect_type(adata$layers, "list")
  expect_equal(names(adata$layers), c("raw", "scaled"))

  # Test uns
  # trackstatus: class=ReticulateAnnData, feature=test_get_uns, status=done
  expect_type(adata$uns, "list")
  expect_equal(adata$uns$metadata, "test")

  # Test obs_names and var_names
  # trackstatus: class=ReticulateAnnData, feature=test_get_obs_names, status=done
  expect_equal(adata$obs_names, paste0("cell", 1:4))
  # trackstatus: class=ReticulateAnnData, feature=test_get_var_names, status=done
  expect_equal(adata$var_names, paste0("gene", 1:5))
})

test_that("ReticulateAnnData slot setting works", {
  skip_if_no_anndata_py()

  adata <- ReticulateAnnData$new(shape = c(3, 4))

  # Set X
  # trackstatus: class=ReticulateAnnData, feature=test_set_X, status=done
  new_x <- matrix(1:12, nrow = 3, ncol = 4)
  adata$X <- new_x
  expect_equal(dim(adata$X), c(3, 4))

  # Set obs
  # trackstatus: class=ReticulateAnnData, feature=test_set_obs, status=done
  new_obs <- data.frame(
    row.names = paste0("cell", 1:3),
    type = c("A", "B", "C")
  )
  adata$obs <- new_obs
  expect_equal(nrow(adata$obs), 3)
  expect_equal(adata$obs$type, c("A", "B", "C"))

  # Set var
  # trackstatus: class=ReticulateAnnData, feature=test_set_var, status=done
  new_var <- data.frame(
    row.names = paste0("gene", 1:4),
    biotype = rep("protein_coding", 4)
  )
  adata$var <- new_var
  expect_equal(nrow(adata$var), 4)
  expect_equal(adata$var$biotype, rep("protein_coding", 4))

  # Set uns
  # trackstatus: class=ReticulateAnnData, feature=test_set_uns, status=done
  new_uns <- list(experiment = "test", date = "2023-01-01")
  adata$uns <- new_uns
  expect_equal(adata$uns$experiment, "test")
})

test_that("ReticulateAnnData obsm, varm, obsp, varp access works", {
  skip_if_no_anndata_py()

  X <- matrix(1:20, nrow = 4, ncol = 5)
  obs <- data.frame(row.names = paste0("cell", 1:4), group = letters[1:4])
  var <- data.frame(row.names = paste0("gene", 1:5), type = rep("gene", 5))

  obsm <- list(
    X_pca = matrix(rnorm(8), nrow = 4, ncol = 2)
  )
  varm <- list(
    PCs = matrix(rnorm(10), nrow = 5, ncol = 2)
  )
  obsp <- list(
    connectivities = matrix(rnorm(16), nrow = 4, ncol = 4)
  )
  varp <- list(
    correlations = matrix(rnorm(25), nrow = 5, ncol = 5)
  )

  adata <- ReticulateAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obsm = obsm,
    varm = varm,
    obsp = obsp,
    varp = varp
  )

  # Test obsm
  # trackstatus: class=ReticulateAnnData, feature=test_get_obsm, status=done
  expect_type(adata$obsm, "list")
  expect_equal(names(adata$obsm), "X_pca")
  expect_equal(dim(adata$obsm$X_pca), c(4, 2))

  # Test varm
  # trackstatus: class=ReticulateAnnData, feature=test_get_varm, status=done
  expect_type(adata$varm, "list")
  expect_equal(names(adata$varm), "PCs")
  expect_equal(dim(adata$varm$PCs), c(5, 2))

  # Test obsp
  # trackstatus: class=ReticulateAnnData, feature=test_get_obsp, status=done
  expect_type(adata$obsp, "list")
  expect_equal(names(adata$obsp), "connectivities")
  expect_equal(dim(adata$obsp$connectivities), c(4, 4))

  # Test varp
  # trackstatus: class=ReticulateAnnData, feature=test_get_varp, status=done
  expect_type(adata$varp, "list")
  expect_equal(names(adata$varp), "correlations")
  expect_equal(dim(adata$varp$correlations), c(5, 5))
})

test_that("ReticulateAnnData obsm, varm, obsp, varp setting works", {
  skip_if_no_anndata_py()

  X <- matrix(1:20, nrow = 4, ncol = 5)
  adata <- ReticulateAnnData$new(X = X)

  # Set obsm
  # trackstatus: class=ReticulateAnnData, feature=test_set_obsm, status=done
  new_obsm <- list(X_pca = matrix(rnorm(8), nrow = 4, ncol = 2))
  adata$obsm <- new_obsm
  expect_equal(names(adata$obsm), "X_pca")
  expect_equal(dim(adata$obsm$X_pca), c(4, 2))

  # Set varm
  # trackstatus: class=ReticulateAnnData, feature=test_set_varm, status=done
  new_varm <- list(PCs = matrix(rnorm(10), nrow = 5, ncol = 2))
  adata$varm <- new_varm
  expect_equal(names(adata$varm), "PCs")
  expect_equal(dim(adata$varm$PCs), c(5, 2))

  # Set obsp
  # trackstatus: class=ReticulateAnnData, feature=test_set_obsp, status=done
  new_obsp <- list(connectivities = matrix(rnorm(16), nrow = 4, ncol = 4))
  adata$obsp <- new_obsp
  expect_equal(names(adata$obsp), "connectivities")
  expect_equal(dim(adata$obsp$connectivities), c(4, 4))

  # Set varp
  # trackstatus: class=ReticulateAnnData, feature=test_set_varp, status=done
  new_varp <- list(correlations = matrix(rnorm(25), nrow = 5, ncol = 5))
  adata$varp <- new_varp
  expect_equal(names(adata$varp), "correlations")
  expect_equal(dim(adata$varp$correlations), c(5, 5))
})

test_that("ReticulateAnnData names setting works", {
  skip_if_no_anndata_py()

  X <- matrix(1:20, nrow = 4, ncol = 5)
  adata <- ReticulateAnnData$new(X = X)

  # Set obs_names
  # trackstatus: class=ReticulateAnnData, feature=test_set_obs_names, status=done
  new_obs_names <- paste0("sample", 1:4)
  adata$obs_names <- new_obs_names
  expect_equal(adata$obs_names, new_obs_names)

  # Set var_names
  # trackstatus: class=ReticulateAnnData, feature=test_set_var_names, status=done
  new_var_names <- paste0("feature", 1:5)
  adata$var_names <- new_var_names
  expect_equal(adata$var_names, new_var_names)
})

test_that("ReticulateAnnData layers setting works", {
  skip_if_no_anndata_py()

  X <- matrix(1:20, nrow = 4, ncol = 5)
  adata <- ReticulateAnnData$new(X = X)

  # Set layers
  # trackstatus: class=ReticulateAnnData, feature=test_set_layers, status=done
  new_layers <- list(
    raw = matrix(21:40, nrow = 4, ncol = 5),
    normalized = matrix(41:60, nrow = 4, ncol = 5)
  )
  adata$layers <- new_layers
  expect_equal(names(adata$layers), c("raw", "normalized"))
  expect_equal(dim(adata$layers$raw), c(4, 5))
  expect_equal(dim(adata$layers$normalized), c(4, 5))
})

test_that("Conversion between ReticulateAnnData and other types works", {
  skip_if_no_anndata_py()

  X <- matrix(1:20, nrow = 4, ncol = 5)
  obs <- data.frame(row.names = paste0("cell", 1:4), group = letters[1:4])
  var <- data.frame(row.names = paste0("gene", 1:5), type = rep("gene", 5))

  # Create ReticulateAnnData
  ret_adata <- ReticulateAnnData$new(X = X, obs = obs, var = var)

  # Convert to InMemoryAnnData
  mem_adata <- ret_adata$as_InMemoryAnnData()
  expect_s3_class(mem_adata, "InMemoryAnnData")
  expect_equal(mem_adata$n_obs(), 4)
  expect_equal(mem_adata$n_vars(), 5)

  # Convert InMemoryAnnData to ReticulateAnnData
  ret_adata2 <- mem_adata$as_ReticulateAnnData()
  expect_s3_class(ret_adata2, "ReticulateAnnData")
  expect_equal(ret_adata2$n_obs(), 4)
  expect_equal(ret_adata2$n_vars(), 5)
})

test_that("automatic py_to_r and r_to_py conversion works", {
  skip_if_no_anndata_py()

  library(reticulate)
  ad <- import("anndata", convert = FALSE)

  # Create Python AnnData
  py_x <- r_to_py(matrix(1:12, nrow = 3, ncol = 4))
  py_adata <- ad$AnnData(X = py_x)

  # Test automatic py_to_r conversion
  r_adata <- py_to_r(py_adata)
  expect_s3_class(r_adata, "ReticulateAnnData")
  expect_equal(r_adata$n_obs(), 3)
  expect_equal(r_adata$n_vars(), 4)

  # Test automatic r_to_py conversion
  py_adata2 <- r_to_py(r_adata)
  expect_true(inherits(py_adata2, "python.builtin.object"))
  expect_equal(as.integer(py_to_r(py_adata2$n_obs)), 3L)
  expect_equal(as.integer(py_to_r(py_adata2$n_vars)), 4L)

  # Test with other AnnData types
  mem_adata <- AnnData(X = matrix(1:20, nrow = 4, ncol = 5))
  py_adata3 <- r_to_py(mem_adata) # Should convert InMemoryAnnData to Python
  expect_true(inherits(py_adata3, "python.builtin.object"))
  expect_equal(as.integer(py_to_r(py_adata3$n_obs)), 4L)
  expect_equal(as.integer(py_to_r(py_adata3$n_vars)), 5L)
})

test_that("automatic py_to_r and r_to_py conversion works", {
  skip_if_no_anndata_py()

  library(reticulate)
  ad <- import("anndata", convert = FALSE)

  # Create Python AnnData
  py_x <- r_to_py(matrix(1:12, nrow = 3, ncol = 4))
  py_adata <- ad$AnnData(X = py_x)

  # Test automatic Python to R conversion via S3 method
  r_adata <- py_to_r(py_adata)
  expect_s3_class(r_adata, "ReticulateAnnData")
  expect_equal(r_adata$n_obs(), 3)
  expect_equal(r_adata$n_vars(), 4)

  # Test automatic R to Python conversion via S3 method
  py_adata2 <- r_to_py(r_adata)
  expect_true(inherits(py_adata2, "python.builtin.object"))
  expect_equal(as.integer(py_to_r(py_adata2$n_obs)), 3L)
  expect_equal(as.integer(py_to_r(py_adata2$n_vars)), 4L)

  # Test with other AnnData types (should also work automatically)
  mem_adata <- AnnData(X = matrix(1:20, nrow = 4, ncol = 5))
  py_adata3 <- r_to_py(mem_adata)
  expect_true(inherits(py_adata3, "python.builtin.object"))
  expect_equal(as.integer(py_to_r(py_adata3$n_obs)), 4L)
  expect_equal(as.integer(py_to_r(py_adata3$n_vars)), 5L)
})

test_that("reticulate utility functions work", {
  library(reticulate)

  # Test conversion utilities with simple objects
  r_mat <- matrix(1:6, 2, 3)
  py_mat <- r_to_py(r_mat)
  r_mat2 <- py_to_r(py_mat)
  expect_equal(dim(r_mat2), c(2, 3))

  # Test with data frame
  r_df <- data.frame(x = 1:3, y = letters[1:3])
  py_df <- r_to_py(r_df)
  r_df2 <- py_to_r(py_df)
  expect_s3_class(r_df2, "data.frame")
  expect_equal(nrow(r_df2), 3)
})
