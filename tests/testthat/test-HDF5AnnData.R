skip_if_not_installed("rhdf5")
skip_if_not_installed("HDF5Array")

requireNamespace("vctrs")

file <- system.file("extdata", "example.h5ad", package = "anndataR")

test_that("opening H5AD works", {
  adata <- HDF5AnnData$new(file, mode = "r")
  expect_true(inherits(adata, "HDF5AnnData"))
})

test_that("reading an HDF5AnnData at a non-root group works", {
  nested_file <- make_nested_h5ad_fixture("mod/rna")

  nested_adata <- HDF5AnnData$new(nested_file, root = "mod/rna", mode = "r")
  root_adata <- HDF5AnnData$new(file, mode = "r")

  expect_equal(nested_adata$X, root_adata$X)
  expect_equal(nested_adata$layers, root_adata$layers)
  expect_equal(nested_adata$obsm, root_adata$obsm)
  expect_equal(nested_adata$varm, root_adata$varm)
  expect_equal(nested_adata$obsp, root_adata$obsp)
  expect_equal(nested_adata$varp, root_adata$varp)
  expect_equal(nested_adata$obs, root_adata$obs)
  expect_equal(nested_adata$var, root_adata$var)
  expect_equal(nested_adata$obs_names, root_adata$obs_names)
  expect_equal(nested_adata$var_names, root_adata$var_names)
})

test_that("reading a missing root group errors", {
  nested_file <- make_nested_h5ad_fixture("mod/rna")

  expect_error(
    HDF5AnnData$new(nested_file, root = "mod/does_not_exist", mode = "r"),
    "does not exist"
  )
})

test_that("writing a fresh AnnData at a non-root group preserves siblings", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5mu")
  rhdf5::h5createFile(h5ad_file)
  rhdf5::h5createGroup(h5ad_file, "mod")
  rhdf5::h5createGroup(h5ad_file, "mod/other_stuff")
  rhdf5::h5write("hello", h5ad_file, "mod/other_stuff/marker")

  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(
    h5ad_file,
    obs = obs,
    var = var,
    root = "mod/rna",
    mode = "w-"
  )

  expect_identical(h5ad$obs_names, as.character(1:10))
  expect_equal(
    rhdf5::h5read(h5ad_file, "mod/other_stuff/marker"),
    "hello",
    ignore_attr = TRUE
  )
})

test_that("nested intermediate root groups are auto-created", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5mu")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)

  expect_no_error(
    HDF5AnnData$new(
      h5ad_file,
      obs = obs,
      var = var,
      root = "mod/rna",
      mode = "w-"
    )
  )
})

test_that("'w-'/'x' fail if the group already exists, not the file", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5mu")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  HDF5AnnData$new(h5ad_file, obs = obs, var = var, root = "mod/rna", mode = "w-")

  expect_error(
    HDF5AnnData$new(
      h5ad_file,
      obs = obs,
      var = var,
      root = "mod/rna",
      mode = "w-"
    ),
    "already exists"
  )

  # A different group in the same (now-existing) file is unaffected
  expect_no_error(
    HDF5AnnData$new(
      h5ad_file,
      obs = obs,
      var = var,
      root = "mod/atac",
      mode = "w-"
    )
  )
})

test_that("'w' at a sub-root recreates only that group", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5mu")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  HDF5AnnData$new(h5ad_file, obs = obs, var = var, root = "mod/rna", mode = "w-")
  HDF5AnnData$new(h5ad_file, obs = obs, var = var, root = "mod/atac", mode = "w-")

  new_obs <- data.frame(row.names = 1:5)
  new_var <- data.frame(row.names = 1:8)
  h5ad <- HDF5AnnData$new(
    h5ad_file,
    obs = new_obs,
    var = new_var,
    root = "mod/rna",
    mode = "w"
  )

  expect_identical(h5ad$obs_names, as.character(1:5))
  # The sibling modality is untouched
  atac <- HDF5AnnData$new(h5ad_file, root = "mod/atac", mode = "r")
  expect_identical(atac$obs_names, as.character(1:10))
})

test_that("'r+'/'a' warn on a non-empty root group, scoped to the group", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5mu")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  HDF5AnnData$new(h5ad_file, obs = obs, var = var, root = "mod/rna", mode = "w-")

  expect_warning(
    HDF5AnnData$new(h5ad_file, root = "mod/rna", mode = "r+"),
    "non-empty group"
  )
  expect_warning(
    HDF5AnnData$new(h5ad_file, root = "mod/rna", mode = "a"),
    "non-empty group"
  )
})

test_that("round-trip write_h5ad()/read_h5ad() at a non-root group works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5mu")
  adata <- InMemoryAnnData$new(
    X = matrix(rnorm(10 * 20), nrow = 10, ncol = 20),
    obs = data.frame(row.names = paste0("Cell", 1:10)),
    var = data.frame(row.names = paste0("Gene", 1:20))
  )

  write_h5ad(adata, h5ad_file, root = "mod/rna", mode = "a")
  roundtrip <- read_h5ad(
    h5ad_file,
    as = "InMemoryAnnData",
    root = "mod/rna"
  )

  expect_equal(roundtrip$X, adata$X, ignore_attr = TRUE)
  expect_equal(roundtrip$obs_names, adata$obs_names)
  expect_equal(roundtrip$var_names, adata$var_names)
})

adata <- HDF5AnnData$new(file, mode = "r")

# GETTERS ----------------------------------------------------------------
# trackstatus: class=HDF5AnnData, feature=test_get_X, status=done
test_that("reading X works", {
  X <- adata$X
  expect_s4_class(X, "dgRMatrix")
  expect_equal(dim(X), c(50, 100))
})

# trackstatus: class=HDF5AnnData, feature=test_get_layers, status=done
test_that("reading layers works", {
  layers <- adata$layers
  expect_true(is.list(layers), "list")
  expect_equal(
    names(layers),
    c("counts", "csc_counts", "dense_X", "dense_counts")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obsm, status=done
test_that("reading obsm works", {
  obsm <- adata$obsm
  expect_true(is.list(obsm), "list")
  expect_equal(
    names(obsm),
    c("X_pca", "X_umap")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_varm, status=done
test_that("reading varm works", {
  varm <- adata$varm
  expect_true(is.list(varm), "list")
  expect_equal(
    names(varm),
    c("PCs")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obsp, status=done
test_that("reading obsp works", {
  obsp <- adata$obsp
  expect_true(is.list(obsp), "list")
  expect_equal(
    names(obsp),
    c("connectivities", "distances")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_varp, status=done
test_that("reading varp works", {
  varp <- adata$varp
  expect_true(is.list(varp), "list")
  expect_equal(
    names(varp),
    c("test_varp")
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obs, status=done
test_that("reading obs works", {
  obs <- adata$obs
  expect_s3_class(obs, "data.frame")
  expect_equal(
    colnames(obs),
    c(
      "Float",
      "FloatNA",
      "Int",
      "IntNA",
      "Bool",
      "BoolNA",
      "n_genes_by_counts",
      "log1p_n_genes_by_counts",
      "total_counts",
      "log1p_total_counts",
      "leiden"
    )
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_var, status=done
test_that("reading var works", {
  var <- adata$var
  expect_s3_class(var, "data.frame")
  expect_equal(
    colnames(var),
    c(
      "String",
      "n_cells_by_counts",
      "mean_counts",
      "log1p_mean_counts",
      "pct_dropout_by_counts",
      "total_counts",
      "log1p_total_counts",
      "highly_variable",
      "means",
      "dispersions",
      "dispersions_norm"
    )
  )
})

# trackstatus: class=HDF5AnnData, feature=test_get_obs_names, status=done
test_that("reading obs names works", {
  obs_names <- adata$obs_names
  expect_vector(obs_names, ptype = character(), size = 50)
})

# trackstatus: class=HDF5AnnData, feature=test_get_var_names, status=done
test_that("reading var names works", {
  var_names <- adata$var_names
  expect_vector(var_names, ptype = character(), size = 100)
})

# SETTERS ----------------------------------------------------------------
test_that("creating empty H5AD works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  expect_silent(HDF5AnnData$new(file = h5ad_file))
})

# trackstatus: class=HDF5AnnData, feature=test_set_X, status=done
test_that("writing X works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(h5ad$X <- X)
})

# trackstatus: class=HDF5AnnData, feature=test_set_layers, status=done
test_that("writing layers works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  X <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  expect_silent(h5ad$layers <- list(layer1 = X, layer2 = X))
})

# trackstatus: class=HDF5AnnData, feature=test_set_obs, status=done
test_that("writing obs works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  obs <- data.frame(
    Letters = LETTERS[1:10],
    Numbers = 1:10,
    row.names = paste0("Row", 1:10)
  )
  h5ad$obs <- obs
  expect_identical(h5ad$obs_names, paste0("Row", 1:10))
})

# trackstatus: class=HDF5AnnData, feature=test_set_var, status=done
test_that("writing var works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  var <- data.frame(
    Letters = LETTERS[1:20],
    Numbers = 1:20,
    row.names = paste0("Row", 1:20)
  )
  h5ad$var <- var
  expect_identical(h5ad$var_names, paste0("Row", 1:20))
})

# trackstatus: class=HDF5AnnData, feature=test_set_obs_names, status=done
test_that("writing obs names works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  h5ad$obs_names <- LETTERS[1:10]
  expect_identical(h5ad$obs_names, LETTERS[1:10])
})

# trackstatus: class=HDF5AnnData, feature=test_set_var_names, status=done
test_that("writing var names works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  h5ad$var_names <- LETTERS[1:20]
  expect_identical(h5ad$var_names, LETTERS[1:20])
})

# trackstatus: class=HDF5AnnData, feature=test_set_obsm, status=done
test_that("writing obsm works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  obsm_x <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  h5ad$obsm <- list(X = obsm_x)
  # obsm should now have rownames added on-the-fly
  expected_obsm_x <- obsm_x
  rownames(expected_obsm_x) <- h5ad$obs_names
  expect_identical(h5ad$obsm$X, expected_obsm_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_varm, status=done
test_that("writing varm works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)
  varm_x <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  h5ad$varm <- list(PCs = varm_x)
  # varm should now have rownames added on-the-fly
  expected_varm_x <- varm_x
  rownames(expected_varm_x) <- h5ad$var_names
  expect_identical(h5ad$varm$PCs, expected_varm_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_obsp, status=done
test_that("writing obsp works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  obsp_x <- matrix(rnorm(10 * 10), nrow = 10, ncol = 10)
  h5ad$obsp <- list(connectivities = obsp_x)
  # obsp should now have dimnames added on-the-fly
  expected_obsp_x <- obsp_x
  dimnames(expected_obsp_x) <- list(h5ad$obs_names, h5ad$obs_names)
  expect_identical(h5ad$obsp$connectivities, expected_obsp_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_varp, status=done
test_that("writing varp works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  varp_x <- matrix(rnorm(20 * 20), nrow = 20, ncol = 20)
  h5ad$varp <- list(connectivities = varp_x)
  # varp should now have dimnames added on-the-fly
  expected_varp_x <- varp_x
  dimnames(expected_varp_x) <- list(h5ad$var_names, h5ad$var_names)
  expect_identical(h5ad$varp$connectivities, expected_varp_x)
})

# trackstatus: class=HDF5AnnData, feature=test_set_uns, status=done
test_that("writing uns works", {
  h5ad_file <- withr::local_tempfile(fileext = ".h5ad")
  obs <- data.frame(row.names = 1:10)
  var <- data.frame(row.names = 1:20)
  h5ad <- HDF5AnnData$new(h5ad_file, obs = obs, var = var)

  h5ad$uns <- list(
    foo = "bar",
    baz = c(1, 2, 3),
    nested = list(
      nested_foo = "nested_bar",
      nested_baz = c(4L, 5L, 6L)
    )
  )
  expect_identical(h5ad$uns$foo, "bar")
  expect_equal(h5ad$uns$baz, c(1, 2, 3), ignore_attr = TRUE)
  expect_identical(h5ad$uns$nested$nested_foo, "nested_bar")
  expect_equal(h5ad$uns$nested$nested_baz, c(4L, 5L, 6L), ignore_attr = TRUE)
})
