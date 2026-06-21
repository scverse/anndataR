skip_if_not_installed("rhdf5")
skip_if_not_installed("HDF5Array")

requireNamespace("vctrs")

filename <- system.file("extdata", "example.h5ad", package = "anndataR")
hdf5_file <- HDF5File$new(filename)

test_that("reading encoding works", {
  encoding <- read_h5ad_encoding(hdf5_file, "obs")
  expect_equal(names(encoding), c("type", "version"))
})

test_that("reading dense matrices works", {
  mat <- read_h5ad_dense_array(hdf5_file, "layers/dense_counts")
  expect_true(is.matrix(mat))
  expect_type(mat, "integer")
  expect_equal(dim(mat), c(50, 100))

  mat <- read_h5ad_dense_array(hdf5_file, "layers/dense_X")
  expect_true(is.matrix(mat))
  expect_type(mat, "double")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading backed dense matrices works", {
  mat <- read_h5ad_dense_array(hdf5_file, "layers/dense_counts", backed = TRUE)
  expect_s4_class(mat, "DelayedMatrix")
  seed <- DelayedArray::seed(mat)
  expect_identical(DelayedArray::type(seed), "integer")
  expect_false(DelayedArray::is_sparse(seed))
  expect_equal(dim(mat), c(50, 100))

  mat <- read_h5ad_dense_array(hdf5_file, "layers/dense_X", backed = TRUE)
  expect_s4_class(mat, "DelayedMatrix")
  seed <- DelayedArray::seed(mat)
  expect_identical(DelayedArray::type(seed), "double")
  expect_false(DelayedArray::is_sparse(seed))
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading sparse matrices works", {
  mat <- read_h5ad_sparse_array(hdf5_file, "layers/csc_counts", type = "csc")
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(50, 100))

  mat <- read_h5ad_sparse_array(hdf5_file, "layers/counts", type = "csr")
  expect_s4_class(mat, "dgRMatrix")
  expect_equal(dim(mat), c(50, 100))
})

test_that("reading backed sparse matrices works", {
  mat <- read_h5ad_sparse_array(
    hdf5_file,
    "layers/csc_counts",
    type = "csc",
    backed = TRUE
  )
  expect_s4_class(mat, "DelayedMatrix")
  seed <- DelayedArray::seed(mat)
  expect_s4_class(seed, "CSR_H5SparseMatrixSeed")
  expect_identical(DelayedArray::type(seed), "double")
  expect_true(DelayedArray::is_sparse(seed))
  expect_equal(dim(mat), c(50, 100))

  mat <- read_h5ad_sparse_array(
    hdf5_file,
    "layers/counts",
    type = "csr",
    backed = TRUE
  )
  expect_s4_class(mat, "DelayedMatrix")
  seed <- DelayedArray::seed(mat)
  expect_s4_class(seed, "CSC_H5SparseMatrixSeed")
  expect_identical(DelayedArray::type(seed), "double")
  expect_true(DelayedArray::is_sparse(seed))
  expect_equal(dim(mat), c(50, 100))
})

# INVARIANT: backed dense reads must equal an INDEPENDENT eager read
test_that("backed dense reads are value-identical to independent eager reads", {
  for (nm in c("layers/dense_counts", "layers/dense_X")) {
    eager <- read_h5ad_dense_array_base(hdf5_file, nm)
    backed <- read_h5ad_dense_array(hdf5_file, nm, backed = TRUE)
    expect_s4_class(backed, "DelayedMatrix")
    expect_equal(as.matrix(backed), eager)
    raw <- rhdf5::h5read(filename, nm) # on-disk orientation is vars x obs
    expect_equal(as.matrix(backed), t(raw), ignore_attr = TRUE)
  }
})

# INVARIANT: backed sparse reads must equal independent eager reads and stay
# sparse.
test_that("backed sparse reads are value-identical to independent eager reads", {
  cases <- list(
    list(name = "layers/csc_counts", type = "csc_matrix"),
    list(name = "layers/counts", type = "csr_matrix")
  )
  for (case in cases) {
    eager <- read_h5ad_sparse_array_base(hdf5_file, case$name, type = case$type)
    backed <- read_h5ad_sparse_array(
      hdf5_file,
      case$name,
      type = case$type,
      backed = TRUE
    )
    expect_s4_class(backed, "DelayedMatrix")
    expect_true(DelayedArray::is_sparse(backed))
    expect_equal(as.matrix(backed), as.matrix(eager))
  }
})

# INVARIANT: a backed HDF5AnnData exposes the same matrices (values, dimnames,
# orientation) as an eager one, for every matrix-valued slot.
test_that("backed HDF5AnnData slots equal eager slots", {
  eager <- HDF5AnnData$new(filename, mode = "r")
  backed <- HDF5AnnData$new(filename, mode = "r", backed = TRUE)
  withr::defer({
    eager$close()
    backed$close()
  })

  expect_s4_class(backed$X, "DelayedMatrix")
  expect_equal(as.matrix(backed$X), as.matrix(eager$X))
  expect_identical(dimnames(backed$X), dimnames(eager$X))

  for (k in eager$layers_keys()) {
    expect_equal(
      as.matrix(backed$layers[[k]]),
      as.matrix(eager$layers[[k]]),
      info = paste("layer", k)
    )
  }
  for (k in names(eager$obsm)) {
    expect_equal(
      as.matrix(backed$obsm[[k]]),
      as.matrix(eager$obsm[[k]]),
      info = paste("obsm", k)
    )
  }
  for (k in names(eager$obsp)) {
    expect_equal(
      as.matrix(backed$obsp[[k]]),
      as.matrix(eager$obsp[[k]]),
      info = paste("obsp", k)
    )
  }
})

# INVARIANT: a backed matrix reads lazily by file path, so it must keep working
# after the source AnnData is closed and garbage-collected.
test_that("backed matrices stay readable after the source AnnData is closed", {
  eager <- HDF5AnnData$new(filename, mode = "r")
  expected <- as.matrix(eager$X)
  eager$close()

  backed <- HDF5AnnData$new(filename, mode = "r", backed = TRUE)
  x_backed <- backed$X
  expect_s4_class(x_backed, "DelayedMatrix")
  backed$close()
  rm(backed)
  gc()

  expect_equal(as.matrix(x_backed), expected)
})

# Regression test: subsetting a backed AnnData must stay lazy AND apply the
# subset.
test_that("subsetting a backed AnnData stays lazy and subsets correctly", {
  backed <- HDF5AnnData$new(filename, mode = "r", backed = TRUE)
  eager <- HDF5AnnData$new(filename, mode = "r")
  withr::defer({
    backed$close()
    eager$close()
  })

  oi <- c(1, 3, 5, 7, 9)
  vi <- 1:10
  vb <- backed[oi, vi]
  ve <- eager[oi, vi]

  expect_s4_class(vb$X, "DelayedMatrix")
  expect_equal(dim(vb$X), c(length(oi), length(vi)))
  expect_equal(as.matrix(vb$X), as.matrix(ve$X))
})

# Regression test: converting a backed AnnData to InMemory must materialize its
# DelayedArrays into ordinary in-memory matrices
test_that("converting a backed AnnData to InMemory materializes its matrices", {
  backed <- HDF5AnnData$new(filename, mode = "r", backed = TRUE)
  withr::defer(backed$close())
  im <- backed$as_InMemoryAnnData()

  expect_false(inherits(im$X, "DelayedArray"))
  expect_equal(as.matrix(im$X), as.matrix(backed$X))
})

# INVARIANT: a backed SingleCellExperiment keeps its assays lazy (DelayedMatrix)
# and value-equal to an eager conversion.
test_that("backed SingleCellExperiment has lazy assays equal to eager ones", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))

  backed <- read_h5ad(filename, as = "HDF5AnnData", backed = TRUE)
  eager <- read_h5ad(filename, as = "HDF5AnnData")
  withr::defer({
    backed$close()
    eager$close()
  })

  sce_b <- backed$as_SingleCellExperiment()
  sce_e <- eager$as_SingleCellExperiment()
  for (a in SummarizedExperiment::assayNames(sce_e)) {
    expect_s4_class(SummarizedExperiment::assay(sce_b, a), "DelayedMatrix")
    expect_equal(
      as.matrix(SummarizedExperiment::assay(sce_b, a)),
      as.matrix(SummarizedExperiment::assay(sce_e, a)),
      info = paste("assay", a)
    )
  }
})

# Regression test (issue #387, reported by LouiseDck): converting a backed SCE
# back to AnnData must preserve shape.
test_that("round-trip of a backed SCE back to AnnData preserves shape", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))

  backed <- read_h5ad(filename, as = "HDF5AnnData", backed = TRUE)
  withr::defer(backed$close())
  sce_backed <- backed$as_SingleCellExperiment()

  ad2 <- as_AnnData(sce_backed)
  expect_equal(ad2$shape(), backed$shape())
})

test_that("reading recarrays works", {
  array_list <- read_h5ad_rec_array(
    hdf5_file,
    "uns/rank_genes_groups/logfoldchanges"
  )
  expect_true(is.list(array_list))
  expect_equal(names(array_list), c("0", "1", "2", "3", "4", "5"))
  for (array in array_list) {
    expect_true(is.vector(array))
    expect_type(array, "double")
    expect_equal(length(array), 100)
  }
})

test_that("reading 1D numeric arrays works", {
  array_1d <- read_h5ad_dense_array(hdf5_file, "obs/Int")
  expect_equal(array_1d, array(0L:49L))

  array_1d <- read_h5ad_dense_array(hdf5_file, "obs/Float")
  expect_equal(array_1d, array(rep(42.42, 50)))
})

test_that("reading 1D sparse numeric arrays works", {
  array_1d <- read_h5ad_sparse_array(hdf5_file, "uns/Sparse1D", type = "csc")
  expect_s4_class(array_1d, "dgCMatrix")
  expect_equal(dim(array_1d), c(1, 6))
})

test_that("reading 1D nullable arrays works", {
  array_1d <- read_h5ad_nullable_integer(hdf5_file, "obs/IntNA")
  expect_vector(array_1d, ptype = integer(), size = 50)
  expect_true(anyNA(array_1d))

  array_1d <- read_h5ad_dense_array(hdf5_file, "obs/FloatNA")
  expected <- array(rep(42.42, 50))
  expected[1] <- NA
  expect_equal(array_1d, expected)

  array_1d <- read_h5ad_nullable_boolean(hdf5_file, "obs/Bool")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_false(anyNA(array_1d))

  array_1d <- read_h5ad_nullable_boolean(hdf5_file, "obs/BoolNA")
  expect_vector(array_1d, ptype = logical(), size = 50)
  expect_true(anyNA(array_1d))
})

test_that("reading string scalars works", {
  scalar <- read_h5ad_string_scalar(hdf5_file, "uns/StringScalar")
  expect_equal(scalar, "A string")
})

test_that("reading numeric scalars works", {
  scalar <- read_h5ad_numeric_scalar(hdf5_file, "uns/IntScalar")
  expect_equal(scalar, 1)
})

test_that("reading string arrays works", {
  array <- read_h5ad_string_array(hdf5_file, "uns/String")
  expect_equal(array, array(paste0("String ", 0L:9L)))

  array <- read_h5ad_string_array(hdf5_file, "uns/String2D")
  expect_true(is.matrix(array))
  expect_type(array, "character")
  expect_equal(dim(array), c(5, 10))
})

test_that("reading mappings works", {
  mapping <- read_h5ad_mapping(hdf5_file, "uns")
  expect_type(mapping, "list")
  expect_type(names(mapping), "character")
})

test_that("reading mapping keys works", {
  keys <- read_h5ad_mapping_keys(hdf5_file, "uns")
  expect_type(keys, "character")
  expect_true(length(keys) > 0)
})

test_that("reading dataframes works", {
  df <- read_h5ad_data_frame(hdf5_file, "obs")
  expect_s3_class(df, "data.frame")
  expect_equal(
    colnames(df),
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

test_that("reading dataframe keys works", {
  keys <- read_h5ad_data_frame_keys(hdf5_file, "obs")

  expect_type(keys$cols, "character")
  expect_identical(
    keys$cols,
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

  expect_type(keys$rows, "character")
  expect_length(keys$rows, 50)
})

test_that("reading H5AD as SingleCellExperiment works", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))

  sce <- read_h5ad(filename, as = "SingleCellExperiment")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("reading H5AD as Seurat works", {
  skip_if_not_installed("Seurat")

  seurat <- read_h5ad(filename, as = "Seurat")
  expect_s4_class(seurat, "Seurat")
})
