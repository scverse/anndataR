skip_if_not_installed("SparseArray")
skip_if_not_installed("BiocGenerics")
skip_if_not_installed("MatrixGenerics")

# A matrix with more than 2^31 - 1 non-zero values needs tens of gigabytes of
# memory, so these tests pretend a dgCMatrix can only hold a few values. This
# runs the real code path on a small matrix. See test-large-sparse-real.R for
# the same checks on a matrix that really is too large.
local_small_nnz_limit <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    sparse_max_nnz = function() 10,
    sparse_block_nnz = function() 7,
    .package = "anndataR",
    .env = env
  )
}

n_obs <- 12L
n_vars <- 8L
obs_names <- paste0("cell", seq_len(n_obs))
var_names <- paste0("gene", seq_len(n_vars))

dense <- local({
  set.seed(1)
  mat <- as.matrix(Matrix::rsparsematrix(n_obs, n_vars, density = 0.4))
  # An empty row and column should not get lost between blocks
  mat[3L, ] <- 0
  mat[, 2L] <- 0
  mat
})
csc <- as(dense, "CsparseMatrix")
csr <- as(dense, "RsparseMatrix")
svt <- as(csc, "SVT_SparseMatrix")

new_adata <- function(X = NULL, layers = NULL) {
  AnnData(
    X = X,
    layers = layers,
    obs = data.frame(row.names = obs_names),
    var = data.frame(row.names = var_names)
  )
}

expect_same_matrix <- function(object, expected = dense) {
  object <- as.matrix(object)
  dimnames(object) <- NULL
  testthat::expect_equal(object, expected)
}

backends <- list(
  h5ad = list(
    ext = ".h5ad",
    write = function(adata, path) write_h5ad(adata, path, mode = "w"),
    read = function(path, ...) read_h5ad(path, ...)
  ),
  zarr_v2 = list(
    ext = ".zarr",
    write = function(adata, path) write_zarr(adata, path, zarr_format = 2L),
    read = function(path, ...) read_zarr(path, ...)
  ),
  zarr_v3 = list(
    ext = ".zarr",
    write = function(adata, path) write_zarr(adata, path, zarr_format = 3L),
    read = function(path, ...) read_zarr(path, ...)
  )
)

skip_if_backend_unavailable <- function(backend) {
  if (backend == "h5ad") {
    testthat::skip_if_not_installed("rhdf5")
  } else {
    testthat::skip_if_not_installed("Rarr", minimum_version = "2.1.45")
  }
}

test_that("sparse_block_bounds() covers every major element exactly once", {
  local_small_nnz_limit()

  # 3 elements of 3 values, an empty one, then a single one above the block size
  pointers <- c(0, 3, 6, 9, 9, 29)
  bounds <- sparse_block_bounds(pointers)

  expect_equal(
    unlist(lapply(bounds, function(b) b[[1L]]:b[[2L]])),
    seq_len(length(pointers) - 1L)
  )
  sizes <- vapply(
    bounds,
    function(b) pointers[[b[[2L]] + 1L]] - pointers[[b[[1L]]]],
    numeric(1L)
  )
  # Only a block of one element can be larger than the block size
  expect_true(all(sizes[sizes > 7] == 20))
  expect_identical(sparse_block_bounds(0), list())
})

for (backend in names(backends)) {
  be <- backends[[backend]]

  test_that(paste(backend, "reads too many non-zero values as SVT"), {
    skip_if_backend_unavailable(backend)
    path <- withr::local_tempfile(fileext = be$ext)
    be$write(new_adata(X = csr, layers = list(csc = csc)), path)

    # Below the limit nothing changes
    adata <- be$read(path)
    expect_s4_class(adata$X, "dgRMatrix")
    expect_s4_class(adata$layers$csc, "dgCMatrix")

    local_small_nnz_limit()
    adata <- be$read(path)
    expect_s4_class(adata$X, "SVT_SparseMatrix")
    expect_s4_class(adata$layers$csc, "SVT_SparseMatrix")
    expect_same_matrix(adata$X)
    expect_same_matrix(adata$layers$csc)
    expect_identical(rownames(adata$X), obs_names)
    expect_identical(colnames(adata$X), var_names)
  })

  test_that(paste(backend, "writes SVT matrices as CSC"), {
    skip_if_backend_unavailable(backend)
    path <- withr::local_tempfile(fileext = be$ext)
    # Without the limit this is a plain dgCMatrix write
    be$write(new_adata(X = svt), path)
    adata <- be$read(path)
    expect_s4_class(adata$X, "dgCMatrix")
    expect_same_matrix(adata$X)

    path <- withr::local_tempfile(fileext = be$ext)
    local({
      local_small_nnz_limit()
      be$write(new_adata(X = svt, layers = list(svt = svt)), path)
    })
    adata <- be$read(path)
    expect_s4_class(adata$X, "dgCMatrix")
    expect_same_matrix(adata$X)
    expect_same_matrix(adata$layers$svt)
  })

  test_that(paste(backend, "SVT writes keep NA and all-zero matrices"), {
    skip_if_backend_unavailable(backend)
    with_na <- dense
    with_na[1L, 1L] <- NA_real_

    path <- withr::local_tempfile(fileext = be$ext)
    local({
      local_small_nnz_limit()
      be$write(
        new_adata(X = as(as(with_na, "CsparseMatrix"), "SVT_SparseMatrix")),
        path
      )
    })
    expect_same_matrix(be$read(path)$X, with_na)
  })

  test_that(paste(backend, "writes SVT matrices of other types as double"), {
    skip_if_backend_unavailable(backend)
    counts <- dense
    counts[counts != 0] <- seq_len(sum(counts != 0))
    storage.mode(counts) <- "integer"

    for (typed in list(counts, dense != 0)) {
      path <- withr::local_tempfile(fileext = be$ext)
      local({
        local_small_nnz_limit()
        be$write(new_adata(X = as(typed, "SVT_SparseMatrix")), path)
      })
      expect_same_matrix(be$read(path)$X, typed + 0)
    }
  })
}

test_that("large H5AD indptr is stored as 64-bit integers", {
  skip_if_not_installed("rhdf5")
  path <- withr::local_tempfile(fileext = ".h5ad")
  local({
    local_small_nnz_limit()
    write_h5ad(new_adata(X = svt), path, mode = "w")
  })

  listing <- rhdf5::h5ls(path, all = TRUE)
  listing <- listing[listing$group == "/X", ]
  rownames(listing) <- listing$name
  expect_identical(listing["indptr", "dtype"], "H5T_STD_I64LE")
  expect_identical(listing["indices", "dtype"], "H5T_STD_I32LE")
  expect_identical(listing["data", "dtype"], "H5T_IEEE_F64LE")
  expect_identical(
    rhdf5::h5readAttributes(path, "X")[["encoding-type"]],
    "csc_matrix"
  )
})

test_that("reading a large H5AD sorts indices and checks them", {
  skip_if_not_installed("rhdf5")
  path <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(new_adata(X = csr), path, mode = "w")

  # Swap the first two indices and values of the first row
  indices <- as.vector(rhdf5::h5read(path, "X/indices"))
  values <- as.vector(rhdf5::h5read(path, "X/data"))
  rhdf5::h5write(
    indices[c(2L, 1L, seq_along(indices)[-(1:2)])],
    path,
    "X/indices"
  )
  rhdf5::h5write(values[c(2L, 1L, seq_along(values)[-(1:2)])], path, "X/data")
  rhdf5::h5closeAll()

  local_small_nnz_limit()
  expect_same_matrix(read_h5ad(path)$X)

  # Indices that don't fit the shape aren't silently accepted
  rhdf5::h5write(indices + 100L, path, "X/indices")
  rhdf5::h5closeAll()
  expect_error(
    read_h5ad_element(HDF5File$new(path), "X", stop_on_error = TRUE),
    "does not fit shape"
  )
})

test_that("a missing SparseArray gives an informative error", {
  skip_if_not_installed("rhdf5")
  path <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(new_adata(X = csr), path, mode = "w")

  local_small_nnz_limit()
  local_mocked_bindings(
    check_requires = function(what, requires, ...) {
      cli::cli_abort("{what} requires {requires}")
    }
  )
  expect_error(
    read_h5ad_element(HDF5File$new(path), "X", stop_on_error = TRUE),
    "requires SparseArray"
  )
})

test_that("py_to_r_matrix() looks at nnz without overflowing", {
  skip_if_no_anndata_py()
  bi <- reticulate::import_builtins(convert = FALSE)
  # Documents why nnz is read as a float. If this starts failing, reticulate no
  # longer overflows and bi$float() in py_to_r_matrix() can go.
  expect_lt(reticulate::py_to_r(bi$int(2.2e9)), 0)
  expect_identical(reticulate::py_to_r(bi$float(bi$int(2.2e9))), 2.2e9)
})

test_that("64-bit Zarr integers round trip above the 32-bit range", {
  skip_if_not_installed("Rarr", minimum_version = "2.1.45")
  values <- c(0, 1, 2^31 - 1, 2^31, 2^31 + 2^32, 3e9, 5e9, 2^40 + 7, 2^52 + 1)

  for (zarr_format in 2:3) {
    store <- withr::local_tempdir()
    Rarr::write_zarr_group(store, "", zarr_version = zarr_format)
    zarr_write_int64(store, "x", values, "zstd", zarr_format)
    expect_identical(
      as.double(Rarr::read_zarr_array(file.path(store, "x"))),
      values
    )
  }
})

test_that("64-bit Zarr integers are compatible with Python", {
  skip_if_not_installed("Rarr", minimum_version = "2.1.45")
  skip_if_no_zarr()
  zarr <- reticulate::import("zarr", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  everything <- reticulate::import_builtins(convert = FALSE)$slice(NULL)
  values <- c(0, 2^31 - 1, 2^31, 3e9, 5e9, 2^40 + 7)
  py_values <- np$asarray(reticulate::np_array(values))$astype("int64")

  for (zarr_format in 2:3) {
    # Written by R, read by Python
    store <- withr::local_tempdir()
    Rarr::write_zarr_group(store, "", zarr_version = zarr_format)
    zarr_write_int64(store, "x", values, "zstd", zarr_format)
    array <- zarr$open_array(file.path(store, "x"), mode = "r")
    expect_identical(reticulate::py_to_r(array$dtype$name), "int64")
    expect_equal(
      reticulate::py_to_r(array$`__getitem__`(everything)$astype("float64")),
      values,
      ignore_attr = TRUE
    )

    # Written by Python, read by R
    path <- file.path(withr::local_tempdir(), "x")
    array <- zarr$create_array(
      path,
      shape = reticulate::tuple(length(values)),
      chunks = reticulate::tuple(2L),
      dtype = "int64",
      zarr_format = zarr_format
    )
    reticulate::py_set_item(array, everything, py_values)
    expect_identical(as.double(Rarr::read_zarr_array(path)), values)
  }

  # Sharded arrays only exist in Zarr v3
  path <- file.path(withr::local_tempdir(), "x")
  array <- zarr$create_array(
    path,
    shape = reticulate::tuple(length(values)),
    chunks = reticulate::tuple(2L),
    shards = reticulate::tuple(4L),
    dtype = "int64",
    zarr_format = 3L
  )
  reticulate::py_set_item(array, everything, py_values)
  expect_identical(as.double(Rarr::read_zarr_array(path)), values)
})

# A sparse matrix with two major elements and 2^31 + 10 non-zero values. Only
# the indptr has values, data and indices have no chunks, so this takes no
# memory or disk space. That is enough to see whether a backed read gets past
# the indptr.
write_large_zarr_sparse <- function(store, type, zarr_format) {
  nnz <- 2^31 + 10
  Rarr::write_zarr_group(store, "", zarr_version = zarr_format)
  Rarr::write_zarr_group(store, "X", zarr_version = zarr_format)
  zarr_write_int64(store, "X/indptr", c(0, 2^30, nnz), "zstd", zarr_format)
  arrays <- c(data = "double", indices = "integer")
  for (array in names(arrays)) {
    Rarr::create_empty_zarr_array(
      file.path(store, "X", array),
      dim = nnz,
      chunk_dim = 2^22,
      data_type = arrays[[array]],
      fill_value = 0,
      zarr_version = zarr_format
    )
  }
  shape <- if (type == "csr_matrix") {
    list(2L, .Machine$integer.max)
  } else {
    list(.Machine$integer.max, 2L)
  }
  Rarr::write_zarr_attributes(
    file.path(store, "X"),
    new.zattrs = list(
      `encoding-type` = type,
      `encoding-version` = "0.1.0",
      shape = shape
    ),
    zarr_version = zarr_format
  )
}

for (type in c("csr_matrix", "csc_matrix")) {
  test_that(paste("backed Zarr reads a large", type), {
    skip_if_not_installed("ZarrArray")

    for (zarr_format in 2:3) {
      store <- withr::local_tempfile(fileext = ".zarr")
      write_large_zarr_sparse(store, type, zarr_format)
      X <- read_zarr_sparse_array(store, "X", backed = TRUE, type = type)
      expect_s4_class(X, "DelayedMatrix")
    }
  })
}

test_that("backed Zarr reads sparse matrices below the limit as before", {
  skip_if_not_installed("ZarrArray")
  store <- withr::local_tempfile(fileext = ".zarr")
  write_zarr(new_adata(X = csr), store)
  expect_same_matrix(read_zarr(store, as = "ZarrAnnData", backed = TRUE)$X)
})

test_that("large SVT writes can be read by Python", {
  skip_if_no_anndata_py()
  ad <- reticulate::import("anndata", convert = FALSE)

  skip_if_not_installed("rhdf5")
  path <- withr::local_tempfile(fileext = ".h5ad")
  local({
    local_small_nnz_limit()
    write_h5ad(new_adata(X = svt), path, mode = "w")
  })
  expect_equal(
    reticulate::py_to_r(ad$read_h5ad(path)$X$toarray()),
    dense,
    ignore_attr = TRUE
  )

  skip_if_not_installed("Rarr", minimum_version = "2.1.45")
  skip_if_no_zarr()
  for (zarr_format in 2:3) {
    path <- withr::local_tempfile(fileext = ".zarr")
    local({
      local_small_nnz_limit()
      write_zarr(new_adata(X = svt), path, zarr_format = zarr_format)
    })
    expect_equal(
      reticulate::py_to_r(ad$read_zarr(path)$X$toarray()),
      dense,
      ignore_attr = TRUE
    )
  }
})

test_that("Python-written large matrices with 64-bit indices can be read", {
  skip_if_no_anndata_py()
  skip_if_not_installed("rhdf5")
  ad <- reticulate::import("anndata", convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  # scipy uses 64-bit indices and indptr for matrices that don't fit in 32 bits
  py_csr <- sp$csr_matrix(reticulate::r_to_py(dense))
  py_csr$indices <- py_csr$indices$astype(np$int64)
  py_csr$indptr <- py_csr$indptr$astype(np$int64)
  py_adata <- ad$AnnData(X = py_csr)

  path <- withr::local_tempfile(fileext = ".h5ad")
  py_adata$write_h5ad(path)
  expect_identical(
    rhdf5::h5ls(path)[rhdf5::h5ls(path)$name == "indptr", "dclass"],
    "INTEGER"
  )

  local_small_nnz_limit()
  expect_same_matrix(read_h5ad(path)$X)

  skip_if_not_installed("Rarr", minimum_version = "2.1.45")
  skip_if_no_zarr()
  path <- withr::local_tempfile(fileext = ".zarr")
  py_adata$write_zarr(path)
  expect_same_matrix(read_zarr(path)$X)
})

test_that("InMemoryAnnData, HDF5AnnData and ZarrAnnData hold SVT matrices", {
  adata <- new_adata(X = svt, layers = list(svt = svt))
  expect_s4_class(adata$X, "SVT_SparseMatrix")
  expect_identical(rownames(adata$X), obs_names)
  expect_same_matrix(adata$layers$svt)

  # Dimnames are checked as for any other matrix
  named <- svt
  rownames(named) <- rev(obs_names)
  expect_error(adata$X <- named, "is not as expected")

  skip_if_not_installed("rhdf5")
  local_small_nnz_limit()
  path <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(new_adata(X = csr), path, mode = "w")
  h5ad <- suppressWarnings(read_h5ad(path, as = "HDF5AnnData", mode = "r+"))
  h5ad$layers <- list(svt = svt)
  expect_same_matrix(h5ad$layers$svt)

  memory <- h5ad$as_InMemoryAnnData()
  expect_s4_class(memory$layers$svt, "SVT_SparseMatrix")
})

test_that("AnnDataView subsets SVT matrices", {
  adata <- new_adata(X = svt, layers = list(svt = svt))
  view <- adata[2:5, c(1L, 3L, 8L)]

  expect_equal(dim(view$X), c(4L, 3L))
  expect_same_matrix(view$X, dense[2:5, c(1L, 3L, 8L)])
  expect_same_matrix(view$layers$svt, dense[2:5, c(1L, 3L, 8L)])
  expect_same_matrix(
    view$as_InMemoryAnnData()$X,
    dense[2:5, c(1L, 3L, 8L)]
  )
})

test_that("to_py_matrix() and to_R_matrix() don't densify", {
  expect_s4_class(to_py_matrix(svt), "SVT_SparseMatrix")
  expect_s4_class(to_R_matrix(svt), "SVT_SparseMatrix")
  expect_same_matrix(to_py_matrix(svt), t(dense))
  expect_same_matrix(to_R_matrix(svt), t(dense))

  skip_if_not_installed("DelayedArray")
  delayed <- DelayedArray::DelayedArray(svt)
  expect_s4_class(to_py_matrix(delayed), "DelayedMatrix")
  expect_same_matrix(to_py_matrix(delayed), t(dense))
})

test_that("is_numeric_matrix() doesn't read a backed matrix", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")

  path <- withr::local_tempfile(fileext = ".h5ad")
  adata <- new_adata()
  adata$obsm <- list(embedding = dense)
  write_h5ad(adata, path, mode = "w")
  embedding <- read_h5ad(path, as = "HDF5AnnData", backed = TRUE)$obsm$embedding
  expect_s4_class(embedding, "DelayedMatrix")

  local({
    local_mocked_bindings(as.matrix = function(...) stop("realized the matrix"))
    expect_true(is_numeric_matrix(embedding))
    expect_true(is_numeric_matrix(svt))
  })

  expect_true(is_numeric_matrix(svt))
  expect_true(is_numeric_matrix(dense))
  expect_true(is_numeric_matrix(csr))
  expect_true(is_numeric_matrix(data.frame(a = 1:2, b = 3:4)))
  expect_false(is_numeric_matrix(data.frame(a = c("x", "y"))))
  expect_false(is_numeric_matrix(matrix("a", 2L, 2L)))
  skip_if_not_installed("DelayedArray")
  expect_false(is_numeric_matrix(DelayedArray::DelayedArray(matrix(
    "a",
    2L,
    2L
  ))))
})

test_that("materialize_matrix() converts SVT to a CsparseMatrix", {
  expect_s4_class(materialize_matrix(svt), "CsparseMatrix")
})

test_that("SVT matrices convert to and from SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  adata <- new_adata(X = svt, layers = list(svt = svt))

  expect_no_warning(sce <- as_SingleCellExperiment(adata))
  expect_identical(dim(sce), c(n_vars, n_obs))
  assay <- SummarizedExperiment::assay(sce, "X")
  expect_s4_class(assay, "SVT_SparseMatrix")
  expect_same_matrix(assay, t(dense))
  expect_same_matrix(SummarizedExperiment::assay(sce, "svt"), t(dense))

  expect_no_warning(back <- as_AnnData(sce, x_mapping = "X"))
  expect_s4_class(back$X, "SVT_SparseMatrix")
  expect_same_matrix(back$X)
  expect_same_matrix(back$layers$svt)
})

test_that("SVT obsp and varp convert to SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  square <- as(
    as(Matrix::rsparsematrix(n_obs, n_obs, density = 0.2), "CsparseMatrix"),
    "SVT_SparseMatrix"
  )
  adata <- new_adata(X = svt)
  adata$obsp <- list(graph = square)

  sce <- as_SingleCellExperiment(adata)
  expect_identical(names(SingleCellExperiment::colPairs(sce)), "graph")
})

test_that("as_Seurat() gives a dgCMatrix when SVT matrices fit", {
  skip_if_not_installed("Seurat")
  adata <- new_adata(X = svt)

  expect_no_message(seurat <- as_Seurat(adata, x_mapping = "counts"))
  layer <- SeuratObject::LayerData(seurat, "counts")
  expect_s4_class(layer, "dgCMatrix")
  expect_same_matrix(layer, t(dense))
})

test_that("as_Seurat() gives a DelayedMatrix for too many non-zero values", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("DelayedArray")
  local_small_nnz_limit()
  adata <- new_adata(X = svt, layers = list(svt = svt))

  expect_message(
    seurat <- as_Seurat(adata, x_mapping = "counts"),
    "DelayedMatrix"
  )
  layer <- SeuratObject::LayerData(seurat, "counts")
  expect_s4_class(layer, "DelayedMatrix")
  expect_same_matrix(layer, t(dense))
  expect_same_matrix(SeuratObject::LayerData(seurat, "svt"), t(dense))
})

test_that("from_Seurat() keeps a DelayedMatrix layer as it is", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("DelayedArray")
  local_mocked_bindings(
    # as.matrix() on a layer that doesn't fit in memory is what should not happen
    to_py_matrix = local({
      original <- to_py_matrix
      function(mat) {
        result <- original(mat)
        expect_s4_class(result, "DelayedMatrix")
        result
      }
    })
  )
  layer <- DelayedArray::DelayedArray(as(t(csc), "SVT_SparseMatrix"))
  dimnames(layer) <- list(var_names, obs_names)
  seurat <- SeuratObject::CreateSeuratObject(
    SeuratObject::CreateAssay5Object(counts = layer)
  )

  adata <- as_AnnData(seurat, x_mapping = "counts")
  expect_s4_class(adata$X, "DelayedMatrix")
  expect_same_matrix(adata$X)
})

test_that("backed matrices can be written to H5AD and Zarr", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("Rarr", minimum_version = "2.1.45")

  path <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(new_adata(X = csr, layers = list(csc = csc)), path, mode = "w")
  backed <- read_h5ad(path, as = "HDF5AnnData", backed = TRUE)
  expect_s4_class(backed$X, "DelayedMatrix")

  for (limit in c("none", "small")) {
    local({
      if (limit == "small") {
        local_small_nnz_limit()
      }
      h5ad_out <- withr::local_tempfile(fileext = ".h5ad")
      zarr_out <- withr::local_tempfile(fileext = ".zarr")
      write_h5ad(backed, h5ad_out, mode = "w")
      write_zarr(backed, zarr_out)
    })
  }

  zarr_out <- withr::local_tempfile(fileext = ".zarr")
  write_zarr(backed, zarr_out)
  adata <- read_zarr(zarr_out)
  expect_same_matrix(adata$X)
  expect_same_matrix(adata$layers$csc)
})

test_that("ReticulateAnnData converts large sparse matrices in both directions", {
  skip_if_no_anndata_py()
  ad <- reticulate::import("anndata", convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  py_adata <- ad$AnnData(
    X = sp$csr_matrix(reticulate::r_to_py(dense)),
    layers = reticulate::dict(csc = sp$csc_matrix(reticulate::r_to_py(dense))),
    obsp = reticulate::dict(
      graph = sp$csr_matrix(reticulate::r_to_py(dense %*% t(dense)))
    )
  )

  local_small_nnz_limit()

  adata <- reticulate::py_to_r(py_adata)
  expect_s4_class(adata$X, "SVT_SparseMatrix")
  expect_s4_class(adata$layers$csc, "SVT_SparseMatrix")
  expect_same_matrix(adata$X)
  expect_same_matrix(adata$layers$csc)
  expect_same_matrix(adata$obsp$graph, dense %*% t(dense))

  adata$X <- svt
  adata$layers <- list(svt = svt)
  expect_same_matrix(
    reticulate::py_to_r(adata$py_anndata()$X$toarray())
  )
  expect_same_matrix(
    reticulate::py_to_r(adata$py_anndata()$layers[["svt"]]$toarray())
  )
  expect_identical(
    reticulate::py_to_r(adata$py_anndata()$X$format),
    "csc"
  )

  memory <- new_adata(X = svt)
  expect_same_matrix(
    reticulate::py_to_r(reticulate::r_to_py(memory)$X$toarray())
  )
})

test_that("reticulate conversion of SVT and DelayedMatrix below the limit", {
  skip_if_no_anndata_py()
  py_svt <- reticulate::r_to_py(svt)
  expect_same_matrix(reticulate::py_to_r(py_svt$toarray()))

  skip_if_not_installed("DelayedArray")
  delayed <- DelayedArray::DelayedArray(svt)
  expect_same_matrix(reticulate::py_to_r(reticulate::r_to_py(
    delayed
  )$toarray()))
  dense_delayed <- DelayedArray::DelayedArray(dense)
  expect_same_matrix(reticulate::py_to_r(reticulate::r_to_py(dense_delayed)))
})
