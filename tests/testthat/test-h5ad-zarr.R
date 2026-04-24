skip_if_not_installed("rhdf5")
skip_if_not_installed("Rarr")

# h5ad file
filename <- system.file("extdata", "example.h5ad", package = "anndataR")
file <- rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY", native = FALSE)

# zarr file
zarr_dir <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(zarr_dir, exdir = td)
store <- file.path(td, "example_v2.zarr")

# helper to compare h5ad and zarr reads for the same path
expect_equal_h5ad_zarr <- function(h5ad_fn, zarr_fn, path, ...) {
  testthat::expect_equal(h5ad_fn(file, path, ...), zarr_fn(store, path, ...))
}

# compare rec arrays of h5ad and zarr
compare_rec_array <- function(rec_array_h5ad, rec_array_zarr, test_fun) {
  test_fun(length(rec_array_h5ad), length(rec_array_zarr[[1]]))
  test_fun(do.call(rbind, rec_array_h5ad), {
    array_list_zarr_mat <- do.call(cbind, rec_array_zarr)
    rownames(array_list_zarr_mat) <-
      paste(0:(nrow(array_list_zarr_mat) - 1))
    array_list_zarr_mat
  })
}

test_that("reading dense matrices is the same for h5ad and zarr", {
  for (path in c("layers/dense_counts", "layers/dense_X")) {
    expect_equal_h5ad_zarr(read_h5ad_dense_array, read_zarr_dense_array, path)
  }
})

test_that("reading sparse matrices is same for h5ad and zarr", {
  sparse_mats <- list(
    list(path = "layers/csc_counts", type = "csc"),
    list(path = "layers/counts", type = "csr")
  )
  for (mat in sparse_mats) {
    expect_equal_h5ad_zarr(
      read_h5ad_sparse_array,
      read_zarr_sparse_array,
      mat$path,
      type = mat$type
    )
  }
})

test_that("reading recarrays is the the same for h5ad and zarr", {
  # h5ad returns a list of 6 arrays of length 100
  array_list_h5ad <- read_h5ad_rec_array(
    file,
    "uns/rank_genes_groups/logfoldchanges"
  )
  # zarr returns a list of 100 arrays of length 6
  array_list_zarr <- read_zarr_rec_array(
    store,
    "uns/rank_genes_groups/logfoldchanges"
  )
  compare_rec_array(array_list_h5ad, array_list_zarr, expect_equal)
})

test_that("reading 1D numeric arrays is the same for h5ad and zarr", {
  for (path in c("obs/Int", "obs/Float")) {
    expect_equal_h5ad_zarr(read_h5ad_dense_array, read_zarr_dense_array, path)
  }
})

test_that("reading 1D sparse numeric arrays is the same for h5ad and zarr", {
  expect_equal_h5ad_zarr(
    read_h5ad_sparse_array,
    read_zarr_sparse_array,
    "uns/Sparse1D",
    type = "csc"
  )
})

test_that("reading 1D nullable arrays is the same for h5ad and zarr", {
  expect_equal_h5ad_zarr(
    read_h5ad_nullable_integer,
    read_zarr_nullable_integer,
    "obs/IntNA"
  )
  expect_equal_h5ad_zarr(
    read_h5ad_dense_array,
    read_zarr_dense_array,
    "obs/FloatNA"
  )
  for (path in c("obs/Bool", "obs/BoolNA")) {
    expect_equal_h5ad_zarr(
      read_h5ad_nullable_boolean,
      read_zarr_nullable_boolean,
      path
    )
  }
})

test_that("reading string scalars is the same for h5ad and zarr", {
  expect_equal_h5ad_zarr(
    read_h5ad_string_scalar,
    read_zarr_string_scalar,
    "uns/StringScalar"
  )
})

test_that("reading numeric scalars is the same for h5ad and zarr", {
  expect_equal_h5ad_zarr(
    read_h5ad_numeric_scalar,
    read_zarr_numeric_scalar,
    "uns/IntScalar"
  )
})

test_that("reading string arrays is the same for h5ad and zarr", {
  for (path in c("uns/String", "uns/String2D")) {
    expect_equal_h5ad_zarr(read_h5ad_string_array, read_zarr_string_array, path)
  }
})

# TODO: Re-enable when recarays are handled consistently, see https://github.com/scverse/anndataR/issues/409
test_that("reading mappings is the same for h5ad and zarr", {
  skip(
    "skipping test for mappings since rec arrays are read differently
       across h5ad and zarr"
  )
  # since rec arrays are read differently across h5ad and zarr,
  # we compare all elements individually
  mapping_h5ad <- read_h5ad_mapping(file, "uns")
  mapping_zarr <- read_zarr_mapping(store, "uns")
  for (nm in names(mapping_h5ad)) {
    if (!nm %in% "rank_genes_groups") {
      expect_equal(mapping_h5ad[[nm]], mapping_zarr[[nm]])
    } else {
      map_ranks_h5ad <- mapping_h5ad$rank_genes_groups
      map_ranks_zarr <- mapping_zarr$rank_genes_groups
      lapply(
        names(map_ranks_h5ad)[!names(map_ranks_h5ad) %in% "params"],
        function(nmr) {
          print(nmr)
          compare_rec_array(
            map_ranks_h5ad[[nmr]],
            map_ranks_zarr[[nmr]],
            expect_equal
          )
        }
      )
    }
  }
})

tmp <- read_zarr_element(store, "uns/neighbors/params/random_state")
tmp2 <- read_h5ad_element(file, "uns/neighbors/params/random_state")

test_that("reading dataframes is the the same for h5ad and zarr", {
  expect_equal_h5ad_zarr(read_h5ad_data_frame, read_zarr_data_frame, "obs")
})

rhdf5::H5Fclose(file)

test_that("reading H5AD as SingleCellExperiment is the same for h5ad and zarr", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("S4Vectors")
  sce_h5ad <- read_h5ad(filename, as = "SingleCellExperiment")
  sce_zarr <- read_zarr(store, as = "SingleCellExperiment")
  # TODO: Update when recarays are handled consistently, see https://github.com/scverse/anndataR/issues/409
  S4Vectors::metadata(sce_zarr) <- S4Vectors::metadata(sce_h5ad)
  expect_equal(sce_h5ad, sce_zarr)
})

test_that("reading H5AD as Seurat is the same for h5ad and zarr", {
  skip_if_not_installed("Seurat")
  sce_h5ad <- read_h5ad(filename, as = "Seurat")
  sce_zarr <- read_zarr(store, as = "Seurat")
  # TODO: Update when recarays are handled consistently, see https://github.com/scverse/anndataR/issues/409
  expect_warning(
    Seurat::Misc(sce_zarr, "rank_genes_groups") <-
      Seurat::Misc(sce_h5ad, "rank_genes_groups"),
    "Overwriting miscellanous"
  )
  # TODO: neighbors/params/random_state and
  # leiden/params/random_state read as 0 in Python anndata but as an empty
  # array in Zarr
  expect_warning(
    Seurat::Misc(sce_zarr, "neighbors") <-
      Seurat::Misc(sce_h5ad, "neighbors"),
    "Overwriting miscellanous"
  )
  expect_warning(
    Seurat::Misc(sce_zarr, "leiden") <-
      Seurat::Misc(sce_h5ad, "leiden"),
    "Overwriting miscellanous"
  )
  # Sort Misc by name to make comparison order-agnostic
  sce_h5ad@misc <- sce_h5ad@misc[sort(names(sce_h5ad@misc))]
  sce_zarr@misc <- sce_zarr@misc[sort(names(sce_zarr@misc))]
  expect_equal(sce_h5ad, sce_zarr)
})
