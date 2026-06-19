skip_if_not_installed("Rarr")

for (zarr_format in c(2, 3)) {
  zarr_zip <- system.file(
    "extdata",
    paste0("example_v", zarr_format, ".zarr.zip"),
    package = "anndataR"
  )
  td <- tempdir(check = TRUE)
  unzip(zarr_zip, exdir = td)
  store <- file.path(td, paste0("example_v", zarr_format, ".zarr"))

  test_that(paste("reading Zarr", zarr_format, "encoding works"), {
    encoding <- read_zarr_encoding(store, "obs")
    expect_equal(names(encoding), c("type", "version"))
  })

  test_that(paste("reading Zarr", zarr_format, "dense matrices works"), {
    mat <- read_zarr_dense_array(store, "layers/dense_counts")
    expect_true(is.matrix(mat))
    expect_type(mat, "integer")
    expect_equal(dim(mat), c(50, 100))

    mat <- read_zarr_dense_array(store, "layers/dense_X")
    expect_true(is.matrix(mat))
    expect_type(mat, "double")
    expect_equal(dim(mat), c(50, 100))
  })

  test_that(paste("reading Zarr", zarr_format, "sparse matrices works"), {
    mat <- read_zarr_sparse_array(store, "layers/csc_counts", type = "csc")
    expect_s4_class(mat, "dgCMatrix")
    expect_equal(dim(mat), c(50, 100))

    mat <- read_zarr_sparse_array(store, "layers/counts", type = "csr")
    expect_s4_class(mat, "dgRMatrix")
    expect_equal(dim(mat), c(50, 100))
  })

  test_that(paste("reading Zarr", zarr_format, "recarrays works"), {
    array_list <- read_zarr_rec_array(
      store,
      "uns/rank_genes_groups/logfoldchanges"
    )
    expect_true(is.list(array_list))
    for (array in array_list) {
      expect_true(is.vector(array))
      # Rarr 2.1.0 introduced a breaking change: structured datatypes (record arrays)
      # now return lists as their internal elements instead of vectors. The
      # Bioconductor release (older Rarr) and devel (>= 2.1.0) therefore return
      # recarrays in different shapes, so the relevant assertions branch on this flag.
      # See https://github.com/scverse/anndataR/issues/409
      if (packageVersion("Rarr") >= "2.1.0") {
        # Rarr >= 2.1.0 (Bioc devel) returns each field as a list of scalars
        expect_type(array, "list")
        expect_type(unlist(array), "double")
      } else {
        # Older Rarr (Bioc release) returns each field as a double vector
        expect_type(array, "double")
      }
      expect_equal(length(array), 6)
    }
  })

  test_that(paste("reading Zarr", zarr_format, "1D numeric arrays works"), {
    array_1d <- read_zarr_dense_array(store, "obs/Int")
    expect_equal(array_1d, array(0L:49L))

    array_1d <- read_zarr_dense_array(store, "obs/Float")
    expect_equal(array_1d, array(rep(42.42, 50)))
  })

  test_that(
    paste("reading Zarr", zarr_format, "1D sparse numeric arrays works"),
    {
      array_1d <- read_zarr_sparse_array(store, "uns/Sparse1D", type = "csc")
      expect_s4_class(array_1d, "dgCMatrix")
      expect_equal(dim(array_1d), c(1, 6))
    }
  )

  test_that(paste("reading Zarr", zarr_format, "1D nullable arrays works"), {
    array_1d <- read_zarr_nullable_integer(store, "obs/IntNA")
    expect_vector(array_1d, ptype = integer(), size = 50)
    expect_true(anyNA(array_1d))

    array_1d <- read_zarr_dense_array(store, "obs/FloatNA")
    expected <- array(rep(42.42, 50))
    expected[1] <- NA
    expect_equal(array_1d, expected)

    array_1d <- read_zarr_nullable_boolean(store, "obs/BoolNA")
    expect_vector(array_1d, ptype = logical(), size = 50)
    expect_true(anyNA(array_1d))
  })

  test_that(paste("reading Zarr", zarr_format, "string scalars works"), {
    scalar <- read_zarr_string_scalar(store, "uns/StringScalar")
    expect_equal(scalar, "A string")
  })

  test_that(paste("reading Zarr", zarr_format, "numeric scalars works"), {
    scalar <- read_zarr_numeric_scalar(store, "uns/IntScalar")
    expect_equal(scalar, 1)
  })

  test_that(paste("reading Zarr", zarr_format, "string arrays works"), {
    array <- read_zarr_string_array(store, "uns/String")
    expect_equal(array, array(paste0("String ", 0L:9L)))

    array <- read_zarr_string_array(store, "uns/String2D")
    expect_true(is.matrix(array))
    expect_type(array, "character")
    expect_equal(dim(array), c(5, 10))
  })

  test_that(paste("reading Zarr", zarr_format, "mappings works"), {
    if (zarr_format == "3") {
      # TODO: Remove when v3 recarray support is implemented
      mapping <- suppressWarnings(read_zarr_mapping(store, "uns"))
    } else {
      mapping <- read_zarr_mapping(store, "uns")
    }
    expect_type(mapping, "list")
    expect_type(names(mapping), "character")
  })

  test_that(paste("reading Zarr", zarr_format, "dataframes works"), {
    df <- read_zarr_data_frame(store, "obs")
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

  test_that(
    paste("reading Zarr", zarr_format, "as SingleCellExperiment works"),
    {
      skip_if_not_installed("SingleCellExperiment")

      if (zarr_format == "v3") {
        # TODO: Remove when v3 recarray support is implemented
        sce <- suppressWarnings(read_zarr(store, as = "SingleCellExperiment"))
      } else {
        sce <- read_zarr(store, as = "SingleCellExperiment")
      }

      expect_s4_class(sce, "SingleCellExperiment")
    }
  )

  test_that(paste("reading Zarr", zarr_format, "as Seurat works"), {
    skip_if_not_installed("SeuratObject")

    if (zarr_format == "v3") {
      # TODO: Remove when v3 recarray support is implemented
      seurat <- suppressWarnings(read_zarr(store, as = "Seurat"))
    } else {
      seurat <- read_zarr(store, as = "Seurat")
    }

    expect_s4_class(seurat, "Seurat")
  })
}
