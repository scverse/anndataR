skip_if_not_installed("Rarr", minimum_version = "2.1.45")
skip_if_not_installed("ZarrArray", minimum_version = "1.1.5")

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

  test_that(paste("reading backed Zarr", zarr_format, "dense matrices works"), {
    mat <- read_zarr_dense_array(store, "layers/dense_counts", backed = TRUE)
    expect_s4_class(mat, "DelayedMatrix")
    seed <- DelayedArray::seed(mat)
    expect_identical(DelayedArray::type(seed), "integer")
    expect_false(DelayedArray::is_sparse(seed))
    expect_equal(dim(mat), c(50, 100))

    mat <- read_zarr_dense_array(store, "layers/dense_X", backed = TRUE)
    expect_s4_class(mat, "DelayedMatrix")
    seed <- DelayedArray::seed(mat)
    expect_identical(DelayedArray::type(seed), "double")
    expect_false(DelayedArray::is_sparse(seed))
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

  test_that(
    paste("reading backed Zarr", zarr_format, "sparse matrices works"),
    {
      mat <- read_zarr_sparse_array(
        store,
        "layers/csc_counts",
        type = "csc",
        backed = TRUE
      )
      expect_s4_class(mat, "DelayedMatrix")
      seed <- DelayedArray::seed(mat)
      expect_s4_class(seed, "CSR_ZarrSparseMatrixSeed")
      expect_identical(DelayedArray::type(seed), "double")
      expect_true(DelayedArray::is_sparse(seed))
      expect_equal(dim(mat), c(50, 100))

      mat <- read_zarr_sparse_array(
        store,
        "layers/counts",
        type = "csr",
        backed = TRUE
      )
      expect_s4_class(mat, "DelayedMatrix")
      seed <- DelayedArray::seed(mat)
      expect_s4_class(seed, "CSC_ZarrSparseMatrixSeed")
      expect_identical(DelayedArray::type(seed), "double")
      expect_true(DelayedArray::is_sparse(seed))
      expect_equal(dim(mat), c(50, 100))
    }
  )

  # INVARIANT: backed dense reads must equal an INDEPENDENT eager read
  test_that("backed dense reads are value-identical to independent eager reads", {
    for (nm in c("layers/dense_counts", "layers/dense_X")) {
      eager <- read_zarr_dense_array_base(store, nm)
      backed <- read_zarr_dense_array(store, nm, backed = TRUE)
      expect_s4_class(backed, "DelayedMatrix")
      expect_equal(as.matrix(backed), eager)
      raw <- Rarr::read_zarr_array(file.path(store, nm))
      expect_equal(as.matrix(backed), raw, ignore_attr = TRUE)
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
      eager <- read_zarr_sparse_array_base(store, case$name, type = case$type)
      backed <- read_zarr_sparse_array(
        store,
        case$name,
        type = case$type,
        backed = TRUE
      )
      expect_s4_class(backed, "DelayedMatrix")
      expect_true(DelayedArray::is_sparse(backed))
      expect_equal(as.matrix(backed), as.matrix(eager))
    }
  })

  # INVARIANT: a backed ZarrAnnData exposes the same matrices (values, dimnames,
  # orientation) as an eager one, for every matrix-valued slot.
  test_that("backed ZarrAnnData slots equal eager slots", {
    eager <- ZarrAnnData$new(store, mode = "r")
    backed <- ZarrAnnData$new(store, mode = "r", backed = TRUE)

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
    eager <- ZarrAnnData$new(store, mode = "r")
    expected <- as.matrix(eager$X)

    backed <- ZarrAnnData$new(store, mode = "r", backed = TRUE)
    x_backed <- backed$X
    expect_s4_class(x_backed, "DelayedMatrix")
    rm(backed)
    gc()

    expect_equal(as.matrix(x_backed), expected)
  })

  # Regression test: subsetting a backed AnnData must stay lazy AND apply the
  # subset.
  test_that("subsetting a backed AnnData stays lazy and subsets correctly", {
    backed <- ZarrAnnData$new(store, mode = "r", backed = TRUE)
    eager <- ZarrAnnData$new(store, mode = "r")

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
    backed <- ZarrAnnData$new(store, mode = "r", backed = TRUE)
    im <- backed$as_InMemoryAnnData()

    expect_false(inherits(im$X, "DelayedArray"))
    expect_equal(as.matrix(im$X), as.matrix(backed$X))
  })

  # INVARIANT: a backed SingleCellExperiment keeps its assays lazy (DelayedMatrix)
  # and value-equal to an eager conversion.
  test_that("backed SingleCellExperiment has lazy assays equal to eager ones", {
    suppressWarnings(skip_if_not_installed("SingleCellExperiment"))

    backed <- read_zarr(store, as = "ZarrAnnData", backed = TRUE)
    eager <- read_zarr(store, as = "ZarrAnnData")

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

  test_that("backed Seurat graphs equal the obsp of an eager AnnData", {
    skip_if_not_installed("Seurat")

    backed <- read_zarr(store, as = "ZarrAnnData", backed = TRUE)
    eager <- read_zarr(store, as = "ZarrAnnData")

    seurat <- backed$as_Seurat()
    for (graph in eager$obsp_keys()) {
      expect_equal(
        as.matrix(seurat[[paste0("RNA_", graph)]]),
        as.matrix(eager$obsp[[graph]]),
        ignore_attr = TRUE,
        info = paste("graph", graph)
      )
    }
  })

  # the graphs in example_v*.zarr are all csr, so write one of each format
  write_graphs_zarr <- function(store, zarr_format) {
    obsp <- Matrix::sparseMatrix(
      i = c(1, 2, 3, 5),
      j = c(2, 3, 1, 6),
      x = c(1, 2, 3, 4),
      dims = c(6, 6)
    )
    varp <- Matrix::sparseMatrix(
      i = c(1, 2, 3),
      j = c(2, 3, 4),
      x = c(1, 2, 3),
      dims = c(4, 4)
    )
    adata <- AnnData(
      X = matrix(as.numeric(1:24), 6, 4),
      obs = data.frame(row.names = paste0("cell", 1:6)),
      var = data.frame(row.names = paste0("gene", 1:4)),
      layers = list(counts = matrix(as.numeric(1:24), 6, 4)),
      obsp = list(
        csc = obsp,
        csr = as(obsp, "RsparseMatrix"),
        dense = as.matrix(obsp)
      ),
      varp = list(
        csc = varp,
        csr = as(varp, "RsparseMatrix"),
        dense = as.matrix(varp)
      )
    )
    write_zarr(adata, store, zarr_format = zarr_format)
  }

  test_that("backed Seurat graphs work for every format", {
    skip_if_not_installed("Seurat")

    store <- withr::local_tempfile(fileext = ".zarr")
    write_graphs_zarr(store, zarr_format)

    backed <- read_zarr(store, as = "ZarrAnnData", backed = TRUE)
    eager <- read_zarr(store, as = "ZarrAnnData")

    seurat <- backed$as_Seurat()
    for (graph in eager$obsp_keys()) {
      expect_equal(
        as.matrix(seurat[[paste0("RNA_", graph)]]),
        as.matrix(eager$obsp[[graph]]),
        ignore_attr = TRUE,
        info = paste("graph", graph)
      )
    }
  })

  test_that("backed SingleCellExperiment pairs work for every format", {
    suppressWarnings(skip_if_not_installed("SingleCellExperiment"))

    store <- withr::local_tempfile(fileext = ".zarr")
    write_graphs_zarr(store, zarr_format)

    backed <- read_zarr(store, as = "ZarrAnnData", backed = TRUE)
    eager <- read_zarr(store, as = "ZarrAnnData")

    sce <- backed$as_SingleCellExperiment()
    for (graph in eager$obsp_keys()) {
      expect_equal(
        as.matrix(SingleCellExperiment::colPair(sce, graph, asSparse = TRUE)),
        as.matrix(eager$obsp[[graph]]),
        ignore_attr = TRUE,
        info = paste("colPair", graph)
      )
    }
    for (graph in eager$varp_keys()) {
      expect_equal(
        as.matrix(SingleCellExperiment::rowPair(sce, graph, asSparse = TRUE)),
        as.matrix(eager$varp[[graph]]),
        ignore_attr = TRUE,
        info = paste("rowPair", graph)
      )
    }
  })

  # Regression test (issue #387, reported by LouiseDck): converting a backed SCE
  # back to AnnData must preserve shape.
  test_that("round-trip of a backed SCE back to AnnData preserves shape", {
    suppressWarnings(skip_if_not_installed("SingleCellExperiment"))

    backed <- read_zarr(store, as = "ZarrAnnData", backed = TRUE)
    sce_backed <- backed$as_SingleCellExperiment()

    ad2 <- as_AnnData(sce_backed)
    expect_equal(ad2$shape(), backed$shape())
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

  test_that(
    paste("reading Zarr", zarr_format, "nullable string arrays works"),
    {
      array <- read_zarr_nullable_string(store, "uns/StringNA")
      expect_vector(array, ptype = character(), size = 10)
      expect_true(anyNA(array))
    }
  )

  test_that(paste("reading Zarr", zarr_format, "mappings works"), {
    mapping <- read_zarr_mapping(store, "uns")
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
      sce <- read_zarr(store, as = "SingleCellExperiment")
      expect_s4_class(sce, "SingleCellExperiment")
    }
  )

  test_that(paste("reading Zarr", zarr_format, "as Seurat works"), {
    skip_if_not_installed("SeuratObject")
    seurat <- read_zarr(store, as = "Seurat")
    expect_s4_class(seurat, "Seurat")
  })

  test_that(paste("reading Zarr", zarr_format, "as InMemoryAnnData works"), {
    adata <- read_zarr(store, as = "InMemoryAnnData", mode = "r")
    expect_equal(class(adata), c("InMemoryAnnData", "AbstractAnnData", "R6"))
  })

  test_that(paste("reading Zarr", zarr_format, " as ZarrAnnData works"), {
    adata <- read_zarr(store, as = "ZarrAnnData", mode = "r")
    expect_equal(class(adata), c("ZarrAnnData", "AbstractAnnData", "R6"))
  })

  test_that(
    paste("reading Zarr", zarr_format, " as backed ZarrAnnData works"),
    {
      adata <- read_zarr(store, as = "ZarrAnnData", backed = TRUE)

      expect_equal(class(adata), c("ZarrAnnData", "AbstractAnnData", "R6"))

      expect_s4_class(adata$X, "DelayedArray")

      for (layer in adata$layers) {
        expect_s4_class(layer, "DelayedArray")
      }

      for (obsm in adata$obsm) {
        expect_s4_class(obsm, "DelayedArray")
      }

      for (varm in adata$varm) {
        expect_s4_class(varm, "DelayedArray")
      }

      for (obsp in adata$obsp) {
        expect_s4_class(obsp, "DelayedArray")
      }

      for (varp in adata$varp) {
        expect_s4_class(varp, "DelayedArray")
      }
    }
  )

  test_that(
    paste("reading Zarr", zarr_format, " as backed SingleCellExperiment works"),
    {
      skip_if_not_installed("SingleCellExperiment")

      sce <- read_zarr(store, as = "SingleCellExperiment", backed = TRUE)

      expect_s4_class(sce, "SingleCellExperiment")

      for (assay_name in SummarizedExperiment::assayNames(sce)) {
        expect_s4_class(
          SummarizedExperiment::assay(sce, assay_name),
          "DelayedArray"
        )
      }

      for (reduced_dim_name in SingleCellExperiment::reducedDimNames(sce)) {
        reduced_dim <- SingleCellExperiment::reducedDim(sce, reduced_dim_name)

        if (inherits(reduced_dim, "LinearEmbeddingMatrix")) {
          expect_s4_class(
            SingleCellExperiment::sampleFactors(reduced_dim),
            "DelayedArray"
          )
          expect_s4_class(
            SingleCellExperiment::featureLoadings(reduced_dim),
            "DelayedArray"
          )
        } else {
          expect_s4_class(reduced_dim, "DelayedArray")
        }
      }

      for (row_pair_name in SingleCellExperiment::rowPairNames(sce)) {
        # SingleCellExperiment does not support DelayedArray for rowPairs
        expect_s4_class(
          SingleCellExperiment::rowPair(sce, row_pair_name),
          "SelfHits"
        )
      }

      for (col_pair_name in SingleCellExperiment::colPairNames(sce)) {
        # SingleCellExperiment does not support DelayedArray for colPairs
        expect_s4_class(
          SingleCellExperiment::colPair(sce, col_pair_name),
          "SelfHits"
        )
      }
    }
  )

  test_that(paste("reading Zarr", zarr_format, " as backed Seurat works"), {
    skip_if_not_installed("SeuratObject")

    seurat <- read_zarr(store, as = "Seurat", backed = TRUE)

    expect_s4_class(seurat, "Seurat")

    for (layer in SeuratObject::Layers(seurat)) {
      expect_s4_class(SeuratObject::LayerData(seurat, layer), "DelayedArray")
    }

    for (reduction in SeuratObject::Reductions(seurat)) {
      # Seurat does not support DelayedArray for embeddings and loadings
      expect_true(is.matrix(SeuratObject::Embeddings(seurat, reduction)))
      expect_true(is.matrix(SeuratObject::Loadings(seurat, reduction)))
    }

    for (graph in SeuratObject::Graphs(seurat)) {
      # Seurat does not support DelayedArray for graphs
      expect_true(inherits(seurat[[graph]], "dgCMatrix"))
    }
  })
}

# unzip the Zarr v3 example to a temp dir that is removed when `env` exits
local_example_v3_store <- function(env = parent.frame()) {
  zarr_zip <- system.file(
    "extdata",
    "example_v3.zarr.zip",
    package = "anndataR"
  )
  td <- withr::local_tempdir(.local_envir = env)
  unzip(zarr_zip, exdir = td)
  file.path(td, "example_v3.zarr")
}

test_that("the Zarr v3 example store is sharded", {
  store <- local_example_v3_store()

  # X is a sparse group, so check the codecs of its data array
  metadata <- jsonlite::read_json(file.path(store, "X", "data", "zarr.json"))
  codecs <- vapply(metadata$codecs, function(codec) codec$name, character(1))
  expect_contains(codecs, "sharding_indexed")
})

test_that("the Zarr v3 example reads like the H5AD example", {
  skip_if_not_installed("rhdf5")
  store <- local_example_v3_store()
  h5ad <- read_h5ad(system.file(
    "extdata",
    "example.h5ad",
    package = "anndataR"
  ))
  zarr <- read_zarr(store)

  for (slot in c("X", "obs", "var", "obs_names", "var_names")) {
    expect_equal(zarr[[slot]], h5ad[[slot]], label = slot)
  }
  # the order of the keys differs between H5AD and Zarr, so compare them by name
  for (slot in c("layers", "obsm", "varm", "obsp", "varp", "uns")) {
    expect_setequal(names(zarr[[slot]]), names(h5ad[[slot]]))
    keys <- setdiff(names(h5ad[[slot]]), "rank_genes_groups")
    expect_gt(length(keys), 0)
    expect_equal(zarr[[slot]][keys], h5ad[[slot]][keys], label = slot)
  }
  # TODO: compare uns/rank_genes_groups when recarrays are read consistently, see
  # https://github.com/scverse/anndataR/issues/409
  expect_equal(
    zarr$uns$rank_genes_groups$params,
    h5ad$uns$rank_genes_groups$params
  )
})

test_that("errors while reading a Zarr element keep the original message", {
  zarr_zip <- system.file(
    "extdata",
    "example_v3.zarr.zip",
    package = "anndataR"
  )
  td <- withr::local_tempdir()
  unzip(zarr_zip, exdir = td)
  store <- file.path(td, "example_v3.zarr")

  local_mocked_bindings(
    read_zarr_dense_array = function(...) stop("value not in {0,...,n}"),
    read_zarr_mapping_keys = function(...) stop("value not in {0,...,n}")
  )
  expect_warning(
    read_zarr_element(store, "layers/dense_counts"),
    "value not in {0,...,n}",
    fixed = TRUE
  )
  expect_error(
    read_zarr_element(store, "layers/dense_counts", stop_on_error = TRUE),
    "value not in {0,...,n}",
    fixed = TRUE
  )
  expect_warning(
    read_zarr_element_keys(store, "uns"),
    "value not in {0,...,n}",
    fixed = TRUE
  )
  expect_error(
    read_zarr_element_keys(store, "uns", stop_on_error = TRUE),
    "value not in {0,...,n}",
    fixed = TRUE
  )
})
