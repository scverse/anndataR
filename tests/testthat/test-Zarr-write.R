skip_if_not_installed("Rarr", minimum_version = "2.1.27")

for (zarr_format in c(2, 3)) {
  # set version
  withr::with_options(
    list("anndataR.zarr_format" = zarr_format),
    {
      test_that("Zarr versioning works", {
        expect_equal(
          getOption("anndataR.zarr_format"),
          zarr_format
        )
      })

      store <- tempfile(fileext = ".zarr")
      if (dir.exists(store)) {
        unlink(store, recursive = TRUE)
      }

      Rarr::write_zarr_group(store, "", zarr_version = zarr_format)
      metafile <- if (zarr_format == 3) "zarr.json" else ".zgroup"

      test_that("Writing Zarr store works", {
        expect_equal(zarr_node_format(store, ""), zarr_format)
      })

      test_that("Writing Zarr store with different version works", {
        store_w_diff <- tempfile(fileext = ".zarr")
        if (dir.exists(store_w_diff)) {
          unlink(store_w_diff, recursive = TRUE)
        }
        zarr_format_diff <- setdiff(c(2, 3), zarr_format)
        Rarr::write_zarr_group(
          store_w_diff,
          "",
          zarr_version = zarr_format_diff
        )
        expect_equal(zarr_node_format(store_w_diff, ""), zarr_format_diff)
      })

      test_that("Writing a Zarr store at a relative path works", {
        withr::with_tempdir({
          ZarrAnnData$new(
            "relative.zarr",
            obs = data.frame(row.names = 1:3),
            var = data.frame(row.names = 1:2),
            mode = "w"
          )
          expect_true(dir.exists("relative.zarr"))
          expect_true(file.exists(file.path("relative.zarr", metafile)))
          expect_equal(zarr_node_format("relative.zarr", ""), zarr_format)
        })
      })

      test_that("Writing a Zarr store with regex characters in the name works", {
        store_regex <- tempfile(pattern = "a+b", fileext = ".zarr")
        ZarrAnnData$new(
          store_regex,
          obs = data.frame(row.names = 1:3),
          var = data.frame(row.names = 1:2),
          mode = "w"
        )
        expect_true(file.exists(file.path(store_regex, metafile)))
        expect_equal(zarr_node_format(store_regex, ""), zarr_format)
        # The store must not be nested inside a directory of the same name
        expect_false(dir.exists(file.path(store_regex, basename(store_regex))))
      })

      test_that("Writing Zarr dense arrays works", {
        array <- matrix(rnorm(20), nrow = 5, ncol = 4)

        expect_silent(write_zarr_element(
          array,
          store,
          "dense_array",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "dense_array"))
        expect_equal(zarr_node_format(store, "dense_array"), zarr_format)
        expect_equal(
          {
            zarr_ov <- Rarr::zarr_overview(
              file.path(store, "dense_array"),
              as_data_frame = TRUE
            )
            zarr_ov$compressor
          },
          "zstd"
        )
        attrs <- Rarr::read_zarr_attributes(file.path(store, "dense_array"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "array")
      })

      test_that("Writing Zarr dense 3D arrays works", {
        value <- array(rnorm(60), dim = c(5, 4, 3))

        expect_silent(
          write_zarr_element(
            value,
            store,
            "dense_3d_array",
            zarr_format = zarr_format
          )
        )
        expect_true(zarr_path_exists(store, "dense_3d_array"))
        expect_equal(zarr_node_format(store, "dense_3d_array"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "dense_3d_array"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "array")
      })

      test_that("Writing Zarr sparse arrays works", {
        array <- matrix(rnorm(20), nrow = 5, ncol = 4)

        csc_array <- as(array, "CsparseMatrix")
        expect_silent(write_zarr_element(
          csc_array,
          store,
          "csc_array",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "csc_array"))
        expect_equal(zarr_node_format(store, "csc_array"), zarr_format)
        expect_true(zarr_path_exists(store, "csc_array/data"))
        expect_true(zarr_path_exists(store, "csc_array/indices"))
        expect_true(zarr_path_exists(store, "csc_array/indptr"))
        attrs <- Rarr::read_zarr_attributes(file.path(store, "csc_array"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "csc_matrix")

        csr_array <- as(array, "RsparseMatrix")
        expect_silent(write_zarr_element(
          csr_array,
          store,
          "csr_array",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "csr_array"))
        expect_equal(zarr_node_format(store, "csr_array"), zarr_format)
        expect_true(zarr_path_exists(store, "csr_array/data"))
        expect_true(zarr_path_exists(store, "csr_array/indices"))
        expect_true(zarr_path_exists(store, "csr_array/indptr"))
        attrs <- Rarr::read_zarr_attributes(file.path(store, "csr_array"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "csr_matrix")
      })

      test_that("Writing dgeMatrix", {
        value <- matrix(rnorm(20), nrow = 5, ncol = 4) |>
          as("dMatrix") |>
          as("generalMatrix") |>
          as("unpackedMatrix")

        expect_silent(
          write_zarr_element(
            value,
            store,
            "dgematrix",
            zarr_format = zarr_format
          )
        )
        expect_true(zarr_path_exists(store, "dgematrix"))
        expect_equal(zarr_node_format(store, "dgematrix"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "dgematrix"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "array")
      })

      test_that("Writing Zarr nullable booleans works", {
        nullable <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
        nullable[5] <- NA

        expect_silent(write_zarr_element(
          nullable,
          store,
          "nullable_bool",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "nullable_bool"))
        expect_equal(zarr_node_format(store, "nullable_bool"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "nullable_bool"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "nullable-boolean")
      })

      test_that("Writing Zarr nullable integers works", {
        nullable <- as.integer(1:5)
        nullable[5] <- NA

        expect_silent(write_zarr_element(
          nullable,
          store,
          "nullable_int",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "nullable_int"))
        expect_equal(zarr_node_format(store, "nullable_int"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "nullable_int"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "nullable-integer")
      })

      test_that("Writing Zarr integer matrices with NA works", {
        value <- matrix(c(1L, NA, 3:6), nrow = 3, ncol = 2)

        expect_silent(write_zarr_element(
          value,
          store,
          "int_matrix_na",
          zarr_format = zarr_format
        ))
        attrs <- Rarr::read_zarr_attributes(file.path(store, "int_matrix_na"))
        expect_equal(attrs[["encoding-type"]], "array")
        expect_equal(
          read_zarr_element(store, "int_matrix_na"),
          matrix(c(1, NaN, 3:6), nrow = 3, ncol = 2)
        )
      })

      test_that("Writing Zarr string arrays works", {
        string <- LETTERS[1:5]

        write_zarr_element(
          string,
          store,
          "string_array",
          zarr_format = zarr_format
        )
        expect_true(zarr_path_exists(store, "string_array"))
        expect_equal(zarr_node_format(store, "string_array"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "string_array"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "string-array")

        string2d <- matrix(LETTERS[1:20], nrow = 5, ncol = 4)

        expect_silent(write_zarr_element(
          string2d,
          store,
          "string_array2D",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "string_array2D"))
        expect_equal(zarr_node_format(store, "string_array2D"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "string_array2D"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "string-array")
      })

      test_that("Writing Zarr categoricals works", {
        categorical <- factor(LETTERS[1:5])

        expect_no_error(write_zarr_element(
          categorical,
          store,
          "categorical",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "categorical"))
        expect_equal(zarr_node_format(store, "categorical"), zarr_format)
        expect_true(zarr_path_exists(store, "categorical/categories"))
        expect_true(zarr_path_exists(store, "categorical/codes"))
        attrs <- Rarr::read_zarr_attributes(file.path(store, "categorical"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "categorical")
      })

      test_that("Writing Zarr string scalars works", {
        string <- "A"

        expect_silent(write_zarr_element(
          string,
          store,
          "string_scalar",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "string_scalar"))
        expect_equal(zarr_node_format(store, "string_scalar"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "string_scalar"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "string")
      })

      test_that("Writing Zarr numeric scalars works", {
        number <- 1.0

        expect_silent(write_zarr_element(
          number,
          store,
          "numeric_scalar",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "numeric_scalar"))
        expect_equal(zarr_node_format(store, "numeric_scalar"), zarr_format)
        attrs <- Rarr::read_zarr_attributes(file.path(store, "numeric_scalar"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "numeric-scalar")
      })

      test_that("Writing Zarr mappings works", {
        mapping <- list(
          array = matrix(rnorm(20), nrow = 5, ncol = 4),
          sparse = as(matrix(rnorm(20), nrow = 5, ncol = 4), "CsparseMatrix"),
          string = LETTERS[1:5],
          numeric = rnorm(5),
          scalar = 2
        )

        expect_silent(write_zarr_element(
          mapping,
          store,
          "mapping",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "mapping"))
        expect_equal(zarr_node_format(store, "mapping"), zarr_format)
        expect_true(zarr_path_exists(store, "mapping/array"))
        expect_true(zarr_path_exists(store, "mapping/sparse"))
        expect_true(zarr_path_exists(store, "mapping/sparse/data"))
        expect_true(zarr_path_exists(store, "mapping/sparse/indices"))
        expect_true(zarr_path_exists(store, "mapping/sparse/indptr"))
        expect_true(zarr_path_exists(store, "mapping/string"))
        expect_true(zarr_path_exists(store, "mapping/numeric"))
        expect_true(zarr_path_exists(store, "mapping/scalar"))
        attrs <- Rarr::read_zarr_attributes(file.path(store, "mapping"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "dict")
      })

      test_that("Writing Zarr data frames works", {
        df <- data.frame(
          Letters = letters[1:5],
          Numbers = 1:5
        )

        expect_silent(write_zarr_element(
          df,
          store,
          "dataframe",
          zarr_format = zarr_format
        ))
        expect_true(zarr_path_exists(store, "dataframe"))
        expect_equal(zarr_node_format(store, "dataframe"), zarr_format)
        expect_true(zarr_path_exists(store, "dataframe/Letters"))
        expect_true(zarr_path_exists(store, "dataframe/Numbers"))
        expect_true(zarr_path_exists(store, "dataframe/_index"))
        attrs <- Rarr::read_zarr_attributes(file.path(store, "dataframe"))
        expect_true(all(
          c("encoding-type", "encoding-version") %in% names(attrs)
        ))
        expect_equal(attrs[["encoding-type"]], "dataframe")
        expect_true(all(c("_index", "column-order") %in% names(attrs)))
        expect_equal(attrs[["_index"]], "_index")
        expect_identical(
          as.character(attrs[["column-order"]]),
          c("Letters", "Numbers")
        )
      })

      test_that("writing Zarr from SingleCellExperiment works", {
        skip_if_not_installed("SingleCellExperiment")
        store <- tempfile(fileext = ".zarr")
        sce <- generate_dataset(format = "SingleCellExperiment")
        write_zarr(sce, store)
        expect_true(dir.exists(store))
        expect_zarr_store_format(store, zarr_format)
      })

      test_that("writing Zarr from SingleCellExperiment works for manual zarr format", {
        skip_if_not_installed("SingleCellExperiment")
        store <- tempfile(fileext = ".zarr")
        sce <- generate_dataset(format = "SingleCellExperiment")
        zarr_format_manual <- setdiff(c(2, 3), zarr_format)
        write_zarr(sce, store, zarr_format = zarr_format_manual)
        expect_true(dir.exists(store))
        expect_zarr_store_format(store, zarr_format_manual)
        expect_equal(getOption("anndataR.zarr_format"), zarr_format)
      })

      test_that("writing Zarr from Seurat works", {
        skip_if_not_installed("SeuratObject")
        store <- tempfile(fileext = ".zarr")
        seurat <- generate_dataset(format = "Seurat")
        write_zarr(seurat, store)
        expect_true(dir.exists(store))
        expect_zarr_store_format(store, zarr_format)
      })

      test_that("writing Zarr from Seurat works for manual zarr format", {
        skip_if_not_installed("SeuratObject")
        store <- tempfile(fileext = ".zarr")
        seurat <- generate_dataset(format = "Seurat")
        zarr_format_manual <- setdiff(c(2, 3), zarr_format)
        write_zarr(seurat, store, zarr_format = zarr_format_manual)
        expect_true(dir.exists(store))
        expect_zarr_store_format(store, zarr_format_manual)
        expect_equal(getOption("anndataR.zarr_format"), zarr_format)
      })

      # a small AnnData covering the element types that are written as nested
      # groups, as those are where a mixed format store would show up first
      make_adata <- function() {
        AnnData(
          X = matrix(as.double(1:6), nrow = 2, ncol = 3),
          obs = data.frame(
            row.names = c("c1", "c2"),
            categorical = factor(c("a", "b")),
            string = c("x", "y")
          ),
          var = data.frame(row.names = c("g1", "g2", "g3")),
          layers = list(double = matrix(as.double(7:12), nrow = 2, ncol = 3)),
          uns = list(scalar = "a string", nested = list(number = 1))
        )
      }

      test_that("the zarr_format argument overrides the option everywhere", {
        store <- tempfile(fileext = ".zarr")
        zarr_format_manual <- setdiff(c(2, 3), zarr_format)
        write_zarr(make_adata(), store, zarr_format = zarr_format_manual)
        expect_zarr_store_format(store, zarr_format_manual)
        expect_equal(getOption("anndataR.zarr_format"), zarr_format)
      })

      test_that("writing to an existing store keeps the format of that store", {
        store <- tempfile(fileext = ".zarr")
        zarr_format_manual <- setdiff(c(2, 3), zarr_format)
        write_zarr(make_adata(), store, zarr_format = zarr_format_manual)

        # the option says one thing, the store on disk says another
        adata <- suppressWarnings(ZarrAnnData$new(store, mode = "r+"))
        adata$obsm <- list(pca = matrix(as.double(1:4), nrow = 2, ncol = 2))

        expect_zarr_store_format(store, zarr_format_manual)
      })

      test_that("asking for a different format than an existing store errors", {
        store <- tempfile(fileext = ".zarr")
        zarr_format_manual <- setdiff(c(2, 3), zarr_format)
        write_zarr(make_adata(), store, zarr_format = zarr_format_manual)

        expect_error(
          suppressWarnings(
            ZarrAnnData$new(store, mode = "r+", zarr_format = zarr_format)
          ),
          "is in Zarr v"
        )
      })

      test_that("mode 'w' rewrites an existing store in the other format", {
        store <- tempfile(fileext = ".zarr")
        zarr_format_manual <- setdiff(c(2, 3), zarr_format)
        write_zarr(make_adata(), store, zarr_format = zarr_format_manual)

        write_zarr(make_adata(), store, mode = "w", zarr_format = zarr_format)

        expect_zarr_store_format(store, zarr_format)
        expect_no_error(read_zarr(store))
      })

      test_that("an invalid zarr_format errors", {
        for (invalid in list(1, 4, "3", NA, c(2, 3))) {
          expect_error(
            write_zarr(
              make_adata(),
              tempfile(fileext = ".zarr"),
              zarr_format = invalid
            ),
            "must be either"
          )
        }
      })

      dir_size <- function(path) {
        files <- list.files(path, recursive = TRUE, full.names = TRUE)
        sum(file.info(files)$size, na.rm = TRUE)
      }

      test_that("writing compressed files works for Zarr", {
        dummy <- generate_dataset(100, 200)
        non_random_X <- matrix(5, 100, 200) # nolint

        adata <- AnnData(
          X = non_random_X,
          obs = dummy$obs,
          var = dummy$var
        )

        store_none <- tempfile(fileext = ".zarr")
        store_compressed <- tempfile(fileext = ".zarr")

        write_zarr(
          adata,
          store_none,
          compression = "none",
          zarr_format = zarr_format
        )

        comp_list <- c(
          "gzip",
          "blosc",
          "zstd",
          "lzma",
          "bz2",
          "zlib",
          "lz4"
        )
        for (comp in comp_list) {
          write_zarr(
            adata,
            store_compressed,
            compression = comp,
            zarr_format = zarr_format
          )
          # measure before removing the store, otherwise the comparison is
          # always against a size of zero
          expect_gt(dir_size(store_none), dir_size(store_compressed))
          unlink(store_compressed, recursive = TRUE)
        }
      })
    }
  )
}

# Integer matrices with NA are written as for H5AD: as double, with NaN for NA
# (#542)
slots <- c("X", "layers", "obsm", "varm", "obsp", "varp")

for (slot in slots) {
  test_that(
    paste0(
      "an integer matrix with NA in ",
      slot,
      " reads back from Zarr as from H5AD"
    ),
    {
      skip_if_not_installed("rhdf5")
      types <- list(
        x_type = if (slot == "X") {
          "integer_matrix_with_nas"
        } else {
          "numeric_matrix"
        },
        layer_types = character(0),
        obs_types = character(0),
        var_types = character(0),
        obsm_types = character(0),
        varm_types = character(0),
        obsp_types = character(0),
        varp_types = character(0),
        uns_types = character(0)
      )
      if (slot != "X") {
        types[[paste0(
          sub("s$", "", slot),
          "_types"
        )]] <- "integer_matrix_with_nas"
      }
      adata <- do.call(
        generate_dataset,
        c(list(n_obs = 10L, n_vars = 20L, format = "AnnData"), types)
      )

      zarr <- withr::local_tempfile(fileext = ".zarr")
      h5ad <- withr::local_tempfile(fileext = ".h5ad")
      write_zarr(adata, zarr)
      write_h5ad(adata, h5ad)

      from_zarr <- read_zarr(zarr)
      from_h5ad <- read_h5ad(h5ad)
      expect_equal(from_zarr[[slot]], from_h5ad[[slot]])
    }
  )
}

test_that("Python anndata reads an integer matrix with NA written to Zarr", {
  skip_if_no_anndata_py()
  adata <- generate_dataset(
    n_obs = 10L,
    n_vars = 20L,
    x_type = "integer_matrix_with_nas",
    layer_types = character(0),
    obs_types = character(0),
    var_types = character(0),
    obsm_types = character(0),
    varm_types = character(0),
    obsp_types = character(0),
    varp_types = character(0),
    uns_types = character(0),
    format = "AnnData"
  )
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_zarr(adata, zarr)

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  x <- reticulate::py_to_r(np$asarray(ad$read_zarr(zarr)$X))
  expect_equal(is.na(x), is.na(adata$X), ignore_attr = TRUE)
  expect_equal(x[!is.na(x)], as.double(adata$X[!is.na(adata$X)]))
})

for (zarr_format in c(2, 3)) {
  test_that(
    paste0(
      "write_zarr() writes length 1 data frame columns and indices as arrays (v",
      zarr_format,
      ")"
    ),
    {
      store <- withr::local_tempfile(fileext = ".zarr")
      adata <- AnnData(
        X = matrix(1, 1, 2),
        obs = data.frame(
          ct = "T cell",
          n = 3L,
          x = 1.5,
          flag = TRUE,
          na_int = NA_integer_,
          na_chr = NA_character_,
          row.names = "c1"
        ),
        var = data.frame(row.names = c("g1", "g2"))
      )
      write_zarr(adata, store, mode = "w", zarr_format = zarr_format)

      encoding <- function(key) {
        read_zarr_encoding(store, key)$type
      }
      expect_equal(encoding("obs/_index"), "string-array")
      expect_equal(encoding("obs/ct"), "string-array")
      expect_equal(encoding("obs/n"), "array")
      expect_equal(encoding("obs/x"), "array")
      expect_equal(encoding("obs/flag"), "array")
      expect_equal(encoding("obs/na_int"), "nullable-integer")
      expect_equal(encoding("obs/na_chr"), "nullable-string-array")

      # A single variable
      adata <- AnnData(
        X = matrix(1, 2, 1),
        obs = data.frame(row.names = c("c1", "c2")),
        var = data.frame(row.names = "g")
      )
      write_zarr(adata, store, mode = "w", zarr_format = zarr_format)
      expect_equal(encoding("var/_index"), "string-array")

      expect_equal(read_zarr(store)$var_names, "g")
    }
  )
}

for (zarr_format in c(2, 3)) {
  test_that(
    paste0(
      "Python can read Zarr stores with length 1 var (v",
      zarr_format,
      ")"
    ),
    {
      skip_if_no_anndata_py()
      skip_if_no_zarr()
      ad <- reticulate::import("anndata", convert = FALSE)
      store <- withr::local_tempfile(fileext = ".zarr")

      adata <- AnnData(
        X = matrix(1, 2, 1),
        obs = data.frame(row.names = c("c1", "c2")),
        var = data.frame(row.names = "g")
      )
      write_zarr(adata, store, mode = "w", zarr_format = zarr_format)
      py_adata <- ad$read_zarr(store)
      expect_equal(reticulate::py_to_r(py_adata$n_vars), 1L)
    }
  )
}
