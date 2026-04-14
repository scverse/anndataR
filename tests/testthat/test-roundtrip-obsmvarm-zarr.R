skip_if_no_anndata_py()
skip_if_no_dummy_anndata()
skip_if_no_zarr()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
zr <- reticulate::import("zarr", convert = FALSE)
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

test_names <- c(
  names(da$matrix_generators),
  names(da$vector_generators)
)

# temporary workaround for
# https://github.com/data-intuitive/dummy-anndata/issues/12
test_names <- setdiff(
  test_names,
  c(
    "categorical",
    "categorical_missing_values",
    "categorical_ordered",
    "categorical_ordered_missing_values",
    "nullable_boolean_array",
    "nullable_integer_array"
  )
)

for (name in test_names) {
  # first generate a python zarr
  adata_py <- da$generate_dataset(
    x_type = NULL,
    obs_types = list(),
    var_types = list(),
    layer_types = list(),
    obsm_types = list(name),
    varm_types = list(name),
    obsp_types = list(),
    varp_types = list(),
    uns_types = list(),
    nested_uns_types = list()
  )

  # create a couple of paths
  file_py <- withr::local_file(
    tempfile(paste0("anndata_py_", name), fileext = ".zarr")
  )
  file_r <- withr::local_file(
    tempfile(paste0("anndata_r_", name), fileext = ".zarr")
  )
  file_r2 <- withr::local_file(
    tempfile(paste0("anndata_r2_", name), fileext = ".zarr")
  )

  # write to file
  adata_py$write_zarr(file_py)
  # Read it back in to get the version as read from disk
  adata_py <- ad$read_zarr(file_py)

  test_that(
    paste0("Reading an AnnData with obsm and varm '", name, "' works"),
    {
      msg <- message_if_known(
        backend = "ZarrAnnData",
        slot = c("obsm", "varm"),
        dtype = name,
        process = "read",
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_zarr(file_py, as = "ZarrAnnData")
      expect_equal(
        adata_r$shape(),
        unlist(reticulate::py_to_r(adata_py$shape))
      )
      expect_equal(
        adata_r$obsm_keys(),
        bi$list(adata_py$obsm$keys())
      )
      expect_equal(
        adata_r$varm_keys(),
        bi$list(adata_py$varm$keys())
      )

      # check that the print output is the same (normalize class names)
      str_r <- capture.output(print(adata_r))
      str_py <- capture.output(print(adata_py))
      str_r <- gsub("[^ ]*AnnData", "AnnData", str_r)
      expect_equal(str_r, str_py)
    }
  )

  test_that(
    paste0(
      "Comparing an anndata with obsm and varm '",
      name,
      "' with reticulate works"
    ),
    {
      msg <- message_if_known(
        backend = "ZarrAnnData",
        slot = c("obsm", "varm"),
        dtype = name,
        process = c("read", "reticulate"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_zarr(file_py, as = "ZarrAnnData")

      # R AnnData now adds dimnames on-the-fly, but Python doesn't preserve them
      # So we need to strip dimnames for comparison
      actual_obsm <- adata_r$obsm[[name]]
      expected_obsm <- py_to_r(py_get_item(adata_py$obsm, name))
      dimnames(actual_obsm) <- NULL
      dimnames(expected_obsm) <- NULL

      expect_equal(
        actual_obsm,
        expected_obsm,
        tolerance = 1e-6
      )

      actual_varm <- adata_r$varm[[name]]
      expected_varm <- py_to_r(py_get_item(adata_py$varm, name))
      dimnames(actual_varm) <- NULL
      dimnames(expected_varm) <- NULL

      expect_equal(
        actual_varm,
        expected_varm,
        tolerance = 1e-6
      )
    }
  )

  gc()

  test_that(
    paste0("Writing an AnnData with obsm and varm '", name, "' works"),
    {
      msg <- message_if_known(
        backend = "ZarrAnnData",
        slot = c("obsm", "varm"),
        dtype = name,
        process = c("read", "write"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_zarr(file_py, as = "InMemoryAnnData")
      write_zarr(adata_r, file_r)

      # read from file
      adata_py2 <- ad$read_zarr(file_r)

      # expect name is one of the keys
      expect_contains(
        bi$list(adata_py2$obsm$keys()),
        name
      )
      expect_contains(
        bi$list(adata_py2$varm$keys()),
        name
      )

      # expect that the objects are the same
      expect_equal_py(
        py_get_item(adata_py2$obsm, name),
        py_get_item(adata_py$obsm, name)
      )
      expect_equal_py(
        py_get_item(adata_py2$varm, name),
        py_get_item(adata_py$varm, name)
      )
    }
  )

  # TODO: is there a way to compare two zarr stores
}
