skip_if_no_anndata()
skip_if_no_dummy_anndata()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
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
  # first generate a python h5ad
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
    tempfile(paste0("anndata_py_", name), fileext = ".h5ad")
  )
  file_r <- withr::local_file(
    tempfile(paste0("anndata_r_", name), fileext = ".h5ad")
  )
  file_r2 <- withr::local_file(
    tempfile(paste0("anndata_r2_", name), fileext = ".h5ad")
  )

  # write to file
  adata_py$write_h5ad(file_py)
  # Read it back in to get the version as read from disk
  adata_py <- ad$read_h5ad(file_py)

  test_that(
    paste0("Reading an AnnData with obsm and varm '", name, "' works"),
    {
      msg <- message_if_known(
        backend = "HDF5AnnData",
        slot = c("obsm", "varm"),
        dtype = name,
        process = "read",
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_h5ad(file_py, as = "HDF5AnnData")
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

      # check that the print output is the same
      str_r <- capture.output(print(adata_r))
      str_py <- capture.output(print(adata_py))
      expect_equal(str_r, str_py)
    }
  )

  # maybe this test simply shouldn't be run if there is a known issue with reticulate
  test_that(
    paste0(
      "Comparing an anndata with obsm and varm '",
      name,
      "' with reticulate works"
    ),
    {
      msg <- message_if_known(
        backend = "HDF5AnnData",
        slot = c("obsm", "varm"),
        dtype = name,
        process = c("read", "reticulate"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_h5ad(file_py, as = "HDF5AnnData")

      expect_equal(
        adata_r$obsm[[name]],
        py_to_r(py_get_item(adata_py$obsm, name)),
        tolerance = 1e-6
      )
      expect_equal(
        adata_r$varm[[name]],
        py_to_r(py_get_item(adata_py$varm, name)),
        tolerance = 1e-6
      )
    }
  )

  gc()

  test_that(
    paste0("Writing an AnnData with obsm and varm '", name, "' works"),
    {
      msg <- message_if_known(
        backend = "HDF5AnnData",
        slot = c("obsm", "varm"),
        dtype = name,
        process = c("read", "write"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_h5ad(file_py, as = "InMemoryAnnData")
      write_h5ad(adata_r, file_r)

      # read from file
      adata_py2 <- ad$read_h5ad(file_r)

      # expect name is one of the keys
      expect_contains(
        bi$list(adata_py2$obsm$keys()),
        name
      )
      expect_contains(
        bi$list(adata_py2$obsm$keys()),
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

  skip_if_no_h5diff()
  # Get all R datatypes that are equivalent to the python datatype (name)
  res <- Filter(function(x) x[[1]] == name, all_equivalences)
  r_datatypes <- vapply(res, function(x) x[[2]], character(1))

  for (r_name in r_datatypes) {
    test_msg <- paste0(
      "Comparing a python generated .h5ad with obsm and varm '",
      name,
      "' with an R generated .h5ad '",
      r_name,
      "' works"
    )
    test_that(test_msg, {
      msg <- message_if_known(
        backend = "HDF5AnnData",
        slot = c("obsm", "varm"),
        dtype = c(name, r_name),
        process = c("h5diff"),
        known_issues = known_issues
      )

      skip_if(!is.null(msg), message = msg)
      # generate an R h5ad
      adata_r <- r_generate_dataset(
        10L,
        20L,
        obsm_types = list(r_name),
        varm_types = list(r_name)
      )
      # TODO: Remove this once issue #286 is fixed https://github.com/scverse/anndataR/issues/286
      if (!(r_name %in% names(vector_generators))) {
        adata_r$varm[[r_name]] <- generate_matrix(20, 20, r_name)
      }
      write_h5ad(adata_r, file_r2, mode = "w")

      # Remove the rhdf5-NA.OK for comparison
      hdf5_clear_rhdf5_attributes(file_r2, paste0("/obsm/", r_name))

      # run h5diff
      res_obsm <- processx::run(
        "h5diff",
        c(
          "-v2",
          file_py,
          file_r2,
          paste0("/obsm/", name),
          paste0("/obsm/", r_name)
        ),
        error_on_status = FALSE
      )
      expect_equal(res_obsm$status, 0, info = res_obsm$stdout)

      # Remove the rhdf5-NA.OK for comparison
      hdf5_clear_rhdf5_attributes(file_r2, paste0("/varm/", r_name))

      res_varm <- processx::run(
        "h5diff",
        c(
          "-v2",
          file_py,
          file_r2,
          paste0("/varm/", name),
          paste0("/varm/", r_name)
        ),
        error_on_status = FALSE
      )
      expect_equal(res_varm$status, 0, info = res_varm$stdout)
    })
  }
}
