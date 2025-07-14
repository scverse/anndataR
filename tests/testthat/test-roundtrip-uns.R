skip_if_no_anndata()
skip_if_not_installed("reticulate")

library(reticulate)
testthat::skip_if_not(
  reticulate::py_module_available("dummy_anndata"),
  message = "Python dummy_anndata module not available for testing"
)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

test_names <- c(
  names(da$matrix_generators),
  names(da$vector_generators),
  names(da$scalar_generators)
)

for (name in test_names) {
  # first generate a python h5ad
  adata_py <- da$generate_dataset(
    x_type = NULL,
    obs_types = list(),
    var_types = list(),
    layer_types = list(),
    obsm_types = list(),
    varm_types = list(),
    obsp_types = list(),
    varp_types = list(),
    uns_types = list(name),
    nested_uns_types = list()
  )

  # create a couple of paths
  file_py <- withr::local_file(
    tempfile(paste0("anndata_py_", name), fileext = ".h5ad")
  )
  file_r <- withr::local_file(
    tempfile(paste0("anndata_r_", name), fileext = ".h5ad")
  )

  # write to file
  adata_py$write_h5ad(file_py)
  # Read it back in to get the version as read from disk
  adata_py <- ad$read_h5ad(file_py)

  test_that(paste0("Reading an AnnData with uns '", name, "' works"), {
    msg <- message_if_known(
      backend = "HDF5AnnData",
      slot = c("uns"),
      dtype = name,
      process = "read",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    adata_r <- read_h5ad(file_py, as = "HDF5AnnData")
    uns_keys <- if (name == "none") {
      list()
    } else {
      names(adata_r$uns)
    }
    expect_equal(
      uns_keys,
      bi$list(adata_py$uns$keys())
    )

    # check that the print output is the same
    str_r <- capture.output(print(adata_r))
    str_py <- capture.output(print(adata_py))
    expect_equal(str_r, str_py)
  })

  # maybe this test simply shouldn't be run if there is a known issue with reticulate
  test_that(
    paste0("Comparing an anndata with uns '", name, "' with reticulate works"),
    {
      skip_if(name == "none", message = "No value to test for 'none'")
      msg <- message_if_known(
        backend = "HDF5AnnData",
        slot = c("uns"),
        dtype = name,
        process = c("read", "reticulate"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- read_h5ad(file_py, as = "HDF5AnnData")

      # Avoid error because {reticulate} can't convert Python Categorical objects
      if (
        name %in% c(
          "categorical",
          "categorical_missing_values",
          "categorical_ordered",
          "categorical_ordered_missing_values"
        )
      ) {
        categorical <- adata_py$uns[[name]]
        categories <- reticulate::py_to_r(categorical$categories)
        codes <- reticulate::py_to_r(categorical$codes)
        ordered <- reticulate::py_to_r(categorical$ordered)
        is_na <- codes == -1L
        codes[is_na] <- 0L
        py_value <- factor(categories[codes + 1], levels = categories, ordered = ordered)
        py_value[is_na] <- NA
      } else {
        py_value <- reticulate::py_to_r(adata_py$uns[[name]])
      }

      expect_equal(
        adata_r$uns[[name]],
        py_value
      )
    }
  )

  test_that(paste0("Writing an AnnData with uns '", name, "' works"), {
    skip_if(name == "none", message = "No value to test for 'none'")

    msg <- message_if_known(
      backend = "HDF5AnnData",
      slot = c("uns"),
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
      bi$list(adata_py2$uns$keys()),
      name
    )

    # expect that the objects are the same
    expect_equal_py(
      py_get_item(adata_py2$uns, name),
      py_get_item(adata_py$uns, name)
    )
  })
}
