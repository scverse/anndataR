skip_if_no_anndata_py()
skip_if_no_dummy_anndata()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

test_names <- c(
  names(da$matrix_generators),
  names(da$vector_generators),
  names(da$scalar_generators)
)

for (fmt in c("h5ad", "zarr")) {
  fmt_config <- get_fmt_config(fmt)

  for (name in test_names) {
    # first generate a python adata
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
      tempfile(paste0("anndata_py_", name), fileext = fmt_config$ext)
    )
    file_r <- withr::local_file(
      tempfile(paste0("anndata_r_", name), fileext = fmt_config$ext)
    )

    # write to file
    adata_py[[fmt_config$py_write_method]](file_py)
    # Read it back in to get the version as read from disk
    adata_py <- ad[[fmt_config$py_read_method]](file_py)

    test_that(
      paste0("Reading an AnnData with uns '", name, "' (", fmt, ") works"),
      {
        msg <- message_if_known(
          backend = fmt_config$backend,
          slot = c("uns"),
          dtype = name,
          process = "read",
          known_issues = known_issues
        )
        skip_if(!is.null(msg), message = msg)

        adata_r <- fmt_config$r_read_fun(file_py, as = fmt_config$backend)

        expect_equal(
          names(adata_r$uns),
          bi$list(adata_py$uns$keys())
        )

        # check that the print output is the same (normalize class names)
        expect_anndata_print_equal(adata_r, adata_py)
      }
    )

    test_that(
      paste0(
        "Comparing an anndata with uns '",
        name,
        "' (",
        fmt,
        ") with reticulate works"
      ),
      {
        msg <- message_if_known(
          backend = fmt_config$backend,
          slot = c("uns"),
          dtype = name,
          process = c("read", "reticulate"),
          known_issues = known_issues
        )
        skip_if(!is.null(msg), message = msg)

        adata_r <- fmt_config$r_read_fun(file_py, as = fmt_config$backend)

        py_value <- convert_py_value(adata_py$uns[[name]], name)

        expect_equal(
          adata_r$uns[[name]],
          py_value
        )
      }
    )

    gc()

    test_that(
      paste0("Writing an AnnData with uns '", name, "' (", fmt, ") works"),
      {
        msg <- message_if_known(
          backend = fmt_config$backend,
          slot = c("uns"),
          dtype = name,
          process = c("read", "write"),
          known_issues = known_issues
        )
        skip_if(!is.null(msg), message = msg)

        adata_r <- fmt_config$r_read_fun(file_py, as = "InMemoryAnnData")
        fmt_config$r_write_fun(adata_r, file_r)

        # read from file
        adata_py2 <- ad[[fmt_config$py_read_method]](file_r)

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
      }
    )
  }
}
