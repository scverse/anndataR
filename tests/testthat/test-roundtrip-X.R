skip_if_no_anndata_py()
skip_if_no_dummy_anndata()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

test_names <- names(da$matrix_generators)

# X must always be 2-dimensional in AnnData
# -> https://github.com/scverse/anndata/blob/2a2c0e3198c298a5c80a73ac343c63203b5ca133/src/anndata/_core/anndata.py#L2164-L2172 # nolint
test_names <- test_names[!grepl("_3d$", test_names)]

for (fmt in c("h5ad", "zarr")) {
  fmt_config <- get_fmt_config(fmt)

  for (name in test_names) {
    # first generate a python adata
    adata_py <- da$generate_dataset(
      x_type = name,
      obs_types = list(),
      var_types = list(),
      layer_types = list(),
      obsm_types = list(),
      varm_types = list(),
      obsp_types = list(),
      varp_types = list(),
      uns_types = list(),
      nested_uns_types = list()
    )

    # create a couple of paths
    file_py <- withr::local_file(
      tempfile(paste0("anndata_py_", name), fileext = fmt_config$ext)
    )
    file_r <- withr::local_file(
      tempfile(paste0("anndata_r_", name), fileext = fmt_config$ext)
    )
    file_r2 <- withr::local_file(
      tempfile(paste0("anndata_r2_", name), fileext = fmt_config$ext)
    )

    # write to file
    adata_py[[fmt_config$py_write_method]](file_py)
    # Read it back in to get the version as read from disk
    adata_py <- ad[[fmt_config$py_read_method]](file_py)

    test_that(
      paste0("Reading an AnnData with X '", name, "' (", fmt, ") works"),
      {
        msg <- message_if_known(
          backend = fmt_config$backend,
          slot = c("X"),
          dtype = name,
          process = "read",
          known_issues = known_issues
        )
        skip_if(!is.null(msg), message = msg)

        adata_r <- fmt_config$r_read_fun(file_py, as = fmt_config$backend)
        expect_equal(
          adata_r$shape(),
          unlist(reticulate::py_to_r(adata_py$shape))
        )

        # check that the print output is the same (normalize class names)
        expect_anndata_print_equal(adata_r, adata_py)
      }
    )

    test_that(
      paste0(
        "Comparing an anndata with X '",
        name,
        "' (",
        fmt,
        ") with reticulate works"
      ),
      {
        msg <- message_if_known(
          backend = fmt_config$backend,
          slot = c("X"),
          dtype = name,
          process = c("read", "reticulate"),
          known_issues = known_issues
        )
        skip_if(!is.null(msg), message = msg)

        adata_r <- fmt_config$r_read_fun(file_py, as = fmt_config$backend)

        # Extract X matrices, removing dimnames for comparison since
        # R AnnData adds dimnames on-the-fly but Python doesn't preserve them
        actual_x <- adata_r$X
        expected_x <- py_to_r(adata_py$X)
        dimnames(actual_x) <- NULL
        dimnames(expected_x) <- NULL

        expect_equal(
          actual_x,
          expected_x,
          tolerance = 1e-6
        )
      }
    )

    gc()

    test_that(
      paste0("Writing an AnnData with X '", name, "' (", fmt, ") works"),
      {
        msg <- message_if_known(
          backend = fmt_config$backend,
          slot = c("X"),
          dtype = name,
          process = c("read", "write"),
          known_issues = known_issues
        )
        skip_if(!is.null(msg), message = msg)

        adata_r <- fmt_config$r_read_fun(file_py, as = "InMemoryAnnData")
        fmt_config$r_write_fun(adata_r, file_r)

        # read from file
        adata_py2 <- ad[[fmt_config$py_read_method]](file_r)

        # expect that the objects are the same
        expect_equal_py(
          adata_py2$X,
          adata_py$X
        )
      }
    )

    if (fmt == "h5ad") {
      skip_if_no_h5diff()
      # Get all R datatypes that are equivalent to the python datatype (name)
      res <- Filter(function(x) x[[1]] == name, matrix_equivalences)
      r_datatypes <- vapply(res, function(x) x[[2]], character(1))

      for (r_name in r_datatypes) {
        test_msg <- paste0(
          "Comparing a python generated .h5ad with X '",
          name,
          "' with an R generated .h5ad '",
          r_name,
          "' works"
        )
        test_that(test_msg, {
          msg <- message_if_known(
            backend = "HDF5AnnData",
            slot = c("X"),
            dtype = c(name, r_name),
            process = c("h5diff"),
            known_issues = known_issues
          )
          skip_if(!is.null(msg), message = msg)

          # generate an R h5ad
          adata_r <- r_generate_dataset(10L, 20L, x_type = list(r_name))
          write_h5ad(adata_r, file_r2, mode = "w")

          # Remove the rhdf5-NA.OK for comparison
          hdf5_clear_rhdf5_attributes(file_r2, "X")

          # run h5diff
          res <- processx::run(
            "h5diff",
            c("-v2", file_py, file_r2, "/X"),
            error_on_status = FALSE
          )

          expect_equal(res$status, 0, info = res$stdout)
        })
      }
    }
  }
}
