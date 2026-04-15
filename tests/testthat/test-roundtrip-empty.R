skip_if_no_anndata_py()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
# anndata >= 0.12 requires opting in to write nullable strings (pd.arrays.StringArray)
# see https://github.com/scverse/anndata/issues/2221
tryCatch(ad$settings$allow_write_nullable_strings <- TRUE, error = function(e) {
  NULL
})
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

name <- "empty"

for (fmt in c("h5ad", "zarr")) {
  fmt_config <- get_fmt_config(fmt)

  # first generate a python adata
  adata_py <- ad$AnnData()

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
    paste0("Reading an AnnData with layer '", name, "' (", fmt, ") works"),
    {
      msg <- message_if_known(
        backend = fmt_config$backend,
        slot = c("none"),
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

  gc()

  test_that(
    paste0("Writing an AnnData with layer '", name, "' (", fmt, ") works"),
    {
      msg <- message_if_known(
        backend = fmt_config$backend,
        slot = c("none"),
        dtype = name,
        process = c("read", "write"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)

      adata_r <- fmt_config$r_read_fun(file_py, as = "InMemoryAnnData")
      fmt_config$r_write_fun(adata_r, file_r)

      # read from file
      adata_py2 <- ad[[fmt_config$py_read_method]](file_r)

      # check that the shape is the same
      expect_equal(
        unlist(reticulate::py_to_r(adata_py2$shape)),
        unlist(reticulate::py_to_r(adata_py$shape))
      )

      # check that the print output is the same
      str_py2 <- capture.output(print(adata_py2))
      str_py <- capture.output(print(adata_py))
      expect_equal(str_py2, str_py)
    }
  )
}
