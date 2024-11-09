skip_if_no_anndata()
skip_if_not_installed("reticulate")

library(reticulate)
testthat::skip_if_not(
  reticulate::py_module_available("dummy_anndata"),
  message = "Python dummy_anndata module not available for testing"
)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
pd <- reticulate::import("pandas", convert = FALSE)
bi <- reticulate::import_builtins()

test_names <- names(da$vector_generators)

for (name in test_names) {
  test_that(paste0("roundtrip with obs and var '", name, "'"), {
    adata_py <- da$generate_dataset(
      x_type = "generate_float_matrix",
      obs_types = list(name),
      var_types = list(name),
      layer_types = list(),
      obsm_types = list(),
      varm_types = list(),
      obsp_types = list(),
      varp_types = list(),
      uns_types = list()
    )
    # remove uns - workaround for https://github.com/data-intuitive/dummy-anndata/issues/2
    adata_py$uns <- bi$dict()
    # TODO: remove X

    # write to file
    filename <- withr::local_file(tempfile(fileext = ".h5ad"))
    adata_py$write_h5ad(filename)

    # read from file
    adata_r <- read_h5ad(filename, to = "InMemoryAnnData")

    # simple checks first
    expect_equal(
      adata_r$shape(),
      unlist(reticulate::py_to_r(adata_py$shape))
    )
    expect_equal(
      adata_r$obs_keys(),
      py_to_r(adata_py$obs_keys())
    )
    expect_equal(
      adata_r$var_keys(),
      py_to_r(adata_py$var_keys())
    )

    # check that the print output is the same
    str_r <- capture.output(print(adata_r))
    str_py <- capture.output(print(adata_py))
    expect_equal(str_r, str_py)

    # if we would test the objects at this stage,
    # we're also testing reticulate's conversion
    # nolint start
    # expect_equal(
    #   adata_r$obs[[name]],
    #   py_to_r(adata_py$obs[[name]]),
    #   tolerance = 1e-6
    # )
    # expect_equal(
    #   adata_r$var[[name]],
    #   py_to_r(adata_py$var[[name]]),
    #   tolerance = 1e-6
    # )
    # nolint end

    # write back to file
    filename2 <- withr::local_file(tempfile(fileext = ".h5ad"))
    write_h5ad(adata_r, filename2)

    # read from file
    adata_py2 <- ad$read_h5ad(filename2)

    # expect slots are unchanged
    zz <- pd$testing$assert_frame_equal(
      adata_py$obs,
      adata_py2$obs,
      check_dtype = FALSE,
      check_exact = FALSE
    )
    expect_null(reticulate::py_to_r(zz))
    pd$testing$assert_frame_equal(
      adata_py$var,
      adata_py2$var,
      check_dtype = FALSE,
      check_exact = FALSE
    )
    expect_null(reticulate::py_to_r(zz))
  })
}
