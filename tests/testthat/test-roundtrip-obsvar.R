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

test_names <- names(da$vector_generators)

for (name in test_names) {
  # first generate a python h5ad
  adata_py <- da$generate_dataset(
    x_type = NULL,
    obs_types = list(name),
    var_types = list(name),
    layer_types = list(),
    obsm_types = list(),
    varm_types = list(),
    obsp_types = list(),
    varp_types = list(),
    uns_types = list(),
    nested_uns_types = list()
  )

  # create a couple of paths
  file_py <- withr::local_file(tempfile(paste0("anndata_py_", name), fileext = ".h5ad"))
  file_r <- withr::local_file(tempfile(paste0("anndata_r_", name), fileext = ".h5ad"))
  file_r2 <- withr::local_file(tempfile(paste0("anndata_r2_", name), fileext = ".h5ad"))

  # write to file
  adata_py$write_h5ad(file_py)

  test_that(paste0("reading an AnnData with obs and var '", name, "' works"), {
    msg <- message_if_known(
      backend = "HDF5AnnData",
      slot = c("obs", "var"),
      dtype = name,
      process = "read",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    adata_r <- read_h5ad(file_py, to = "HDF5AnnData")
    expect_equal(
      adata_r$shape(),
      unlist(reticulate::py_to_r(adata_py$shape))
    )
    expect_equal(
      adata_r$obs_keys(),
      bi$list(adata_py$obs_keys())
    )
    expect_equal(
      adata_r$var_keys(),
      bi$list(adata_py$var_keys())
    )

    # check that the print output is the same
    str_r <- capture.output(print(adata_r))
    str_py <- capture.output(print(adata_py))
    expect_equal(str_r, str_py)
  })

  # maybe this test simply shouldn't be run if there is a known issue with reticulate
  test_that(paste0("Comparing an anndata with obs and var '", name, "' with reticulate works"), {
    msg <- message_if_known(
      backend = "HDF5AnnData",
      slot = c("obs", "var"),
      dtype = name,
      process = c("read", "reticulate"),
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    adata_r <- read_h5ad(file_py, to = "HDF5AnnData")

    expect_equal(
      adata_r$obs[[name]],
      py_to_r(adata_py$obs)[[name]],
      tolerance = 1e-6
    )
    expect_equal(
      adata_r$var[[name]],
      py_to_r(adata_py$var)[[name]],
      tolerance = 1e-6
    )
  })

  test_that(paste0("Writing an AnnData with obs and var '", name, "' works"), {
    msg <- message_if_known(
      backend = "HDF5AnnData",
      slot = c("obsp", "varp"),
      dtype = name,
      process = c("read", "write"),
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    adata_r <- read_h5ad(file_py, to = "InMemoryAnnData")
    write_h5ad(adata_r, file_r)

    # read from file
    adata_py2 <- ad$read_h5ad(file_r)

    # expect name is one of the keys
    expect_contains(
      bi$list(adata_py2$obs$keys()),
      name
    )
    expect_contains(
      bi$list(adata_py2$var$keys()),
      name
    )

    # expect that the objects are the same
    expect_equal_py(adata_py2$obs, adata_py$obs)
    expect_equal_py(adata_py2$var, adata_py$var)
  })

  skip_if_no_h5diff()
  # Get all R datatypes that are equivalent to the python datatype (name)
  res <- Filter(function(x) x[[1]] == name, vector_equivalences)
  r_datatypes <- sapply(res, function(x) x[[2]])

  for (r_name in r_datatypes){
    test_msg <- paste0("Comparing a python generated .h5ad with obs and var '", name,
                       "' with an R generated .h5ad '", r_name, "' works")
    test_that(test_msg, {
      msg <- message_if_known(
        backend = "HDF5AnnData",
        slot = c("obs", "var"),
        dtype = c(name, r_name),
        process = c("h5diff"),
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)
      # generate an R h5ad
      adata_r <- r_generate_dataset(10L, 20L, obs_types = list(r_name), var_types = list(r_name))
      write_h5ad(adata_r, file_r2)

      # run h5diff
      res_obs <- processx::run("h5diff",
                               c("-v", file_py, file_r2, paste0("/obs/", name), paste0("/obs/", r_name)),
                               error_on_status = FALSE)
      expect_equal(res_obs$status, 0, info = res_obs$stdout)

      res_var <- processx::run("h5diff",
                               c("-v", file_py, file_r2, paste0("/var/", name), paste0("/var/", r_name)),
                               error_on_status = FALSE)
      expect_equal(res_var$status, 0, info = res_var$stdout)
    })
  }
}
