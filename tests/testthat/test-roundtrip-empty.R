skip_if_no_anndata()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

# first generate a python h5ad
adata_py <- ad$AnnData()

name <- "empty"

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

test_that(paste0("Reading an AnnData with layer '", name, "' works"), {
  msg <- message_if_known(
    backend = "HDF5AnnData",
    slot = c("none"),
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

  # check that the print output is the same
  str_r <- capture.output(print(adata_r))
  str_py <- capture.output(print(adata_py))
  expect_equal(str_r, str_py)
})

gc()

test_that(paste0("Writing an AnnData with layer '", name, "' works"), {
  msg <- message_if_known(
    backend = "HDF5AnnData",
    slot = c("none"),
    dtype = name,
    process = c("read", "write"),
    known_issues = known_issues
  )
  skip_if(!is.null(msg), message = msg)

  adata_r <- read_h5ad(file_py, as = "InMemoryAnnData")
  write_h5ad(adata_r, file_r)

  # read from file
  adata_py2 <- ad$read_h5ad(file_r)

  # check that the print output is the same
  expect_equal(
    unlist(reticulate::py_to_r(adata_py2$shape)),
    unlist(reticulate::py_to_r(adata_py$shape))
  )

  # check that the print output is the same
  str_py2 <- capture.output(print(adata_py2))
  str_py <- capture.output(print(adata_py))
  expect_equal(str_py2, str_py)
})
