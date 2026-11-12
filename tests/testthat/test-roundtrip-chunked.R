skip_if_no_anndata_py()
skip_if_no_dummy_anndata()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
bi <- reticulate::import_builtins()

py_keys <- function(mapping) py_to_r(bi$list(mapping$keys()))
chunk_dim <- function(store, name) {
  overview <- Rarr::zarr_overview(file.path(store, name), as_data_frame = TRUE)
  unlist(overview$chunk_dim)
}
slots <- c("X", "obs", "var", "layers", "obsm", "varm", "obsp", "varp", "uns")
mapping_slots <- c("layers", "obsm", "varm", "obsp", "varp")

# A tiny chunk_size splits every array into many chunks. Reading such a store
# must give exactly what a store written as single chunks (chunk_size = NULL)
# gives, in R and in Python.
for (fmt in c("zarrv2", "zarrv3")) {
  fmt_config <- get_fmt_config(fmt)

  withr::with_options(list("anndataR.zarr_format" = fmt_config$zarr_format), {
    ad$settings$zarr_write_format <- fmt_config$zarr_format

    file_py <- withr::local_file(tempfile("anndata_py_", fileext = ".zarr"))
    file_single <- withr::local_file(tempfile("single_", fileext = ".zarr"))
    file_chunked <- withr::local_file(tempfile("chunked_", fileext = ".zarr"))

    # every element type of dummy_anndata
    da$generate_dataset()$write_zarr(file_py)
    adata_r <- read_zarr(file_py, as = "InMemoryAnnData")
    suppressWarnings({
      write_zarr(adata_r, file_single, chunk_size = NULL)
      write_zarr(adata_r, file_chunked, chunk_size = 4L)
    })

    test_that(paste0("chunk_size splits arrays into chunks (", fmt, ")"), {
      for (name in c("obs/_index", "var/_index")) {
        expect_lt(chunk_dim(file_chunked, name), chunk_dim(file_single, name))
      }
    })

    test_that(
      paste0("A chunked store reads as a single-chunk one in R (", fmt, ")"),
      {
        single <- read_zarr(file_single, as = "InMemoryAnnData")
        chunked <- read_zarr(file_chunked, as = "InMemoryAnnData")
        for (slot in slots) {
          expect_equal(chunked[[slot]], single[[slot]], info = slot)
        }
      }
    )

    test_that(
      paste0(
        "A chunked store reads as a single-chunk one in Python (",
        fmt,
        ")"
      ),
      {
        single <- ad$read_zarr(file_single)
        chunked <- ad$read_zarr(file_chunked)
        expect_equal_py(chunked$X, single$X)
        expect_equal_py(chunked$obs, single$obs)
        expect_equal_py(chunked$var, single$var)
        for (slot in mapping_slots) {
          # the order of the keys follows the store's directory listing
          keys <- sort(py_keys(single[[slot]]))
          expect_equal(sort(py_keys(chunked[[slot]])), keys, info = slot)
          for (key in keys) {
            expect_equal_py(chunked[[slot]]$get(key), single[[slot]]$get(key))
          }
        }
      }
    )
  })
}
