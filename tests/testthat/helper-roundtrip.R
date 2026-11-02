#' Get format config for roundtrip tests
#'
#' Get a list of backend-specific values for a given file format
#'
#' @param fmt Either `"h5ad"` or `"zarr"`
#' @return A named list with elements: `backend`, `ext`, `r_read_fun`,
#'   `r_write_fun`, `py_read_method`, `py_write_method`
get_fmt_config <- function(fmt = c("h5ad", "zarrv2", "zarrv3")) {
  fmt <- match.arg(fmt)

  if (fmt %in% c("zarrv2", "zarrv3")) {
    skip_if_no_zarr() # nolint: object_usage_linter
    list(
      backend = "ZarrAnnData",
      ext = ".zarr",
      r_read_fun = read_zarr,
      r_write_fun = write_zarr,
      py_read_method = "read_zarr",
      py_write_method = "write_zarr",
      zarr_format = if (fmt == "zarrv2") 2L else 3L
    )
  } else {
    list(
      backend = "HDF5AnnData",
      ext = ".h5ad",
      r_read_fun = read_h5ad,
      r_write_fun = write_h5ad,
      py_read_method = "read_h5ad",
      py_write_method = "write_h5ad"
    )
  }
}

#' Normalise Python AnnData print output
#'
#' Since anndata 0.13.2 the Python repr lists `X` among the layers as
#' `None (.X)`. anndataR does not, so this entry is removed from the `layers:`
#' line, and the line is dropped when no other layer is left.
#'
#' @param lines Character vector of printed lines
#' @return The lines without the `X` layer entry
normalise_py_anndata_print <- function(lines) {
  # remove the x entry and its separator from the layers line
  is_layers <- grepl("^ *layers: ", lines)
  lines[is_layers] <- sub(
    ", None \\(\\.X\\)|None \\(\\.X\\)(, )?",
    "",
    lines[is_layers]
  )
  lines[!grepl("^ *layers: *$", lines)]
}

#' Formats for the roundtrip tests
#'
#' All formats Python anndata can write. anndata main dropped
#' `anndata.settings.zarr_write_format` and only writes Zarr v3
#' (scverse/anndata#2628); Zarr v2 is then recorded as a skipped test and left
#' out, so the loop over the formats doesn't need to know.
#'
#' @param ad The Python `anndata` module
#' @return A character vector of formats accepted by `get_fmt_config()`
roundtrip_formats <- function(ad) {
  fmts <- c("h5ad", "zarrv2", "zarrv3")
  if (reticulate::py_has_attr(ad$settings, "zarr_write_format")) {
    return(fmts)
  }
  testthat::test_that("Python anndata can write Zarr v2", {
    testthat::skip("Python anndata cannot write Zarr v2")
  })
  setdiff(fmts, "zarrv2")
}

#' Set the Zarr format Python anndata writes
#'
#' Does nothing for H5AD (`zarr_format = NULL`) or when the setting no longer
#' exists, in which case `roundtrip_formats()` has already left out Zarr v2.
#'
#' @param ad The Python `anndata` module
#' @param zarr_format The Zarr format to write, `2L`, `3L` or `NULL`
py_set_zarr_write_format <- function(ad, zarr_format) {
  if (
    !is.null(zarr_format) &&
      reticulate::py_has_attr(ad$settings, "zarr_write_format")
  ) {
    ad$settings$zarr_write_format <- zarr_format
  }
  invisible(zarr_format)
}

#' Expect AnnData print output to match
#'
#' Compares the print output of an R AnnData object with a Python AnnData
#' object, normalising backend class names and the Python `X` layer entry
#' before comparing.
#'
#' @param adata_r An R AnnData object
#' @param adata_py A Python AnnData object
expect_anndata_print_equal <- function(adata_r, adata_py) {
  str_r <- capture.output(print(adata_r))
  str_py <- normalise_py_anndata_print(capture.output(print(adata_py)))

  # Normalise class names in R output to match Python output
  str_r <- gsub("[^ ]*AnnData", "AnnData", str_r)

  testthat::expect_equal(str_r, str_py)
}
