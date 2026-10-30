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
