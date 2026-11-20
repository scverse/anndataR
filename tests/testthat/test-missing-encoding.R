known_issues <- read_known_issues()

uns <- list(a = 1L, nested = list(citation = "some citation"))
adata <- AnnData(
  obs = data.frame(row.names = "cell"),
  var = data.frame(row.names = "gene"),
  uns = uns
)
missing <- "uns/nested/citation"

# Write `adata` and remove the encoding metadata of `missing`, like
# `uns/citation` in some CELLxGENE files (#583)
write_without_encoding <- function(format) {
  if (format == "h5ad") {
    path <- withr::local_tempfile(
      fileext = ".h5ad",
      .local_envir = parent.frame()
    )
    adata$write_h5ad(path)
    for (attr in c("encoding-type", "encoding-version")) {
      rhdf5::h5deleteAttribute(path, missing, attr)
    }
    return(path)
  }

  zarr_format <- if (format == "zarr_v2") 2L else 3L
  path <- withr::local_tempfile(
    fileext = ".zarr",
    .local_envir = parent.frame()
  )
  adata$write_zarr(path, zarr_format = zarr_format)
  element <- file.path(path, missing)
  if (zarr_format == 2L) {
    unlink(file.path(element, ".zattrs"))
  } else {
    metadata <- jsonlite::read_json(file.path(element, "zarr.json"))
    metadata$attributes <- NULL
    jsonlite::write_json(
      metadata,
      file.path(element, "zarr.json"),
      auto_unbox = TRUE,
      null = "null"
    )
  }
  path
}

formats <- list(
  h5ad = list(backend = "HDF5AnnData", read = read_h5ad, pkg = "rhdf5"),
  zarr_v2 = list(backend = "ZarrAnnData", read = read_zarr, pkg = "Rarr"),
  zarr_v3 = list(backend = "ZarrAnnData", read = read_zarr, pkg = "Rarr")
)

for (format in names(formats)) {
  fmt <- formats[[format]]

  test_that(
    paste(
      "an element without encoding metadata is named in the error,",
      format
    ),
    {
      skip_if_not_installed(fmt$pkg)
      path <- write_without_encoding(format)

      expect_error(
        read_elem(path, missing),
        paste0("\"", missing, "\".*has no.*encoding-type")
      )
      expect_warning(
        fmt$read(path, as = "InMemoryAnnData"),
        paste0("\"", missing, "\"")
      )
    }
  )

  test_that(
    paste("an element without encoding metadata is read,", format),
    {
      msg <- message_if_known(
        backend = fmt$backend,
        slot = "uns",
        dtype = "no_encoding",
        process = "read",
        known_issues = known_issues
      )
      skip_if(!is.null(msg), message = msg)
      skip_if_not_installed(fmt$pkg)
      path <- write_without_encoding(format)

      # Python anndata reads these elements with an OldFormatWarning
      expect_warning(result <- fmt$read(path, as = "InMemoryAnnData"))
      expect_identical(result$uns$a, uns$a)
      expect_identical(
        as.vector(result$uns$nested$citation),
        uns$nested$citation
      )
    }
  )
}
