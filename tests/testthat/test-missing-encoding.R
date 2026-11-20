known_issues <- read_known_issues()

uns <- list(
  a = 1L,
  nested = list(citation = "some citation", kept = "kept")
)
obs <- data.frame(
  row.names = c("cell1", "cell2"),
  missing = c(1.5, 2.5),
  kept = c("x", "y")
)
adata <- AnnData(
  obs = obs,
  var = data.frame(row.names = "gene"),
  uns = uns
)
missing <- c("uns/nested/citation", "obs/missing")

# Write `adata` and remove the encoding metadata of the `missing` elements, like
# `uns/citation` in some CELLxGENE files (#583)
write_without_encoding <- function(format) {
  if (format == "h5ad") {
    path <- withr::local_tempfile(
      fileext = ".h5ad",
      .local_envir = parent.frame()
    )
    adata$write_h5ad(path)
    for (element in missing) {
      for (attr in c("encoding-type", "encoding-version")) {
        rhdf5::h5deleteAttribute(path, element, attr)
      }
    }
    return(path)
  }

  zarr_format <- if (format == "zarr_v2") 2L else 3L
  path <- withr::local_tempfile(
    fileext = ".zarr",
    .local_envir = parent.frame()
  )
  adata$write_zarr(path, zarr_format = zarr_format)
  for (element in file.path(path, missing)) {
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
      "an element without encoding metadata is skipped with a warning,",
      format
    ),
    {
      skip_if_not_installed(fmt$pkg)
      path <- write_without_encoding(format)

      warnings <- capture_warnings(
        result <- fmt$read(path, as = "InMemoryAnnData")
      )
      expect_length(warnings, 2)
      for (element in missing) {
        expect_match(
          warnings,
          paste0("Skipping element ", element),
          all = FALSE
        )
      }

      # Its parent and siblings are kept
      expect_identical(result$uns$a, uns$a)
      expect_identical(result$uns$nested, uns$nested["kept"])
      expect_identical(result$obs, obs["kept"])
    }
  )

  test_that(
    paste(
      "read_elem() skips nested elements without encoding metadata,",
      format
    ),
    {
      skip_if_not_installed(fmt$pkg)
      path <- write_without_encoding(format)

      expect_warning(
        result <- read_elem(path, "uns"),
        "Skipping element uns/nested/citation"
      )
      expect_identical(result$nested, uns$nested["kept"])

      # Reading the element itself is still an error, naming the element
      expect_error(
        read_elem(path, "uns/nested/citation"),
        "\"uns/nested/citation\".*has no.*encoding-type"
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
      expect_identical(result$uns, uns)
      expect_identical(result$obs, obs)
    }
  )
}
