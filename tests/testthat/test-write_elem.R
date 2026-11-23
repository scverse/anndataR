skip_if_not_installed("rhdf5")

adata <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")

# Every top-level slot and one level below
elem_names <- c(
  "X",
  "obs",
  "var",
  "layers",
  "obsm",
  "varm",
  "obsp",
  "varp",
  "uns",
  unlist(lapply(c("layers", "obsm", "uns"), function(slot) {
    paste0(slot, "/", names(adata[[slot]]))
  }))
)

ref_h5ad <- withr::local_tempfile(fileext = ".h5ad")
write_h5ad(adata, ref_h5ad)

# Write each element to its own group so nested names don't clash with the
# slot they are in
elem_target <- function(name) {
  paste0("elem/", gsub("/", "__", name, fixed = TRUE))
}

# TODO: Remove once scalar NA strings are written to Zarr correctly (#581)
drop_scalar_na_string <- function(value) {
  if (is.list(value) && !is.data.frame(value)) {
    value$scalar_character_with_nas <- NULL
  }
  value
}

test_that("write_elem() roundtrips every element for H5AD", {
  path <- withr::local_tempfile(fileext = ".h5ad")

  for (name in elem_names) {
    value <- read_elem(ref_h5ad, name)
    write_elem(path, elem_target(name), value)
    expect_equal(read_elem(path, elem_target(name)), value, info = name)
  }
})

for (zarr_format in c(2L, 3L)) {
  test_that(
    paste0("write_elem() roundtrips every element for Zarr v", zarr_format),
    {
      skip_if_not_installed("Rarr")
      withr::local_options(anndataR.zarr_format = zarr_format)

      path <- withr::local_tempfile(fileext = ".zarr")

      for (name in setdiff(elem_names, "uns/scalar_character_with_nas")) {
        value <- drop_scalar_na_string(read_elem(ref_h5ad, name))
        write_elem(path, elem_target(name), value)
        expect_equal(read_elem(path, elem_target(name)), value, info = name)
      }
      expect_identical(get_zarr_format(path), zarr_format)
    }
  )
}

test_that("write_elem() keeps other elements of an existing H5AD file", {
  path <- withr::local_tempfile(fileext = ".h5ad")
  write_h5ad(adata, path)
  before <- read_h5ad(path)

  write_elem(path, "uns/new", 1:3)

  after <- read_h5ad(path)
  expect_equal(as.vector(after$uns$new), 1:3)
  expect_equal(after$X, before$X)
  expect_equal(after$obs, before$obs)
  expect_equal(after$layers, before$layers)
  expect_equal(after$uns[names(before$uns)], before$uns)
})

test_that("write_elem() keeps other elements of an existing Zarr store", {
  skip_if_not_installed("Rarr")

  for (zarr_format in c(2L, 3L)) {
    info <- paste0("Zarr v", zarr_format)
    path <- withr::local_tempfile(fileext = ".zarr")
    write_zarr(adata, path, zarr_format = zarr_format)
    before <- read_zarr(path)

    # The store format is used even if the option is different
    withr::with_options(list(anndataR.zarr_format = 5L - zarr_format), {
      write_elem(path, "uns/new", 1:3)
    })

    expect_identical(
      zarr_node_format(path, "uns/new"),
      zarr_format,
      info = info
    )
    expect_equal(as.vector(read_elem(path, "uns/new")), 1:3, info = info)
    after <- read_zarr(path)
    expect_equal(after$X, before$X, info = info)
    expect_equal(after$obs, before$obs, info = info)
    expect_equal(after$layers, before$layers, info = info)
    expect_equal(after$uns[names(before$uns)], before$uns, info = info)
  }
})

test_that("write_elem() creates parent groups", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_elem(h5ad, "/a/b/c", 1:3)
  expect_equal(as.vector(read_elem(h5ad, "a/b/c")), 1:3)
  expect_equal(as.vector(read_elem(h5ad, "a")$b$c), 1:3)

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_elem(zarr, "/a/b/c", 1:3)
  expect_equal(as.vector(read_elem(zarr, "a/b/c")), 1:3)
  expect_equal(as.vector(read_elem(zarr, "a")$b$c), 1:3)
})

test_that("write_elem() parent groups can be read back with the AnnData", {
  for (format in c("h5ad", "zarr")) {
    if (format == "zarr") {
      skip_if_not_installed("Rarr")
    }
    path <- withr::local_tempfile(fileext = paste0(".", format))
    small <- AnnData(
      obs = data.frame(row.names = "cell"),
      var = data.frame(row.names = "gene"),
      uns = list(a = 1L)
    )
    if (format == "h5ad") {
      write_h5ad(small, path)
    } else {
      write_zarr(small, path)
    }

    write_elem(path, "uns/results/scores", c(0.1, 0.2))

    read <- if (format == "h5ad") read_h5ad else read_zarr
    expect_no_warning(result <- read(path, as = "InMemoryAnnData"))
    expect_equal(result$uns$a, 1L, info = format)
    expect_equal(
      as.vector(result$uns$results$scores),
      c(0.1, 0.2),
      info = format
    )
  }
})

for (format in c("h5ad", "zarr_v2", "zarr_v3")) {
  test_that(
    paste("write_elem() checks where an element would be written,", format),
    {
      small <- generate_dataset(n_obs = 3L, n_vars = 2L, format = "AnnData")
      if (format == "h5ad") {
        path <- withr::local_tempfile(fileext = ".h5ad")
        write_h5ad(small, path)
      } else {
        skip_if_not_installed("Rarr")
        path <- withr::local_tempfile(fileext = ".zarr")
        zarr_format <- if (format == "zarr_v2") 2L else 3L
        write_zarr(small, path, zarr_format = zarr_format)
      }

      # Nothing can be written inside an array
      expect_error(
        write_elem(path, "X/foo", 1:3),
        "\"X\" is an array, not a group"
      )
      expect_error(
        write_elem(path, "X/foo", small),
        "\"X\" is an array, not a group"
      )

      # An AnnData can't be written inside an encoded element
      expect_error(
        write_elem(path, "uns/sub", small),
        "\"uns\" is an encoded \"dict\" element"
      )

      # An AnnData can only replace an AnnData or an empty group
      expect_error(
        write_elem(path, "uns", small, overwrite = TRUE),
        "\"uns\" is not an AnnData"
      )

      # Nothing was written
      read <- if (format == "h5ad") read_h5ad else read_zarr
      result <- read(path, as = "InMemoryAnnData")
      expect_equal(result$X, small$X)
      expect_setequal(names(result$uns), names(small$uns))

      # A plain group can hold an AnnData, as in a .h5mu file
      write_elem(path, "mod/rna", small)
      expect_equal(read_elem(path, "mod/rna")$X, small$X)
    }
  )
}

test_that("write_elem() only replaces existing elements with overwrite = TRUE", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  write_elem(h5ad, "x", 1:3)
  expect_error(write_elem(h5ad, "x", 4:6), "already exists")
  expect_equal(as.vector(read_elem(h5ad, "x")), 1:3)
  write_elem(h5ad, "x", letters, overwrite = TRUE)
  expect_equal(as.vector(read_elem(h5ad, "x")), letters)

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  write_elem(zarr, "x", 1:3)
  expect_error(write_elem(zarr, "x", 4:6), "already exists")
  expect_equal(as.vector(read_elem(zarr, "x")), 1:3)
  write_elem(zarr, "x", letters, overwrite = TRUE)
  expect_equal(as.vector(read_elem(zarr, "x")), letters)
})

test_that("write_elem() follows the anndataR.write_null option", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  withr::with_options(list(anndataR.write_null = TRUE), {
    write_elem(h5ad, "null", NULL)
  })
  expect_null(read_elem(h5ad, "null"))

  withr::with_options(list(anndataR.write_null = FALSE), {
    write_elem(h5ad, "not_written", NULL)
  })
  expect_error(read_elem(h5ad, "not_written"), "not found")
})

test_that("write_elem() passes options in ... to the writers", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  mat <- matrix(as.numeric(1:200), nrow = 10)

  write_elem(h5ad, "gzip", mat, compression = "gzip", chunk_size = NULL)
  expect_equal(read_elem(h5ad, "gzip"), mat)
  expect_error(
    write_elem(h5ad, "bad", mat, compression = "bad"),
    "should be one of"
  )

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  zarr_chunk_dim <- function(name) {
    unlist(
      Rarr::zarr_overview(
        file.path(zarr, name),
        as_data_frame = TRUE
      )$chunk_dim
    )
  }

  write_elem(zarr, "gzip", mat, compression = "gzip")
  expect_equal(read_elem(zarr, "gzip"), mat)

  write_elem(zarr, "single", mat, chunk_size = NULL)
  expect_equal(zarr_chunk_dim("single"), dim(mat))

  # 400 bytes is 50 doubles, so the 200 values are split
  write_elem(zarr, "chunked", mat, chunk_size = 400)
  expect_true(prod(zarr_chunk_dim("chunked")) < length(mat))
  expect_equal(read_elem(zarr, "chunked"), mat)
})

test_that("write_elem() rejects an invalid chunk_size before writing to Zarr", {
  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  mat <- matrix(as.numeric(1:6), nrow = 2)
  write_elem(zarr, "x", mat)

  expect_error(
    write_elem(zarr, "x", mat * 2, overwrite = TRUE, chunk_size = -1),
    "chunk_size"
  )
  expect_equal(read_elem(zarr, "x"), mat)
})

test_that("write_elem() errors for unsupported arguments in ...", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  expect_error(write_elem(h5ad, "x", 1:3, compresion = "gzip"), "compresion")
  expect_error(
    write_elem(h5ad, "x", 1:3, FALSE, "auto", "gzip"),
    "must be named"
  )
  expect_error(write_elem(h5ad, "x", 1:3, zarr_format = 2), "zarr_format")
  expect_error(
    write_elem(h5ad, "x", 1:3, stop_on_error = FALSE),
    "stop_on_error"
  )

  skip_if_not_installed("Rarr")
  zarr <- withr::local_tempfile(fileext = ".zarr")
  expect_error(write_elem(zarr, "x", 1:3, zarr_format = 2), "zarr_format")
})

test_that("write_elem() errors for non-AnnData roots and unsupported values", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  expect_error(write_elem(h5ad, "/", 1:3), "root group")
  expect_error(write_elem(h5ad, "", 1:3), "root group")
  expect_error(write_elem(h5ad, "x", function() NULL), "not supported")
})

test_that("write_elem() validates its arguments", {
  dir <- withr::local_tempdir()
  h5ad <- file.path(dir, "test.h5ad")

  expect_error(write_elem(file.path(dir, "test.txt"), "x", 1:3), "format")
  expect_error(write_elem(h5ad, NULL, 1:3), "single string")
  expect_error(write_elem(h5ad, "x", 1:3, overwrite = NA), "overwrite")
  expect_error(write_elem(h5ad, "x", 1:3, format = "loom"), "must be one of")
  expect_false(file.exists(h5ad))
})

test_that("write_elem() creates a Zarr store in an empty directory", {
  skip_if_not_installed("Rarr")
  dir <- withr::local_tempdir()

  empty_zarr <- file.path(dir, "empty.zarr")
  dir.create(empty_zarr)
  write_elem(empty_zarr, "x", 1:3)
  expect_equal(as.vector(read_elem(empty_zarr, "x")), 1:3)

  empty <- file.path(dir, "empty")
  dir.create(empty)
  write_elem(empty, "x", 1:3, format = "zarr")
  expect_equal(as.vector(read_elem(empty, "x")), 1:3)
})

test_that("a failed write_elem() leaves no new file or store behind", {
  dir <- withr::local_tempdir()

  for (ext in c("h5ad", "zarr")) {
    if (ext == "zarr") {
      skip_if_not_installed("Rarr")
    }
    path <- file.path(dir, paste0("new.", ext))

    expect_error(write_elem(path, "/", 1), "root group")
    expect_false(file.exists(path))

    expect_error(write_elem(path, "x", new.env()))
    expect_false(file.exists(path))

    # An existing file or store is kept
    write_elem(path, "x", 1:3)
    expect_error(write_elem(path, "y", new.env()))
    expect_equal(as.vector(read_elem(path, "x")), 1:3)
  }
})

test_that("write_elem() returns the path invisibly", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")

  expect_invisible(write_elem(h5ad, "x", 1:3))
  expect_identical(write_elem(h5ad, "y", 1:3), h5ad)
})

# Mappings are read back sorted by name and dense Matrix objects are read back
# as base matrices, so normalise both before comparing
normalise_slot <- function(x) {
  if (!is.list(x) || is.data.frame(x)) {
    return(x)
  }
  lapply(x[sort(names(x))], function(value) {
    if (inherits(value, "denseMatrix")) as.matrix(value) else value
  })
}

expect_adata_equal <- function(actual, expected, info = NULL) {
  for (slot in c("X", "obs", "var", "layers", "obsm", "obsp")) {
    testthat::expect_equal(
      normalise_slot(actual[[slot]]),
      normalise_slot(expected[[slot]]),
      info = paste(info, slot)
    )
  }
}

small_adata <- generate_dataset(
  n_obs = 3L,
  n_vars = 4L,
  layer_types = character(0),
  obsm_types = character(0),
  varm_types = character(0),
  obsp_types = character(0),
  varp_types = character(0),
  uns_types = character(0),
  format = "AnnData"
)

for (fmt in c("h5ad", "zarrv2", "zarrv3")) {
  ext <- if (fmt == "h5ad") ".h5ad" else ".zarr"
  zarr_format <- if (fmt == "zarrv2") 2L else 3L

  test_that(paste0("write_elem() writes AnnData objects to groups for ", fmt), {
    if (fmt != "h5ad") {
      skip_if_not_installed("Rarr")
    }
    withr::local_options(anndataR.zarr_format = zarr_format)
    path <- withr::local_tempfile(fileext = ext)

    write_elem(path, "other", 1:3)
    write_elem(path, "mod/rna", adata, compression = "gzip")

    result <- read_elem(path, "mod/rna")
    expect_s3_class(result, if (fmt == "h5ad") "HDF5AnnData" else "ZarrAnnData")
    expect_adata_equal(result, adata, info = fmt)
    expect_equal(as.vector(read_elem(path, "other")), 1:3)

    # Existing groups are only replaced with overwrite = TRUE
    expect_error(write_elem(path, "mod/rna", small_adata), "already exists")
    write_elem(path, "mod/rna", small_adata, overwrite = TRUE)
    expect_adata_equal(read_elem(path, "mod/rna"), small_adata, info = fmt)
    expect_equal(as.vector(read_elem(path, "other")), 1:3)
  })

  test_that(
    paste0("write_elem() writes AnnData objects to the root for ", fmt),
    {
      if (fmt != "h5ad") {
        skip_if_not_installed("Rarr")
      }
      withr::local_options(anndataR.zarr_format = zarr_format)
      path <- withr::local_tempfile(fileext = ext)

      # A new path is the same as write_h5ad()/write_zarr()
      write_elem(path, "/", adata)
      expect_adata_equal(read_elem(path, "/"), adata, info = fmt)
      result <- if (fmt == "h5ad") read_h5ad(path) else read_zarr(path)
      expect_adata_equal(result, adata, info = fmt)

      # A non-empty file is only replaced with overwrite = TRUE
      write_elem(path, "other", 1:3)
      expect_error(write_elem(path, "/", small_adata), "not empty")
      write_elem(path, "/", small_adata, overwrite = TRUE)
      expect_adata_equal(read_elem(path, "/"), small_adata, info = fmt)
      expect_error(read_elem(path, "other"), "not found")
    }
  )
}
