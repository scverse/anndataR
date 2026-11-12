skip_if_no_anndata_py()
skip_if_no_dummy_anndata()
skip_if_no_py_module("h5py")

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
pd <- reticulate::import("pandas", convert = FALSE)
h5py <- reticulate::import("h5py", convert = FALSE)
bi <- reticulate::import_builtins(convert = FALSE)
# anndata >= 0.12 requires opting in to write nullable strings (pd.arrays.StringArray)
tryCatch(ad$settings$allow_write_nullable_strings <- TRUE, error = function(e) {
  NULL
})

known_issues <- read_known_issues()

test_names <- c(
  names(da$matrix_generators),
  names(da$vector_generators),
  names(da$scalar_generators),
  "nullable_string_array"
)

# Open a file or store with h5py/zarr, as Python read_elem()/write_elem() expect
py_open <- function(path, fmt_config, mode) {
  if (fmt_config$backend == "HDF5AnnData") {
    h5py$File(path, mode)
  } else {
    zarr <- reticulate::import("zarr", convert = FALSE)
    zarr$open_group(path, mode = mode, zarr_format = fmt_config$zarr_format)
  }
}

# Close a handle from py_open(), Zarr groups don't need closing
py_close <- function(handle) {
  if (reticulate::py_has_attr(handle, "close")) {
    handle$close()
  }
}

# Write with Python write_elem()
#
# Python write_elem() doesn't allow "/" in keys for every element type, so
# create the parent group first and write the last part of the name to it
py_write_elem <- function(path, fmt_config, name, value) {
  handle <- py_open(path, fmt_config, "a")
  on.exit(py_close(handle))
  group <- if (dirname(name) == ".") {
    handle
  } else {
    handle$require_group(dirname(name))
  }
  ad$io$write_elem(group, basename(name), value)
}

# Read with Python read_elem()
py_read_elem <- function(path, fmt_config, name) {
  handle <- py_open(path, fmt_config, "r")
  on.exit(py_close(handle))
  ad$io$read_elem(py_get_item(handle, name))
}

# Generate a Python value with the same generators as test-roundtrip-uns.R
py_generate_value <- function(name) {
  if (name == "nullable_string_array") {
    # dummy_anndata has no `nullable_string_array` generator
    return(pd$array(
      reticulate::r_to_py(list(NULL, "a", "b", NULL, "c")),
      dtype = "string"
    ))
  }

  adata_py <- da$generate_dataset(
    x_type = NULL,
    obs_types = list(),
    var_types = list(),
    layer_types = list(),
    obsm_types = list(),
    varm_types = list(),
    obsp_types = list(),
    varp_types = list(),
    uns_types = list(name),
    nested_uns_types = list()
  )
  py_get_item(adata_py$uns, name)
}

for (fmt in roundtrip_formats(ad)) {
  fmt_config <- get_fmt_config(fmt)

  withr::with_options(
    list("anndataR.zarr_format" = fmt_config$zarr_format),
    {
      py_set_zarr_write_format(ad, fmt_config$zarr_format)

      for (name in test_names) {
        elem_name <- paste0("elem/", name)

        file_py <- withr::local_file(
          tempfile(paste0("elem_py_", name), fileext = fmt_config$ext)
        )
        file_r <- withr::local_file(
          tempfile(paste0("elem_r_", name), fileext = fmt_config$ext)
        )

        # Write with Python and read it back to get the value as read from disk
        py_write_elem(file_py, fmt_config, elem_name, py_generate_value(name))
        py_value <- py_read_elem(file_py, fmt_config, elem_name)

        test_that(
          paste0(
            "read_elem() reads '",
            name,
            "' written by Python (",
            fmt,
            ")"
          ),
          {
            msg <- message_if_known(
              backend = fmt_config$backend,
              slot = "uns",
              dtype = name,
              process = c("read", "reticulate"),
              known_issues = known_issues
            )
            skip_if(!is.null(msg), message = msg)

            expect_equal(
              read_elem(file_py, elem_name),
              convert_py_value(py_value, name)
            )
          }
        )

        test_that(
          paste0(
            "Python reads '",
            name,
            "' written by write_elem() (",
            fmt,
            ")"
          ),
          {
            msg <- message_if_known(
              backend = fmt_config$backend,
              slot = "uns",
              dtype = name,
              process = c("read", "write"),
              known_issues = known_issues
            )
            skip_if(!is.null(msg), message = msg)

            write_elem(file_r, elem_name, read_elem(file_py, elem_name))

            expect_equal_py(
              py_read_elem(file_r, fmt_config, elem_name),
              py_value
            )
          }
        )

        gc()
      }

      test_that(
        paste0(
          "read_elem() reads AnnData groups written by Python (",
          fmt,
          ")"
        ),
        {
          file_py <- withr::local_file(
            tempfile("elem_py_adata", fileext = fmt_config$ext)
          )
          adata_py <- da$generate_dataset(
            n_obs = 5L,
            n_vars = 4L,
            x_type = "float_matrix",
            obs_types = list("dense_array"),
            var_types = list("dense_array"),
            layer_types = list(),
            obsm_types = list(),
            varm_types = list(),
            obsp_types = list(),
            varp_types = list(),
            uns_types = list(),
            nested_uns_types = list()
          )

          py_write_elem(file_py, fmt_config, "mod/rna", adata_py)
          adata_py <- py_read_elem(file_py, fmt_config, "mod/rna")

          adata_r <- read_elem(file_py, "mod/rna")
          expect_s3_class(adata_r, fmt_config$backend)
          expect_anndata_print_equal(adata_r, adata_py)
          expect_equal(adata_r$X, py_to_r(adata_py$X), ignore_attr = TRUE)
          expect_equal(adata_r$obs_names, py_to_r(adata_py$obs_names$to_list()))
        }
      )

      test_that(
        paste0(
          "Python reads AnnData groups written by write_elem() (",
          fmt,
          ")"
        ),
        {
          file_r <- withr::local_file(
            tempfile("elem_r_adata", fileext = fmt_config$ext)
          )
          adata_r <- AnnData(
            X = matrix(as.numeric(1:20), nrow = 5),
            obs = data.frame(value = 1:5, row.names = paste0("cell", 1:5)),
            var = data.frame(value = 1:4, row.names = paste0("gene", 1:4)),
            layers = list(counts = matrix(as.numeric(20:1), nrow = 5))
          )

          write_elem(file_r, "mod/rna", adata_r)

          adata_py <- py_read_elem(file_r, fmt_config, "mod/rna")
          expect_true(py_to_r(bi$isinstance(adata_py, ad$AnnData)))
          expect_anndata_print_equal(adata_r, adata_py)
          expect_equal(py_to_r(adata_py$X), adata_r$X, ignore_attr = TRUE)
          expect_equal(py_to_r(adata_py$obs_names$to_list()), adata_r$obs_names)
        }
      )
    }
  )
}
