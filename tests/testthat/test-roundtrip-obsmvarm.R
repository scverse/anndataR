skip_if_no_anndata_py()
skip_if_no_dummy_anndata()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
da <- reticulate::import("dummy_anndata", convert = FALSE)
bi <- reticulate::import_builtins()

known_issues <- read_known_issues()

test_names <- c(
  names(da$matrix_generators),
  names(da$vector_generators)
)

# temporary workaround for
# https://github.com/data-intuitive/dummy-anndata/issues/12
test_names <- setdiff(
  test_names,
  c(
    "categorical",
    "categorical_missing_values",
    "categorical_ordered",
    "categorical_ordered_missing_values",
    "nullable_boolean_array",
    "nullable_integer_array"
  )
)

# Types to test, each paired with the obsm/varm key it is stored under.
# dummy_anndata collapses every requested `df_<type>` into a single element
# named "dataframe", so for data frames the key is not the name of the type.
# The dummy-anndata#12 types excluded above are included as data frames
# because they round trip correctly as columns of one, as they do for obs/var.
test_elements <- c(
  lapply(test_names, function(type) list(type = type, key = type)),
  lapply(names(da$vector_generators), function(type) {
    list(type = paste0("df_", type), key = "dataframe")
  })
)

for (fmt in c("h5ad", "zarrv2", "zarrv3")) {
  fmt_config <- get_fmt_config(fmt)

  withr::with_options(
    list("anndataR.zarr_format" = fmt_config$zarr_format),
    {
      if (grepl("zarr", fmt_config$ext, fixed = TRUE)) {
        ad$settings$zarr_write_format <- fmt_config$zarr_format
      }

      for (element in test_elements) {
        name <- element$type
        key <- element$key

        # first generate a python adata
        adata_py <- da$generate_dataset(
          x_type = NULL,
          obs_types = list(),
          var_types = list(),
          layer_types = list(),
          obsm_types = list(name),
          varm_types = list(name),
          obsp_types = list(),
          varp_types = list(),
          uns_types = list(),
          nested_uns_types = list()
        )

        # create a couple of paths
        file_py <- withr::local_file(
          tempfile(paste0("anndata_py_", name), fileext = fmt_config$ext)
        )
        file_r <- withr::local_file(
          tempfile(paste0("anndata_r_", name), fileext = fmt_config$ext)
        )
        file_r2 <- withr::local_file(
          tempfile(paste0("anndata_r2_", name), fileext = fmt_config$ext)
        )

        # write to file
        adata_py[[fmt_config$py_write_method]](file_py)
        # Read it back in to get the version as read from disk
        adata_py <- ad[[fmt_config$py_read_method]](file_py)

        test_that(
          paste0(
            "Reading an AnnData with obsm and varm '",
            name,
            "' (",
            fmt,
            ") works"
          ),
          {
            msg <- message_if_known(
              backend = fmt_config$backend,
              slot = c("obsm", "varm"),
              dtype = name,
              process = "read",
              known_issues = known_issues
            )
            skip_if(!is.null(msg), message = msg)

            adata_r <- fmt_config$r_read_fun(file_py, as = fmt_config$backend)
            expect_equal(
              adata_r$shape(),
              unlist(reticulate::py_to_r(adata_py$shape))
            )
            expect_equal(
              adata_r$obsm_keys(),
              bi$list(adata_py$obsm$keys())
            )
            expect_equal(
              adata_r$varm_keys(),
              bi$list(adata_py$varm$keys())
            )

            # check that the print output is the same (normalize class names)
            expect_anndata_print_equal(adata_r, adata_py)
          }
        )

        test_that(
          paste0(
            "Comparing an anndata with obsm and varm '",
            name,
            "' (",
            fmt,
            ") with reticulate works"
          ),
          {
            msg <- message_if_known(
              backend = fmt_config$backend,
              slot = c("obsm", "varm"),
              dtype = name,
              process = c("read", "reticulate"),
              known_issues = known_issues
            )
            skip_if(!is.null(msg), message = msg)

            adata_r <- fmt_config$r_read_fun(file_py, as = fmt_config$backend)

            # R AnnData now adds dimnames on-the-fly, but Python doesn't preserve them
            # So we need to strip dimnames for comparison. A data.frame always has
            # row names, and `py_to_r()` keeps the pandas index as the row.names
            # attribute, so the data frame is rebuilt to make them comparable
            strip_names <- function(x) {
              if (is.data.frame(x)) {
                x <- as.data.frame(as.list(x), check.names = FALSE)
                rownames(x) <- NULL
              } else {
                dimnames(x) <- NULL
              }
              x
            }

            actual_obsm <- strip_names(adata_r$obsm[[key]])
            expected_obsm <- strip_names(
              py_to_r(py_get_item(adata_py$obsm, key))
            )

            expect_equal(
              actual_obsm,
              expected_obsm,
              tolerance = 1e-6
            )

            actual_varm <- strip_names(adata_r$varm[[key]])
            expected_varm <- strip_names(
              py_to_r(py_get_item(adata_py$varm, key))
            )

            expect_equal(
              actual_varm,
              expected_varm,
              tolerance = 1e-6
            )
          }
        )

        gc()

        test_that(
          paste0(
            "Writing an AnnData with obsm and varm '",
            name,
            "' (",
            fmt,
            ") works"
          ),
          {
            msg <- message_if_known(
              backend = fmt_config$backend,
              slot = c("obsm", "varm"),
              dtype = name,
              process = c("read", "write"),
              known_issues = known_issues
            )
            skip_if(!is.null(msg), message = msg)

            adata_r <- fmt_config$r_read_fun(file_py, as = "InMemoryAnnData")
            fmt_config$r_write_fun(adata_r, file_r)

            # read from file
            adata_py2 <- ad[[fmt_config$py_read_method]](file_r)

            # expect the element is one of the keys
            expect_contains(
              bi$list(adata_py2$obsm$keys()),
              key
            )
            expect_contains(
              bi$list(adata_py2$varm$keys()),
              key
            )

            # expect that the objects are the same. Python `anndata` requires
            # the index of a data frame in obsm/varm to be the obs_names or
            # var_names of the parent object and refuses to read the file at
            # all if it is not, so for data frames this also checks the index
            expect_equal_py(
              py_get_item(adata_py2$obsm, key),
              py_get_item(adata_py$obsm, key)
            )
            expect_equal_py(
              py_get_item(adata_py2$varm, key),
              py_get_item(adata_py$varm, key)
            )
          }
        )

        if (fmt == "h5ad") {
          skip_if_no_h5diff()
          # Get all R datatypes that are equivalent to the python datatype (name)
          res <- Filter(function(x) x[[1]] == name, all_equivalences)
          r_datatypes <- vapply(res, function(x) x[[2]], character(1))

          for (r_name in r_datatypes) {
            test_msg <- paste0(
              "Comparing a python generated .h5ad with obsm and varm '",
              name,
              "' with an R generated .h5ad '",
              r_name,
              "' works"
            )
            test_that(test_msg, {
              msg <- message_if_known(
                backend = "HDF5AnnData",
                slot = c("obsm", "varm"),
                dtype = c(name, r_name),
                process = c("h5diff"),
                known_issues = known_issues
              )

              skip_if(!is.null(msg), message = msg)
              # generate an R h5ad
              adata_r <- r_generate_dataset(
                10L,
                20L,
                obsm_types = list(r_name),
                varm_types = list(r_name)
              )
              write_h5ad(adata_r, file_r2, mode = "w")

              # When the R element is a data frame, its index must be the
              # obs_names/var_names of the parent object. This cannot be
              # checked against {file_py} because the Python and R generators
              # use different observation and variable names, so compare the
              # index against the parent names within {file_r2}. Do this before
              # the attributes are cleared below, so that both datasets still
              # carry the same attributes.
              for (slot in c("obsm", "varm")) {
                encoding <- rhdf5::h5readAttributes(
                  file_r2,
                  paste0("/", slot, "/", r_name)
                )[["encoding-type"]]

                if (!identical(encoding, "dataframe")) {
                  next
                }

                axis <- if (slot == "obsm") "obs" else "var"
                res_index <- processx::run(
                  "h5diff",
                  c(
                    "-v2",
                    file_r2,
                    file_r2,
                    paste0("/", axis, "/_index"),
                    paste0("/", slot, "/", r_name, "/_index")
                  ),
                  error_on_status = FALSE
                )
                expect_equal(res_index$status, 0, info = res_index$stdout)
              }

              hdf5_file_r2 <- HDF5File$new(file_r2)

              # Remove the rhdf5-NA.OK for comparison
              hdf5_clear_rhdf5_attributes(
                hdf5_file_r2,
                paste0("/obsm/", r_name)
              )

              # run h5diff
              res_obsm <- processx::run(
                "h5diff",
                c(
                  "-v2",
                  file_py,
                  file_r2,
                  paste0("/obsm/", name),
                  paste0("/obsm/", r_name)
                ),
                error_on_status = FALSE
              )
              expect_equal(res_obsm$status, 0, info = res_obsm$stdout)

              # Remove the rhdf5-NA.OK for comparison
              hdf5_clear_rhdf5_attributes(
                hdf5_file_r2,
                paste0("/varm/", r_name)
              )

              res_varm <- processx::run(
                "h5diff",
                c(
                  "-v2",
                  file_py,
                  file_r2,
                  paste0("/varm/", name),
                  paste0("/varm/", r_name)
                ),
                error_on_status = FALSE
              )
              expect_equal(res_varm$status, 0, info = res_varm$stdout)
            })
          }
        }
      }
    }
  )
}
