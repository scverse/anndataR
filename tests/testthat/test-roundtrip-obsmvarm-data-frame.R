skip_if_no_anndata_py()

library(reticulate)

ad <- reticulate::import("anndata", convert = FALSE)
pd <- reticulate::import("pandas", convert = FALSE)

# Python anndata requires the index of a data frame in obsm/varm to be equal to
# the obs_names/var_names of the parent object and refuses to read the whole
# file if it is not. These tests guard against writing such a file.

n_obs <- 10L
n_vars <- 20L
obs_names <- paste0("Cell", seq_len(n_obs))
var_names <- paste0("Gene", seq_len(n_vars))

r_adata <- function(obsm, varm) {
  AnnData(
    X = matrix(runif(n_obs * n_vars), n_obs, n_vars),
    obs = data.frame(row.names = obs_names),
    var = data.frame(row.names = var_names),
    obsm = obsm,
    varm = varm
  )
}

for (fmt in c("h5ad", "zarrv2", "zarrv3")) {
  fmt_config <- get_fmt_config(fmt)
  # The Zarr backend needs {Rarr}, which `skip_if_no_zarr()` does not check
  skip_without <- if (fmt == "h5ad") "rhdf5" else "Rarr"

  withr::with_options(
    list("anndataR.zarr_format" = fmt_config$zarr_format),
    {
      if (grepl("zarr", fmt_config$ext, fixed = TRUE)) {
        ad$settings$zarr_write_format <- fmt_config$zarr_format
      }

      test_that(
        paste0(
          "Writing an AnnData with obsm and varm data frames (",
          fmt,
          ") can be read by Python"
        ),
        {
          skip_if_not_installed(skip_without)

          file_r <- withr::local_file(
            tempfile("anndata_r_obsmvarm_df", fileext = fmt_config$ext)
          )

          adata_r <- r_adata(
            obsm = list(
              df = data.frame(
                A = runif(n_obs),
                B = runif(n_obs),
                row.names = obs_names
              )
            ),
            varm = list(
              df = data.frame(C = runif(n_vars), row.names = var_names)
            )
          )
          fmt_config$r_write_fun(adata_r, file_r)

          adata_py <- ad[[fmt_config$py_read_method]](file_r)

          # The index must match the parent names, not the automatic 1:nrow
          # sequence that R adds when row names are stripped
          expect_equal(
            as.character(py_to_r(
              py_get_item(adata_py$obsm, "df")$index$values
            )),
            obs_names
          )
          expect_equal(
            as.character(py_to_r(
              py_get_item(adata_py$varm, "df")$index$values
            )),
            var_names
          )

          expect_equal(
            as.character(py_to_r(
              py_get_item(adata_py$obsm, "df")$columns$values
            )),
            c("A", "B")
          )
          expect_equal(
            py_to_r(py_get_item(adata_py$obsm, "df")$to_numpy()),
            as.matrix(adata_r$obsm$df),
            ignore_attr = TRUE,
            tolerance = 1e-6
          )
        }
      )

      test_that(
        paste0(
          "Writing an AnnData with obsm data frames without row names (",
          fmt,
          ") can be read by Python"
        ),
        {
          skip_if_not_installed(skip_without)

          file_r <- withr::local_file(
            tempfile("anndata_r_obsmvarm_nodf", fileext = fmt_config$ext)
          )

          # Row names are optional, the index is taken from the parent object
          adata_r <- r_adata(
            obsm = list(df = data.frame(A = runif(n_obs))),
            varm = list(df = data.frame(C = runif(n_vars)))
          )
          fmt_config$r_write_fun(adata_r, file_r)

          adata_py <- ad[[fmt_config$py_read_method]](file_r)
          expect_equal(
            as.character(py_to_r(
              py_get_item(adata_py$obsm, "df")$index$values
            )),
            obs_names
          )
          expect_equal(
            as.character(py_to_r(
              py_get_item(adata_py$varm, "df")$index$values
            )),
            var_names
          )
        }
      )

      test_that(
        paste0(
          "Reading and rewriting a Python AnnData with obsm and varm data ",
          "frames (",
          fmt,
          ") works"
        ),
        {
          skip_if_not_installed(skip_without)

          file_py <- withr::local_file(
            tempfile("anndata_py_obsmvarm_df", fileext = fmt_config$ext)
          )
          file_r <- withr::local_file(
            tempfile("anndata_r2_obsmvarm_df", fileext = fmt_config$ext)
          )

          adata_py <- ad$AnnData(
            X = matrix(runif(n_obs * n_vars), n_obs, n_vars),
            obs = pd$DataFrame(index = obs_names),
            var = pd$DataFrame(index = var_names)
          )
          py_set_item(
            adata_py$obsm,
            "df",
            pd$DataFrame(
              list(A = runif(n_obs)),
              index = obs_names
            )
          )
          py_set_item(
            adata_py$varm,
            "df",
            pd$DataFrame(
              list(C = runif(n_vars)),
              index = var_names
            )
          )
          adata_py[[fmt_config$py_write_method]](file_py)

          adata_r <- fmt_config$r_read_fun(file_py, as = "InMemoryAnnData")
          expect_identical(rownames(adata_r$obsm$df), obs_names)
          expect_identical(rownames(adata_r$varm$df), var_names)

          fmt_config$r_write_fun(adata_r, file_r)

          adata_py2 <- ad[[fmt_config$py_read_method]](file_r)
          expect_equal_py(
            py_get_item(adata_py2$obsm, "df"),
            py_get_item(adata_py$obsm, "df")
          )
          expect_equal_py(
            py_get_item(adata_py2$varm, "df"),
            py_get_item(adata_py$varm, "df")
          )
        }
      )
    }
  )
}
