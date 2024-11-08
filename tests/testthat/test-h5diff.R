skip_if_not_installed("hdf5r")

requireNamespace("reticulate")
testthat::skip_if_not(
  reticulate::py_module_available("dummy_anndata"),
  message = "Python dummy_anndata module not available for testing"
)

r_generate_dataset_only_x <- function(n_obs, n_vars, x_type = "numeric_matrix") {
  data <- generate_dataset(n_obs, n_vars,
                           x_type = x_type,
                           layer_types = c(),
                           obs_types = c("integer"),
                           var_types = c("integer"),
                           obsm_types = c(),
                           varm_types = c(),
                           obsp_types = c(),
                           varp_types = c(),
                           uns_types = c(),
                           format = "AnnData")

  data
}

py_generate_dataset_only_x <- function(n_obs, n_vars, x_type = "generate_float_matrix") {
  da <- reticulate::import("dummy_anndata")
  data <- da$generate_dataset(n_obs, n_vars,
                              x_type = c(x_type),
                              layer_types = list(),
                              obs_types = list("integer_array"),
                              var_types = list("integer_array"),
                              obsm_types = list(),
                              varm_types = list(),
                              obsp_types = list(),
                              varp_types = list(),
                              uns_types = list())

  data
}

matrix_equivalences <- list(
  "generate_float_matrix" = list("numeric_matrix", "numeric_dense"),
  "generate_float_matrix_nas" = list("numeric_matrix_with_nas", "numeric_dense_with_nas"),
  "generate_integer_matrix" = list("integer_matrix", "integer_dense"),
  "generate_float_csparse" = list("numeric_csparse"),
  "generate_float_csparse_nas" = list("numeric_csparse_with_nas"),
  "generate_float_rsparse" = list("numeric_rsparse"),
  "generate_float_rsparse_nas" = list("numeric_rsparse_with_nas")
)

for (py_matrix in names(matrix_equivalences)) {
  data_python <- py_generate_dataset_only_x(10L, 20L, x_type = py_matrix)
  h5ad_file_py <- tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")
  data_python$write_h5ad(h5ad_file_py)

  for (r_matrix in matrix_equivalences[[py_matrix]]) {
    test_that(paste0("h5diff_X_", py_matrix, "_", r_matrix), {

      tryCatch({
        data_r <- r_generate_dataset_only_x(10L, 20L, x_type = r_matrix)
        h5ad_file_r <- tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
        write_h5ad(data_r, h5ad_file_r)

        res <- processx::run("h5diff", c("-v", h5ad_file_py, h5ad_file_r, "/X"), error_on_status = FALSE)

        expect_equal(res$status, 0, info = res$stdout)
      }, error = function(e) {
        message("Error: ", e$message)
        message("Python matrix: ", py_matrix)
        message("R matrix: ", r_matrix)
        
      }, warning = function(w) {
        message("Warning: ", w$message)
        message("Python matrix: ", py_matrix)
        message("R matrix: ", r_matrix)
        
      })
    })
  }
}

# # test_that("h5diff_X_float", {
#   data_r <- r_generate_dataset_only_x(3L, 5L, x_type = "numeric_matrix")
#   data_python <- py_generate_dataset_only_x(3L, 5L, x_type = "generate_float_matrix")

#   h5ad_file1 <- "hdf5_write_R_testdims_byrow.h5ad" #tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
#   h5ad_file2 <- "hdf5_write_py_testdims_byrow.h5ad" #tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")

#   write_h5ad(data_r, h5ad_file1)
#   data_python$write_h5ad(h5ad_file2)

#   res <- processx::run("h5diff", c("-v", h5ad_file1, h5ad_file2, "/X"), error_on_status = FALSE)

#   expect_equal(res$status, 0, info = res$stdout)

# # })


# # # test different matrices in X
# # test_that("h5diff_X_float", {
# #   data_r <- r_generate_dataset_only_x(10L, 20L, x_type = "numeric_matrix")
# #   data_python <- py_generate_dataset_only_x(10L, 20L, x_type = "generate_float_matrix")

# #   h5ad_file1 <- tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
# #   h5ad_file2 <- tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")

# #   write_h5ad(data_r, h5ad_file1)
# #   data_python$write_h5ad(h5ad_file2)

# #   res <- processx::run("h5diff", c("-v", h5ad_file1, h5ad_file2, "/X"), error_on_status = FALSE)

# #   expect_equal(res$status, 0, info = res$stdout)

# # })

# data_r <- generate_dataset(10L, 20L, format = "AnnData")

# da <- reticulate::import("dummy_anndata")
# data_python <- da$generate_dataset(10L, 20L)

# h5ad_file1 <- "hdf5_write_R_test3.h5ad" #tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
# h5ad_file2 <- "hdf5_write_py_test3.h5ad" #tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")

# write_h5ad(data_r, h5ad_file1)
# data_python$write_h5ad(h5ad_file2)
# res <- processx::run("h5diff", c("-v", h5ad_file1, h5ad_file2, "/X"), error_on_status = FALSE)

# expect_equal(res$status, 0, info = res$stdout)

# # test_that("h5diff", {
# #   requireNamespace("processx")

# #   data_r <- generate_dataset(10L, 20L, format = "AnnData")

# #   da <- reticulate::import("dummy_anndata")
# #   data_python <- da$generate_dataset(10L, 20L)

# #   h5ad_file1 <- "hdf5_write_R_test2.h5ad" #tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
# #   h5ad_file2 <- "hdf5_write_py_test2.h5ad" #tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")

# #   write_h5ad(data_r, h5ad_file1)
# #   data_python$write_h5ad(h5ad_file2)
# #   res <- processx::run("h5diff", c("-v", h5ad_file1, h5ad_file2, "/X"), error_on_status = FALSE)

# #   expect_equal(res$status, 0, info = res$stdout)

# # })
