# skip_if_not_installed("hdf5r")

# requireNamespace("reticulate")
# testthat::skip_if_not(
#   reticulate::py_module_available("dummy_anndata"),
#   message = "Python dummy_anndata module not available for testing"
# )

matrix_equivalences <- list(
  "float_matrix" = list("numeric_matrix"), #, "numeric_dense"), #numeric dense does dgematrix
  "float_matrix_nas" = list("numeric_matrix_with_nas"), #, "numeric_dense_with_nas"), #numeric dense does dgematrix
  "integer_matrix" = list("integer_matrix"),# , "integer_dense"),
  "float_csparse" = list("numeric_csparse"),
  "float_csparse_nas" = list("numeric_csparse_with_nas"),
  "float_rsparse" = list("numeric_rsparse"),
  "float_rsparse_nas" = list("numeric_rsparse_with_nas")
)

vector_equivalence <- list(
  "categorical" = list("factor"),
  "categorical_ordered" = list("factor_ordered"),
  "categorical_missing_values" = list("factor_with_nas"),
  "categorical_ordered_missing_values" = list("factor_ordered_with_nas"),
  "string_array" = list("character"),
  "dense_array" = list("numeric"),
  "integer_array" = list("integer"),
  "boolean_array" = list("logical"),
  "nullable_integer_array" = list("integer_with_nas"),
  "nullable_boolean_array" = list("logical_with_nas")
)

# TODO: check if processx is available
# TODO: check if h5diff is available --> hdf5-tools

da <- reticulate::import("dummy_anndata")

check_arg <- function(args, name, falseval) {
  if (name %in% names(args)) {
    args[[name]]
  } else {
    falseval
  }
}

py_generate_dataset <- function(n_obs, n_vars, write=FALSE, ...){
  args <- list(...)

  data <- da$generate_dataset(n_obs, n_vars,
                              x_type = check_arg(args, "x_type", NULL),
                              layer_types = check_arg(args, "layer_types", character()),
                              obs_types = ifelse("obs_types" %in% names(args), args$obs_types, list("integer_array")),
                              var_types = ifelse("var_types" %in% names(args), args$var_types, list("integer_array")),
                              obsm_types = check_arg(args, "obsm_types", character()),
                              varm_types = check_arg(args, "varm_types", character()),
                              obsp_types = check_arg(args, "obsp_types", character()),
                              varp_types = check_arg(args, "varp_types", character()),
                              uns_types = check_arg(args, "uns_types", character()),
                              nested_uns_types = check_arg(args, "nested_uns_types", character()))

  if (write) { 
    py_write_dataset(data)
  }
  data
}

py_write_dataset <- function(dataset, file=NULL){
  if (is.null(file)) {
    file <- tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")
  }
  dataset$write_h5ad(file)
}

r_generate_dataset <- function(n_obs, n_vars, write=FALSE, ...){
  args <- list(...)

  data <- generate_dataset(n_obs, n_vars,
                           x_type = check_arg(args, "x_type", "numeric_matrix"),
                           layer_types = check_arg(args, "layer_types", character()),
                           obs_types = ifelse("obs_types" %in% names(args), args$obs_types, "integer"),
                           var_types = ifelse("var_types" %in% names(args), args$var_types, "integer"),
                           obsm_types = check_arg(args, "obsm_types", character()),
                           varm_types = check_arg(args, "varm_types", character()),
                           obsp_types = check_arg(args, "obsp_types", character()),
                           varp_types = check_arg(args, "varp_types", character()),
                           uns_types = check_arg(args, "uns_types", character()),
                           format = "AnnData")
  if (write) {
    r_write_dataset(data)
  }

  data
}

r_write_dataset <- function(dataset, file=NULL){
  if (is.null(file)) {
    file <- tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
  }
  write_h5ad(dataset, file)
  file
}

# Test obs & var
# - there can be vectors in there

for (py_vector in names(vector_equivalence)) {
  data_python <- py_generate_dataset(10L, 20L, obs_types = list(py_vector))
  py_location <- py_write_dataset(data_python)

  for (r_vector in vector_equivalences[[py_matrix]]) {

      data_r <- r_generate_dataset(10L, 20L, obs_types = list(r_vector))
      r_location <- r_write_dataset(data_r)

      tryCatch({
        res <- processx::run("h5diff", c("-v", h5ad_file_py, h5ad_file_r, "/obs"), error_on_status = FALSE)
      }, error = function(e) {
        message("Error: ", e$message)
        message("Python matrix: ", py_matrix)
        message("R matrix: ", r_matrix)
      }, warning = function(w) {
        message("Warning: ", w$message)
        message("Python matrix: ", py_matrix)
        message("R matrix: ", r_matrix)
        
      })

  }
}



# for (py_matrix in names(matrix_equivalences)) {
#   data_python <- py_generate_dataset_only_x(10L, 20L, x_type = py_matrix)
#   h5ad_file_py <- tempfile(pattern = "hdf5_write_py_", fileext = ".h5ad")
#   data_python$write_h5ad(h5ad_file_py)

#   for (r_matrix in matrix_equivalences[[py_matrix]]) {
#     test_that(paste0("h5diff_X_", py_matrix, "_", r_matrix), {

#       tryCatch({
#         data_r <- r_generate_dataset_only_x(10L, 20L, x_type = r_matrix)
#         h5ad_file_r <- tempfile(pattern = "hdf5_write_R_", fileext = ".h5ad")
#         write_h5ad(data_r, h5ad_file_r)

#         res <- processx::run("h5diff", c("-v", h5ad_file_py, h5ad_file_r, "/X"), error_on_status = FALSE)

#         expect_equal(res$status, 0, info = res$stdout)
#       }, error = function(e) {
#         message("Error: ", e$message)
#         message("Python matrix: ", py_matrix)
#         message("R matrix: ", r_matrix)
        
#       }, warning = function(w) {
#         message("Warning: ", w$message)
#         message("Python matrix: ", py_matrix)
#         message("R matrix: ", r_matrix)
        
#       })
#     })
#   }
# }

# Mismatch <- R6::R6Class("Mismatch",
#   public = list(
#     Rgenerated = NULL,
#     Pygenerated = NULL,
#     errormsg = NULL,
#     initialize = function(message) {
#       self$message <- message
#     },
#     message = NULL
#   )
# )

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
