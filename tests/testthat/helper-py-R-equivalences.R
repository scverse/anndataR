# python, R
matrix_equivalences <- list(
    c("float_matrix", "numeric_matrix"),
    # c("float_matrix", "numeric_dense"),
    c("float_matrix_nas", "numeric_matrix_with_nas"),
    # c("float_matrix_nas", "numeric_dense_with_nas"),
    c("integer_matrix", "integer_matrix"),
    # c("integer_matrix", "integer_dense"),
    c("float_csparse", "numeric_csparse"),
    c("float_csparse_nas", "numeric_csparse_with_nas"),
    c("float_rsparse", "numeric_rsparse"),
    c("float_rsparse_nas", "numeric_rsparse_with_nas")
)

# python, R
vector_equivalences <- list(
    c("categorical", "factor"),
    c("categorical_ordered", "factor_ordered"),
    c("categorical_missing_values", "factor_with_nas"),
    c("categorical_ordered_missing_values", "factor_ordered_with_nas"),
    c("string_array", "character"),
    c("dense_array", "numeric"),
    c("integer_array", "integer"),
    c("boolean_array", "logical"),
    c("nullable_integer_array", "integer_with_nas"),
    c("nullable_boolean_array", "logical_with_nas")
)

all_equivalences <- c(matrix_equivalences, vector_equivalences)

check_arg <- function(args, name, falseval) {
  if (name %in% names(args)) {
    args[[name]][[1]]
  } else {
    falseval
  }
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