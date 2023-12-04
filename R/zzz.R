# The defaults of the generate_dataset() function are set here.
# This is done here to ensure all helper functions are available
# by the time the defaults are set.
formals(generate_dataset) <- list(
  n_obs = 10L,
  n_vars = 20L,
  x_type = names(matrix_generators)[[1]],
  layer_types = names(matrix_generators),
  obs_types = names(vector_generators),
  var_types = names(vector_generators),
  obsm_types = c(names(matrix_generators), names(vector_generators)),
  varm_types = c(names(matrix_generators), names(vector_generators)),
  obsp_types = names(matrix_generators),
  varp_types = names(matrix_generators),
  format = c("list", "AnnData", "SingleCellExperiment", "Seurat")
)