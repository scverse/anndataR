#' Generate a dataset
#'
#' Generate a dataset with different types of columns and layers
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#' @param format Object type to output, one of "list", "SingleCellExperiment",
#' or "Seurat"
#' @param ... Arguments passed to generate_dataset_as_list
#'
#' @return Object containing the generated dataset as defined by `output`
#' 
#' @noRd
#'
#' @examples
#' dummy <- generate_dataset()
generate_dataset <- function(
  n_obs = 10L,
  n_vars = 20L,
  format = c("list", "SingleCellExperiment", "Seurat"),
  ...
) {
  format <- match.arg(format)

  fun <- switch(format,
    "list" = generate_dataset_as_list,
    "SingleCellExperiment" = generate_dataset_as_sce,
    "Seurat" = generate_dataset_as_seurat
  )

  fun(n_obs = n_obs, n_vars = n_vars, ...)
}

#' Dummy data list
#'
#' Generate a dummy dataset as a list
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#'
#' @return A list with the generated dataset
#' 
#' @noRd
generate_dataset_as_list <- function(
  n_obs = 10L,
  n_vars = 20L,
  x_type = names(matrix_generators)[[1]],
  layer_types = names(matrix_generators),
  obs_types = names(vector_generators),
  var_types = names(vector_generators),
  obsm_types = c(names(matrix_generators), names(vector_generators)),
  varm_types = c(names(matrix_generators), names(vector_generators)),
  obsp_types = names(matrix_generators),
  varp_types = names(matrix_generators)
) {
  # generate X
  X <- generate_matrix(n_obs, n_vars, x_type)

  # generate layers
  layers <- lapply(layer_types, generate_matrix, n_obs = n_obs, n_vars = n_vars)
  names(layers) <- layer_types

  # generate obs
  obs <- generate_dataframe(n_obs, obs_types)

  # generate var
  var <- generate_dataframe(n_vars, var_types)

  # generate obs_names
  obs_names <- paste0("cell", seq_len(n_obs))

  # generate var_names
  var_names <- paste0("gene", seq_len(n_vars))

  # generate obsm
  obsm <- lapply(obsm_types, function(obsm_type) {
    if (obsm_type %in% names(vector_generators)) {
      generate_dataframe(n_obs, obsm_type)
    } else {
      generate_matrix(n_obs, n_vars = 10L, obsm_type)
    }
  })
  names(obsm) <- obsm_types

  # generate varm
  varm <- lapply(varm_types, function(varm_type) {
    if (varm_type %in% names(vector_generators)) {
      generate_dataframe(n_vars, varm_type)
    } else {
      generate_matrix(n_vars, n_obs = 10L, varm_type)
    }
  })
  names(varm) <- varm_types

  # generate obsp
  obsp <- lapply(obsp_types, generate_matrix, n_obs = n_obs, n_vars = n_obs)
  names(obsp) <- obsp_types

  # generate varp
  varp <- lapply(varp_types, generate_matrix, n_obs = n_vars, n_vars = n_vars)
  names(varp) <- varp_types

  # generate uns by combining other classes
  uns <- list(
    integer = 1L,
    numeric = 1,
    character = "a",
    factor = factor("a"),
    logical = TRUE,
    integer_na = NA_integer_,
    numeric_na = NA_real_,
    character_na = NA_character_,
    factor_na = NA_character_,
    logical_na = NA,
    list = list(1L, 1, "a", factor("a"), TRUE)
  )
  vectors_for_uns <- lapply(names(vector_generators), generate_vector, n = 10L)
  names(vectors_for_uns) <- paste0("vec_", names(vector_generators))
  obsm_for_uns <- obsm
  names(obsm_for_uns) <- paste0("obsm_", names(obsm_for_uns))

  uns <- c(
    uns,
    vectors_for_uns,
    obsm_for_uns
  )

  list(
    X = X,
    obs = obs,
    obs_names = obs_names,
    obsm = obsm,
    obsp = obsp,
    var = var,
    var_names = var_names,
    varm = varm,
    varp = varp,
    layers = layers,
    uns = uns
  )
}

#' Dummy SingleCellExperiment
#'
#' Generate a dummy dataset as a SingleCellExperiment object
#'
#' @param ... Parameters passed to `generate_dataset_as_list`
#'
#' @return SingleCellExperiment containing the generated data
#' 
#' @noRd
generate_dataset_as_sce <- function(...) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop(
      "Creating a SingleCellExperiment requires the 'SingleCellExperiment'",
      "package to be installed"
    )
  }

  dummy <- generate_dataset_as_list(...)

  assays_list <- c(
    list(X = dummy$X),
    dummy$layers
  )
  assays_list <- lapply(assays_list, Matrix::t)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    rowData = dummy$var,
    colData = dummy$obs
  )
  colnames(sce) <- dummy$obs_names
  rownames(sce) <- dummy$var_names

  return(sce)
}

#' Dummy Seurat
#'
#' Generate a dummy dataset as a Seurat object
#'
#' @param ... Parameters passed to `generate_dataset_as_list`
#'
#' @return Seurat containing the generated data
generate_dataset_as_seurat <- function(...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Creating a Seurat requires the 'SeuratObject' package to be installed"
    )
  }

  dummy <- generate_dataset_as_list(...)

  X <- t(dummy$layers[["integer_csparse"]])
  colnames(X) <- dummy$obs_names
  rownames(X) <- dummy$var_names

  seurat <- SeuratObject::CreateSeuratObject(X)

  X2 <- Matrix::t(dummy$layers[["numeric_csparse"]])
  colnames(X2) <- dummy$obs_names
  rownames(X2) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "data", X2)

  X3 <- Matrix::t(dummy$layers[["numeric_matrix"]])
  colnames(X3) <- dummy$obs_names
  rownames(X3) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "scale.data", X3)

  seurat <- SeuratObject::AddMetaData(seurat, dummy$obs)

  return(seurat)
}
