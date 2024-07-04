#' Generate a dataset
#'
#' Generate a dataset with different types of columns and layers
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#' @param x_type Type of matrix to generate for X
#' @param layer_types Types of matrices to generate for layers
#' @param obs_types Types of vectors to generate for obs
#' @param var_types Types of vectors to generate for var
#' @param obsm_types Types of matrices to generate for obsm
#' @param varm_types Types of matrices to generate for varm
#' @param obsp_types Types of matrices to generate for obsp
#' @param varp_types Types of matrices to generate for varp
#' @param uns_types Types of objects to generate for uns
#' @param example If `TRUE`, the types will be overridden to a small set of
#'   types. This is useful for documentations.
#' @param format Object type to output, one of "list", "AnnData",
#'   "SingleCellExperiment", or "Seurat".
#'
#' @return Object containing the generated dataset as defined by `output`
#'
#' @export
#'
#' @examples
#' dummy <- generate_dataset()
#' \dontrun{
#' dummy <- generate_dataset(format = "AnnData")
#' dummy <- generate_dataset(format = "SingleCellExperiment")
#' dummy <- generate_dataset(format = "Seurat")
#' }
generate_dataset <- function(
    n_obs = 10L,
    n_vars = 20L,
    x_type = "numeric_matrix",
    layer_types = c(
      "numeric_matrix", "numeric_dense", "numeric_csparse", "numeric_rsparse", "numeric_matrix_with_nas", #
      "numeric_dense_with_nas", "numeric_csparse_with_nas", "numeric_rsparse_with_nas", "integer_matrix",
      "integer_dense", "integer_csparse", "integer_rsparse", "integer_matrix_with_nas", "integer_dense_with_nas",
      "integer_csparse_with_nas", "integer_rsparse_with_nas"
    ),
    obs_types = c(
      "character", "integer", "factor", "factor_ordered", "logical", "numeric", "character_with_nas",
      "integer_with_nas", "factor_with_nas", "factor_ordered_with_nas", "logical_with_nas", "numeric_with_nas"
    ),
    var_types = c(
      "character", "integer", "factor", "factor_ordered", "logical", "numeric", "character_with_nas",
      "integer_with_nas", "factor_with_nas", "factor_ordered_with_nas", "logical_with_nas", "numeric_with_nas"
    ),
    obsm_types = c(
      "numeric_matrix", "numeric_dense", "numeric_csparse", "numeric_rsparse", "numeric_matrix_with_nas",
      "numeric_dense_with_nas", "numeric_csparse_with_nas", "numeric_rsparse_with_nas", "integer_matrix",
      "integer_dense", "integer_csparse", "integer_rsparse", "integer_matrix_with_nas", "integer_dense_with_nas",
      "integer_csparse_with_nas", "integer_rsparse_with_nas", "character", "integer", "factor", "factor_ordered",
      "logical", "numeric", "character_with_nas", "integer_with_nas", "factor_with_nas", "factor_ordered_with_nas",
      "logical_with_nas", "numeric_with_nas"
    ),
    varm_types = c(
      "numeric_matrix", "numeric_dense", "numeric_csparse", "numeric_rsparse", "numeric_matrix_with_nas",
      "numeric_dense_with_nas", "numeric_csparse_with_nas", "numeric_rsparse_with_nas", "integer_matrix",
      "integer_dense", "integer_csparse", "integer_rsparse", "integer_matrix_with_nas", "integer_dense_with_nas",
      "integer_csparse_with_nas", "integer_rsparse_with_nas", "character", "integer", "factor", "factor_ordered",
      "logical", "numeric", "character_with_nas", "integer_with_nas", "factor_with_nas", "factor_ordered_with_nas",
      "logical_with_nas", "numeric_with_nas"
    ),
    obsp_types = c(
      "numeric_matrix", "numeric_dense", "numeric_csparse", "numeric_rsparse", "numeric_matrix_with_nas",
      "numeric_dense_with_nas", "numeric_csparse_with_nas", "numeric_rsparse_with_nas", "integer_matrix",
      "integer_dense", "integer_csparse", "integer_rsparse", "integer_matrix_with_nas", "integer_dense_with_nas",
      "integer_csparse_with_nas", "integer_rsparse_with_nas"
    ),
    varp_types = c(
      "numeric_matrix", "numeric_dense", "numeric_csparse", "numeric_rsparse", "numeric_matrix_with_nas",
      "numeric_dense_with_nas", "numeric_csparse_with_nas", "numeric_rsparse_with_nas", "integer_matrix",
      "integer_dense", "integer_csparse", "integer_rsparse", "integer_matrix_with_nas", "integer_dense_with_nas",
      "integer_csparse_with_nas", "integer_rsparse_with_nas"
    ),
    uns_types = c(
      "scalar_character", "scalar_integer", "scalar_factor", "scalar_factor_ordered", "scalar_logical",
      "scalar_numeric", "scalar_character_with_nas", "scalar_integer_with_nas", "scalar_factor_with_nas",
      "scalar_factor_ordered_with_nas", "scalar_logical_with_nas", "scalar_numeric_with_nas", "vec_character",
      "vec_integer", "vec_factor", "vec_factor_ordered", "vec_logical", "vec_numeric",
      "vec_character_with_nas", "vec_integer_with_nas", "vec_factor_with_nas",
      "vec_factor_ordered_with_nas", "vec_logical_with_nas", "vec_numeric_with_nas",
      "df_character", "df_integer", "df_factor", "df_factor_ordered",
      "df_logical", "df_numeric", "df_character_with_nas", "df_integer_with_nas",
      "df_factor_with_nas", "df_factor_ordered_with_nas", "df_logical_with_nas",
      "df_numeric_with_nas", "mat_numeric_matrix", "mat_numeric_dense", "mat_numeric_csparse", "mat_numeric_rsparse",
      "mat_numeric_matrix_with_nas", "mat_numeric_dense_with_nas", "mat_numeric_csparse_with_nas",
      "mat_numeric_rsparse_with_nas", "mat_integer_matrix", "mat_integer_dense", "mat_integer_csparse",
      "mat_integer_rsparse", "mat_integer_matrix_with_nas", "mat_integer_dense_with_nas",
      "mat_integer_csparse_with_nas", "mat_integer_rsparse_with_nas"
    ),
    example = FALSE,
    format = c("list", "AnnData", "SingleCellExperiment", "Seurat")) {
  format <- match.arg(format)

  if (example) {
    x_type <- "numeric_matrix"
    layer_types <- c("numeric_matrix", "numeric_dense", "numeric_csparse")
    obs_types <- c("character", "integer", "factor")
    var_types <- c("character", "integer", "factor")
    obsm_types <- c("numeric_matrix", "numeric_dense", "numeric_csparse")
    varm_types <- c("numeric_matrix", "numeric_dense", "numeric_csparse")
    obsp_types <- c("numeric_matrix", "numeric_dense", "numeric_csparse")
    varp_types <- c("numeric_matrix", "numeric_dense", "numeric_csparse")
    uns_types <- c("scalar_character", "vec_integer", "df_logical", "mat_numeric_matrix")
  }

  dataset_list <- .generate_dataset_as_list(
    n_obs = n_obs,
    n_vars = n_vars,
    x_type = x_type,
    layer_types = layer_types,
    obs_types = obs_types,
    var_types = var_types,
    obsm_types = obsm_types,
    varm_types = varm_types,
    obsp_types = obsp_types,
    varp_types = varp_types,
    uns_types = uns_types
  )

  conversion_fun <- switch(format,
    "list" = identity,
    "SingleCellExperiment" = .generate_dataset_as_sce,
    "Seurat" = .generate_dataset_as_seurat,
    "AnnData" = .generate_dataset_as_anndata
  )

  return(conversion_fun(dataset_list))
}

#' Generate a dummy dataset as a list
#'
#' @inheritParams generate_dataset
#'
#' @return A list with the generated dataset
#'
#' @noRd
.generate_dataset_as_list <- function(
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
    uns_types = c(
      paste0("scalar_", names(vector_generators)),
      paste0("vec_", names(vector_generators)),
      paste0("df_", names(vector_generators)),
      paste0("mat_", names(matrix_generators))
    )) {
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
  rownames(obs) <- obs_names

  # generate var_names
  var_names <- paste0("gene", seq_len(n_vars))
  rownames(var) <- var_names

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
      generate_matrix(n_vars, n_vars = 10L, varm_type)
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
  uns <- lapply(uns_types, function(uns_type) {
    if (uns_type == "list") {
      # TODO: add multiple data types here
      list(1L, 1, "a", factor("a"), TRUE, list(nested = list(list = c(1, 2, 3))))
    } else if (uns_type == "scalar_character_with_nas") {
      NA_character_
    } else if (uns_type == "scalar_integer_with_nas") {
      NA_integer_
    } else if (uns_type == "scalar_factor_with_nas") {
      factor(NA_character_)
    } else if (uns_type == "scalar_factor_ordered_with_nas") {
      factor(NA_character_, ordered = TRUE)
    } else if (uns_type == "scalar_logical_with_nas") {
      NA_real_
    } else if (grepl("scalar_", uns_type)) {
      generate_vector(1L, gsub("scalar_", "", uns_type))
    } else if (grepl("vec_", uns_type)) {
      generate_vector(10L, gsub("vec_", "", uns_type))
    } else if (grepl("df_", uns_type)) {
      generate_dataframe(10L, gsub("df_", "", uns_type))
    } else if (grepl("mat_", uns_type)) {
      generate_matrix(10L, 10L, gsub("mat_", "", uns_type))
    } else {
      stop("Unknown uns type: '", uns_type, "'")
    }
  })
  names(uns) <- uns_types

  # return list
  list(
    X = X,
    obs = obs,
    obsm = obsm,
    obsp = obsp,
    var = var,
    varm = varm,
    varp = varp,
    layers = layers,
    uns = uns
  )
}

#' Convert a dummy dataset to a SingleCellExperiment object
#'
#' @param dataset_list Output of `.generate_dataset_as_list()`
#'
#' @return SingleCellExperiment containing the generated data
#'
#' @noRd
.generate_dataset_as_sce <- function(dataset_list) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop(
      "Creating a SingleCellExperiment requires the 'SingleCellExperiment'",
      "package to be installed"
    )
  }

  assays_list <- c(
    list(X = dataset_list$X),
    dataset_list$layers
  )
  assays_list <- lapply(assays_list, Matrix::t)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    rowData = dataset_list$var,
    colData = dataset_list$obs
  )
  colnames(sce) <- dataset_list$obs_names
  rownames(sce) <- dataset_list$var_names

  # TODO: add obsm, varm, obsp, varp, uns?

  return(sce)
}

#' Convert a dummy dataset to a Seurat object
#'
#' @param dataset_list Output of `.generate_dataset_as_list()`
#'
#' @return Seurat containing the generated data
#'
#' @noRd
.generate_dataset_as_seurat <- function(dataset_list) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Creating a Seurat requires the 'SeuratObject' package to be installed"
    )
  }

  X <- t(dataset_list$layers[["integer_csparse"]])
  colnames(X) <- dataset_list$obs_names
  rownames(X) <- dataset_list$var_names

  seurat <- SeuratObject::CreateSeuratObject(X)

  X2 <- Matrix::t(dataset_list$layers[["numeric_csparse"]])
  colnames(X2) <- dataset_list$obs_names
  rownames(X2) <- dataset_list$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "data", X2)

  X3 <- Matrix::t(dataset_list$layers[["numeric_matrix"]])
  colnames(X3) <- dataset_list$obs_names
  rownames(X3) <- dataset_list$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "scale.data", X3)

  # TODO: Seurat v5 now supports more than just these three layers

  seurat <- SeuratObject::AddMetaData(seurat, dataset_list$obs)

  # TODO: add obsm, varm, obsp, varp, uns?

  return(seurat)
}

#' Convert a dummy dataset to an AnnData object
#'
#' @param list Output of `.generate_dataset_as_list()`
#'
#' @return SingleCellExperiment containing the generated data
#'
#' @noRd
.generate_dataset_as_anndata <- function(dataset_list) { # nolint
  AnnData(
    X = dataset_list$X,
    obs = dataset_list$obs,
    obsm = dataset_list$obsm,
    obsp = dataset_list$obsp,
    var = dataset_list$var,
    varm = dataset_list$varm,
    varp = dataset_list$varp,
    layers = dataset_list$layers,
    uns = dataset_list$uns
  )
}
