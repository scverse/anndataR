#' Dummy data
#'
#' Generate a dummy dataset
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#' @param output Object type to output, one of "list", "SingleCellExperiment",
#' or "Seurat"
#'
#' @return Object containing the generated dataset as defined by `output`
#'
#' @examples
#' dummy <- dummy_data()
dummy_data <- function(
    n_obs = 10L,
    n_vars = 20L,
    output = c(
      "list", "SingleCellExperiment", "Seurat"
    )) {
  output <- match.arg(output)

  switch(output,
    "list" = dummy_list(n_obs = n_obs, n_vars = n_vars),
    "SingleCellExperiment" = dummy_SingleCellExperiment(
      n_obs = n_obs, n_vars = n_vars
    ),
    "Seurat" = dummy_Seurat(n_obs = n_obs, n_vars = n_vars)
  )
}

#' Dummy data list
#'
#' Generate a dummy dataset as a list
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#'
#' @return A list with the generated dataset
dummy_list <- function(n_obs = 10L, n_vars = 20L) {
  # generate X
  X <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)

  # generate layers in dense, sparse, row compressed and column compressed
  layers <- list(
    dense = as.matrix(X),
    CsparseMatrix = as(X, "CsparseMatrix"),
    RsparseMatrix = as(X, "RsparseMatrix")
  )

  # generate obs with different types (character, integer, factors, etc.)
  obs <- data.frame(
    character = paste0("cell", seq_len(n_obs)),
    integer = seq_len(n_obs),
    factor = factor(paste0("cell", seq_len(n_obs))),
    factor_ordered = factor(paste0("cell", seq_len(n_obs)), ordered = TRUE),
    logical = sample(c(TRUE, FALSE), n_obs, replace = TRUE),
    numeric = runif(n_obs),
    character_with_nas = c(paste0("cell", seq_len(n_obs - 1)), NA),
    integer_with_nas = c(seq_len(n_obs - 1), NA),
    factor_with_nas = c(
      factor(paste0("cell", seq_len(n_obs - 1))),
      NA_character_
    ),
    factor_ordered_with_nas = c(
      factor(paste0("cell", seq_len(n_obs - 1)), ordered = TRUE),
      NA_character_
    ),
    logical_with_nas = c(sample(c(TRUE, FALSE), n_obs - 1, replace = TRUE), NA),
    numeric_with_nas = c(runif(n_obs - 1), NA)
  )

  # generate var with different types (character, integer, factors, etc.)
  var <- data.frame(
    character = paste0("gene", seq_len(n_vars)),
    integer = seq_len(n_vars),
    factor = factor(paste0("gene", seq_len(n_vars))),
    factor_ordered = factor(paste0("gene", seq_len(n_vars)), ordered = TRUE),
    logical = sample(c(TRUE, FALSE), n_vars, replace = TRUE),
    numeric = runif(n_vars),
    character_with_nas = c(paste0("gene", seq_len(n_vars - 1)), NA),
    integer_with_nas = c(seq_len(n_vars - 1), NA),
    factor_with_nas = c(
      factor(paste0("gene", seq_len(n_vars - 1))),
      NA_character_
    ),
    factor_ordered_with_nas = c(
      factor(paste0("gene", seq_len(n_vars - 1)), ordered = TRUE),
      NA_character_
    ),
    logical_with_nas = c(sample(c(TRUE, FALSE), n_vars - 1, replace = TRUE), NA),
    numeric_with_nas = c(runif(n_vars - 1), NA)
  )

  # generate obs_names
  obs_names <- paste0("cell", seq_len(n_obs))

  # generate var_names
  var_names <- paste0("gene", seq_len(n_vars))

  list(
    X = X,
    obs = obs,
    obs_names = obs_names,
    var = var,
    var_names = var_names,
    layers = layers
  )
}

#' Dummy SingleCellExperiment
#'
#' Generate a dummy dataset as a SingleCellExperiment object
#'
#' @param ... Parameters passed to `dummy_list`
#'
#' @return SingleCellExperiment containing the generated data
dummy_SingleCellExperiment <- function(...) { # nolint
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop(
      "Creating a SingleCellExperiment requires the 'SingleCellExperiment'",
      "package to be installed"
    )
  }

  dummy <- dummy_data(...)

  assays_list <- c(
    list(X = dummy$X),
    dummy$layers
  )
  assays_list <- lapply(assays_list, t)

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
#' @param ... Parameters passed to `dummy_list`
#'
#' @return Seurat containing the generated data
dummy_Seurat <- function(...) { # nolint
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Creating a Seurat requires the 'SeuratObject' package to be installed"
    )
  }

  dummy <- dummy_data(...)

  X <- t(dummy$layers[["dense"]])
  colnames(X) <- dummy$obs_names
  rownames(X) <- dummy$var_names

  seurat <- SeuratObject::CreateSeuratObject(X)

  X2 <- as(t(dummy$layers[["CsparseMatrix"]]), "CsparseMatrix")
  colnames(X2) <- dummy$obs_names
  rownames(X2) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "data", X2)

  # seurat doesn't support RsparseMatrices
  X3 <- as.matrix(t(dummy$layers[["RsparseMatrix"]]))
  colnames(X3) <- dummy$obs_names
  rownames(X3) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "scale.data", X3)

  seurat <- SeuratObject::AddMetaData(seurat, dummy$obs)

  return(seurat)
}
