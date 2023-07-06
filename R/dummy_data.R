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
#' @export
#' 
#' @examples
#' dummy <- dummy_data()
dummy_data <- function(n_obs = 10L, n_vars = 20L,
                            output = c("list", "SingleCellExperiment",
                                       "Seurat")) {
  output <- match.arg(output)
  
  switch (output,
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

  # generate layers
  layers <- list(
    X2 = X * 2,
    X3 = X * 3
  )

  # generate obs
  obs <- data.frame(
    cell_type = sample(c("tcell", "bcell"), n_obs, replace = TRUE),
    cluster = sample.int(3, n_obs, replace = TRUE)
  )

  # generate var
  var <- data.frame(
    geneinfo = sample(c("a", "b", "c"), n_vars, replace = TRUE)
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
dummy_SingleCellExperiment <- function(...) {
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop(
      "Creating a SingleCellExperiment requires the 'SingleCellExperiment'",
      "package to be installed"
    )
  }
  
  dummy <- dummy_data(...)
  
  assays_list = c(
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
dummy_Seurat <- function(...) {
  
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Creating a Seurat requires the 'SeuratObject' package to be installed"
    )
  }
  
  dummy <- dummy_data(...)

  X <- t(dummy$X)
  colnames(X) <- dummy$obs_names
  rownames(X) <- dummy$var_names
  
  seurat <- SeuratObject::CreateSeuratObject(X)
  
  X2 <- t(dummy$layers$X2)
  colnames(X2) <- dummy$obs_names
  rownames(X2) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "data", X2)
  
  X3 <- as.matrix(t(dummy$layers$X3))
  colnames(X3) <- dummy$obs_names
  rownames(X3) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(seurat, "scale.data", X3)
  
  seurat <- SeuratObject::AddMetaData(seurat, dummy$obs)
  
  return(seurat)
}
