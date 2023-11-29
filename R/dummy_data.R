#' Dummy data
#'
#' Generate a dummy dataset in a selected format.
#'
#' @param n_obs Number of observations to generate.
#' @param n_var Number of variables to generate.
#' @param n_obsm Number of embedding dimensions.
#' @param n_varm Number of loading dimensions. 
#' @param output Object type to output, one of "list", "SingleCellExperiment",
#' or "Seurat".
#' @param density Of the matrices on a scale from 0-1.
#'  Passed to \link[Matrix]{rsparsematrix}.
#' @param ... Additional arguments passed to subfunctions
#'  (determined by \code{output}).
#' 
#'
#' @returns Object containing the generated dataset as defined by `output`
#'
#' @export
#' @examples
#' dummy <- dummy_data("InMemoryAnnData")
dummy_data <- function(output = c(
                         "list", 
                         "SingleCellExperiment",
                         "Seurat",
                         "InMemoryAnnData",
                         "HDF5AnnData"
                       ),
                       n_obs = 10L, 
                       n_var = 20L,
                       n_obsm = n_obs %/% 2,
                       n_varm = n_obsm, 
                       density=.1, 
                       ...) {
  output <- match.arg(output)  

  switch(output,
    "list" = dummy_list(n_obs = n_obs, 
                        n_var = n_var,
                        n_obsm = n_obsm,
                        n_varm = n_varm,
                        density = density),
    "SingleCellExperiment" = dummy_SingleCellExperiment(
      n_obs = n_obs,
      n_var = n_var, 
      n_obsm = n_obsm,
      n_varm = n_varm,
      density = density,
      ...
    ),
    "Seurat" = dummy_Seurat(n_obs = n_obs, 
                            n_var = n_var,
                            n_obsm = n_obsm,
                            n_varm = n_varm,
                            density = density,
                            ...),
    "InMemoryAnnData" = dummy_AnnData(n_obs = n_obs, 
                                      n_var = n_var,
                                      n_obsm = n_obsm,
                                      n_varm = n_varm,
                                      density = density,
                                      output_class = "InMemoryAnnData",
                                      ...),
    "HDF5AnnData" = dummy_AnnData(n_obs = n_obs,
                                  n_var = n_var, 
                                  n_obsm = n_obsm,
                                  n_varm = n_varm,
                                  density = density,
                                  output_class = "HDF5AnnData",
                                  ...)
  )
}

#' Dummy data list
#'
#' Generate a dummy dataset as a list
#'
#' @inheritParams dummy_data
#' @return A list with the generated dataset.
#' 
#' @keywords internal
#' @importFrom Matrix rsparsematrix
#' @importFrom stats runif
dummy_list <- function(n_obs = 10L, 
                       n_var = 20L,
                       n_obsm = n_obs %/% 2,
                       n_varm = n_obsm,
                       density = .1) {
  # generate X
  X <- Matrix::rsparsematrix(nrow = n_obs, 
                             ncol = n_var, 
                             density = density) 
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
    geneinfo = sample(c("a", "b", "c"), n_var, replace = TRUE)
  ) 
  # generate obsm
  obsm <- list(
    "X_pca"=matrix(runif(n_obs*n_obsm), n_obs, n_obsm),
    "X_umap"=matrix(runif(n_obs*n_obsm), n_obs, n_obsm)
    )
  # generate varm
  varm <- list(
    "PCs"=matrix(runif(n_var*n_varm), n_var, n_varm)
  )
  # generate obs_names
  obs_names <- paste0("cell", seq_len(n_obs))
  # generate var_names
  var_names <- paste0("gene", seq_len(n_var))
  # generate obsp
  obsp <- list(
    obs_graph=Matrix::rsparsematrix(n_obs, n_obs, density=density)
  )
  # generate varp
  varp <- list(
    var_graph=Matrix::rsparsematrix(n_var, n_var, density=density)
  )
  ## Return
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
    layers = layers
  )
}

#' Dummy SingleCellExperiment
#'
#' Generate a dummy dataset as a SingleCellExperiment object
#'
#' @inheritDotParams dummy_list
#'
#' @returns SingleCellExperiment containing the generated data
#' 
#' @keywords internal
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
    colData = dummy$obs,
    reducedDims = dummy$obsm
  )
  colnames(sce) <- dummy$obs_names
  rownames(sce) <- dummy$var_names 

  return(sce)
}

#' Dummy Seurat
#'
#' Generate a dummy dataset as a Seurat object
#'
#' @inheritDotParams dummy_list
#'
#' @returns Seurat containing the generated data
#' 
#' @keywords internal
#' @importFrom Matrix t
dummy_Seurat <- function(...,
                         assay="RNA") { #nolint 
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Creating a Seurat requires the 'SeuratObject' package to be installed"
    )
  }

  dummy <- dummy_data(...)

  X <- Matrix::t(dummy$X)
  colnames(X) <- dummy$obs_names
  rownames(X) <- dummy$var_names

  seurat <- SeuratObject::CreateSeuratObject(X, 
                                             assay = assay)

  X2 <- t(dummy$layers$X2)
  colnames(X2) <- dummy$obs_names
  rownames(X2) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(object = seurat,
                                       layer =  "data",
                                       new.data = X2,
                                       assay = assay)

  X3 <- as.matrix(Matrix::t(dummy$layers$X3))
  colnames(X3) <- dummy$obs_names
  rownames(X3) <- dummy$var_names
  seurat <- SeuratObject::SetAssayData(object = seurat, 
                                       layer = "scale.data",
                                       new.data = X3, 
                                       assay = assay)

  seurat <- SeuratObject::AddMetaData(object = seurat,
                                      metadata = dummy$obs)

  ## Add DimReduc objects
  drl <- to_DimReduc(obj = dummy, assay=assay)
  seurat <- add_DimReduc(obj = seurat, 
                         drl = drl)
  ## Add obsp/varp as graphs
  seurat@graphs <- dummy[c("obsp","varp")]
  return(seurat)
}



#' Dummy anndata
#'
#' Generate a dummy dataset as a \link[anndataR]{InMemoryAnnData} or 
#' \link[anndataR]{HDF5AnnData} object
#' 
#' @param output_class Name of the AnnData class. 
#' Must be one of:
#' \itemize{
#' \item{"InMemoryAnnData": }{
#' Produces \link[anndataR]{InMemoryAnnData} (default).}
#' \item{"HDF5AnnData": }{
#' Produces \link[anndataR]{HDF5AnnData}.}
#' }
#' @param file Path to save object to 
#' (only used when `output_class="HDF5AnnData"`).
#' 
#' @inheritDotParams dummy_list
#'
#' @returns \link[anndataR]{InMemoryAnnData} or 
#' \link[anndataR]{HDF5AnnData} containing the generated data. 
#' 
#' @keywords internal
dummy_AnnData <- function(output_class=c("InMemoryAnnData",
                                         "HDF5AnnData"),
                          file=tempfile(fileext = ".h5ad"),
                          ...) { #nolint 
  
  dummy <- dummy_data(...)
  output_class <- output_class[1]
  generator <- get_generator(output_class)
  if(output_class=="HDF5AnnData"){
    ad <- generator$new(
      X = dummy$X,
      # layers = dummy$layers,
      obs = dummy$obs,
      obs_names = dummy$obs_names,
      obsm = dummy$obsm,
      obsp = dummy$obsp,
      var = dummy$var,
      var_names = dummy$var_names,
      varm = dummy$varm,
      varp = dummy$varp,
      file=file
    )
  } else {
    ad <- generator$new(
      X = dummy$X,
      # layers = dummy$layers,
      obs = dummy$obs,
      obs_names = dummy$obs_names,
      obsm = dummy$obsm,
      obsp = dummy$obsp,
      var = dummy$var,
      var_names = dummy$var_names,
      varm = dummy$varm,
      varp = dummy$varp,
    )
  } 
  ## Set layers
  # Do this in a second step as it seems to cause errors when setting 
  # during initialisation
  ad$layers <- dummy$layers
  ## Return
  return(ad)
}

