#' Convert between \code{AnnData} and \code{DimReducObject}
#'
#' Create a named list of \link[SeuratObject]{DimReduc} objects from an
#'  \pkg{AnnData}. Alternatively, can supply arbitrary matrices or
#'  lists of matrices to the \code{obsm}/\code{varm} manually.
#' @source Some internal code adapted from 
#' \href{https://github.com/PMBio/MuDataSeurat/blob/main/R/ReadH5MU.R}{
#' MuDataSeurat}.
#' @param obj An object of class \link[anndataR]{InMemoryAnnData} or 
#' \link[anndataR]{HDF5AnnData}.
#' @param obsm A named list of observation (cell) embedding matrices.
#' @param varm A named list of variable (gene) loading matrices.
#' @param obs_names Names of all observations.
#' @param var_names Names of all variables.
#' @param key_map A key:value mapping indicating pairs of obsm/varm slots in
#'  the \code{obj} object. For some common reductions,
#'  there are conventional names for the loadings slots
#' @inheritDotParams SeuratObject::CreateDimReducObject
#' @returns A named list of \link[SeuratObject]{CreateDimReducObject} objects.
#'
#' @export
#' @examples
#' adata <- dummy_data("InMemoryAnnData")
#' drl <- to_DimReduc(obj=adata)
#' drl
to_DimReduc <- function(obj = NULL,
                        obsm = obj$obsm,
                        varm = obj$varm,
                        obs_names = obj$obs_names,
                        var_names = obj$var_names,
                        key_map= list("X_pca" = "PCs", 
                                      "X_umap" = NULL,
                                      "X_mofa" = "LFs"), 
                        ...) {
  
  requireNamespace("SeuratObject") 
  
  ## Need to at least have obsm
  if(is.null(obsm)) {
    messager("obsm is NULL; returning NULL.")
    return(NULL)
  }
  if(!is.list(obsm) || is.null(names(obsm))){
    messager("obsm name not provided; setting name to 'obsm'.")
    obsm <- list(obsm=obsm)
  }
  if(!is.list(varm) || is.null(names(varm))){
    messager("varm name not provided; setting name to 'varm'.")
    varm <- list(varm=varm)
  }  
  
  nms_msg <- function(X,
                      X_name,
                      X_type="loadings",
                      dim=2){
    dim_nm <- if(dim==1) "rownames" else "colnames"
    messager("No",dim_nm,"present in",paste0(X_type,"; setting to"),
             shQuote(paste0(X_name,"_[1:",ncol(X),"]"))
    ) 
  }
  # Add embeddings
  lapply(stats::setNames(names(obsm),
                         names(obsm)),
         function(key){
    messager("Preparing DimReduc for:",shQuote(key))              
    ## Prepare obsm
    obsm_nm2 <- gsub('X_', '', key)
    if(!key %in% names(key_map)){
      key_map[[key]] <- key
    }
    embeddings <- obsm[[key]]
    if(is.null(rownames(embeddings))){ 
      if(is.null(obs_names)) stop("Must provide `obs_names`.")
      rownames(embeddings) <- obs_names
    }
    if(is.null(colnames(embeddings))){ 
      nms_msg(X = embeddings,
              X_name = obsm_nm2, 
              X_type = "embeddings")
      colnames(embeddings) <- paste(obsm_nm2,seq(ncol(embeddings)),sep="_")  
    }
    embeddings <- as.matrix(embeddings)
    ## Prepare matching varm
    loadings <- matrix()
    value <- key_map[[key]]
    if (key %in% names(key_map) &&
        !is.null(value)) {
      if (value %in% names(varm)) {
        loadings <- varm[[value]]
        if(is.null(rownames(loadings))){
          if(is.null(var_names)) stop("Must provide `var_names`.")
          rownames(loadings) <- var_names
        }
        if(is.null(colnames(loadings))){
          nms_msg(X = embeddings,
                  X_name = obsm_nm2, 
                  X_type = "embeddings")
          colnames(loadings) <- paste(obsm_nm2,seq(ncol(loadings)),sep="_")  
        }
      }
    }
    ## Prepare stdev
    emb_stdev <- numeric()
    if ("uns" %in% names(obj)) {
      if (obsm_nm2 %in% names(obj[["uns"]])) {
        if ("variance" %in% names(obj[["uns"]][[obsm_nm2]])) {
          emb_stdev <- sqrt(obj[["uns"]][[obsm_nm2]][["variance"]])
        }
      }
    }
    SeuratObject::CreateDimReducObject(
      embeddings = as.matrix(embeddings),
      loadings = as.matrix(loadings),
      key = paste0(obsm_nm2, "_"), 
      stdev = emb_stdev,
      ...
    )
  })
}

#' Convert between \code{AnnData} and \code{DimReducObject}
#'
#' Create  obsm and varm AnnData objects from a named list of 
#' \link[SeuratObject]{DimReduc} object.
#' @param obj A named list of \link[SeuratObject]{DimReduc} objects.
#' Alternatively, can be a \link[SeuratObject]{SeuratObject} from which a 
#' named list of \link[SeuratObject]{DimReduc} objects will be extracted.
#' @param rm_rownames Remove rownames from all matrices.
#' @param rm_colnames Remove colnames from all matrices.
#' @param drop_null Drop elements that are NULL.
#' @returns A nested named list of \code{obsm} and \code{varm} matrices.
#' @export
#' @examples 
#' adata <- dummy_data("InMemoryAnnData")
#' obj <- to_DimReduc(obj=adata)
#' obsm_varm <- from_DimReduc(obj)
from_DimReduc <- function(obj,
                          rm_rownames=TRUE,
                          rm_colnames=TRUE,
                          drop_null=TRUE
                          ){
  
  requireNamespace("SeuratObject")

  ## Handle multiple object types
  if(methods::is(obj,"DimReduc")){
    stop("Must provide DimReduc objects as a named list.")
  } else if(SeuratObject::IsNamedList(obj)){
    if(!all(mapply(methods::is,"DimReduc"),obj)){
      stop("All elements of named list must be of class 'DimReduc'.")
    }
    drl <- obj
  } else if(methods::is(obj,"Seurat")){
    drl <- obj@reductions
  }
  ## Prepare obsm
  obsm <- lapply(drl, function(dr){
    x <- SeuratObject::Embeddings(dr)
    if(rm_rownames) rownames(x) <- NULL
    if(rm_colnames) colnames(x) <- NULL
    if(all(dim(x)==c(1,1)) && all(is.na(x))) return(NULL)
    x
  })
  if(isTRUE(drop_null)) obsm <- Filter(Negate(is.null), obsm)
  ## Prepare varm
  varm <- lapply(drl, function(dr){
    x <- SeuratObject::Loadings(dr)
    if(rm_rownames) rownames(x) <- NULL
    if(rm_colnames) colnames(x) <- NULL
    if(all(dim(x)==c(1,1)) && all(is.na(x))) return(NULL)
    x
  })
  if(isTRUE(drop_null)) varm <- Filter(Negate(is.null), varm)
  ## Return 
  return(
    list(
      obsm = obsm,
      varm = varm
    )
  )
}
 

add_DimReduc <- function(obj,
                         drl){
  if(length(drl)>0){
    for(nm in names(drl)){
      obj[[nm]] <- drl[[nm]]
    }
  }
  return(obj)
}