#' @rdname Seurat
#'
#' @title Convert Between AnnData and Seurat
#'
#' @description `to_Seurat()` converts an AnnData object to a Seurat object.
#' NOTE: Only 1 assay per object is currently. 
#' @param obj An AnnData object.
#' @inheritDotParams SeuratObject::CreateSeuratObject
#'
#' @importFrom Matrix t
#'
#' @export
#' @examples
#' obj <- dummy_data(output="InMemoryAnnData")
#' to_Seurat(obj) 
to_Seurat <- function(obj, 
                      use_layer=NULL,
                      ...) { # nolint
  requireNamespace("SeuratObject")
  
  stopifnot(inherits(obj, "AbstractAnnData"))
  use_layer <- use_layer[1]
  # translate var_names
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  var_names_ <- .toseurat_check_obsvar_names(obj$var_names, "var_names")
  
  # translate obs_names
  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  obs_names_ <- .toseurat_check_obsvar_names(obj$obs_names, "obs_names")
  
  # translate var
  # trackstatus: class=Seurat, feature=get_var, status=done
  var_ <- obj$var
  rownames(var_) <- var_names_
  
  # translate obs
  # trackstatus: class=Seurat, feature=get_obs, status=done
  obs_ <-
    if (ncol(obj$obs) > 0) {
      ob <- obj$obs
      rownames(ob) <- obs_names_
      ob
    } else {
      NULL
    }
  
  # translate X
  # trackstatus: class=Seurat, feature=get_X, status=wip
  seurat_slots <- c("counts","data","scale.data")
  x_ <-
    if(!is.null(use_layer) && 
       use_layer %in% obj$layers_keys()){
      seurat_slots <- seurat_slots[seurat_slots!=use_layer]
      Matrix::t(obj$layers[[use_layer]])
    } else if (is.null(obj$X) &&
               !is.null(obj$layers_keys())){
      layer1 <- obj$layers_keys()[1]
      seurat_slots <- seurat_slots[seurat_slots!=layer1]
      message("No X detected. Using first layer as counts: ",layer1)
      Matrix::t(obj$layers[[layer1]])
    } else if (!is.null(obj$X)) {
      Matrix::t(obj$X)
    } else {
      stop("obj must contain at least one matrix in either `X` or `layers`.")
    }
  dimnames(x_) <- list(var_names_, obs_names_) 
    x_assay <- SeuratObject::CreateAssayObject(counts = x_)
  # seurat <- SeuratObject::SetAssayData(seurat, "scale.data", X3)
  
  # create seurat object
  if (ncol(var_) > 0) {
    # don't add var metadata if the data frame does not contain any columns
    x_assay <- SeuratObject::AddMetaData(x_assay, metadata = var_)
  }
  seurat_obj <- SeuratObject::CreateSeuratObject(x_assay, 
                                                 meta.data = obs_, 
                                                 ...)
  
  # add layers
  # trackstatus: class=Seurat, feature=get_layers, status=wip
  #### Store as assay data matrices ####
  if(all(seurat_slots %in% obj$layers_keys())){
    message("Setting layers as AssayData matrices.")
    for(key in seurat_slots){
      seurat_obj <- SeuratObject::SetAssayData(object = seurat_obj, 
                                               slot = key,
                                               new.data = t(obj$layers[[key]])) 
    }
  #### Store as multiple assays ####
  } else {
    message("Setting layers as new assays.")
    for (key in obj$layers_keys()) {
      layer_ <- t(obj$layers[[key]])
      dimnames(layer_) <- list(var_names_, obs_names_)
      seurat_obj[[key]] <- SeuratObject::CreateAssayObject(counts = layer_)
    }
  } 
  #### Add DimReducObjects ####
  if(all(c("obsm_keys","varm_keys") %in% names(obj)) &&
     !is.null(obj$obsm_keys()) &&
     !is.null(obj$varm_keys()) ){  
    for(key in obj$obsm_keys()){
      assay <- SeuratObject::Assays(seurat_obj)[1]
      seurat_obj@reductions[[key]] <- to_DimReduc(obsm=obj$obsm[[key]],
                                                  varm=obj$varm[[key]],
                                                  stdev=obj$stdev[[key]],
                                                  assay=assay,
                                                  key=key)
    } 
  }  
  #### Return ####
  seurat_obj
}

.toseurat_check_obsvar_names <- function(names, label) {
  if (any(grepl("_", names))) {
    # mimic seurat behaviour
    warning(wrap_message(
      "'", label, "' ",
      "cannot have underscores ('_') when converting to Seurat, ",
      "replacing with dashes ('-')"
    ))
    names <- gsub("_", "-", names)
  }
  
  names
}

#' @rdname Seurat
#'
#' @description `from_Seurat()` converts a Seurat object to an AnnData object.
#'
#' @param seurat_obj An object inheriting from Seurat.
#'
#' @inheritParams dummy_anndata
#' 
#' @returns \link[anndataR]{InMemoryAnnData} or \link[anndataR]{HDF5AnnData}
#'
#' @export
# TODO: Add parameter to choose which how counts, data and scaled.data are translated into X and layers
# TODO: add tests with Seurat objects not created by anndataR
from_Seurat <- function(seurat_obj, 
                        output_class = c("InMemoryAnnData",
                                         "HDF5AnnData"), 
                        ...) { # nolint
  
  stopifnot(inherits(seurat_obj, "Seurat"))
  
  # get obs_names
  # trackstatus: class=Seurat, feature=set_obs_names, status=done
  obs_names <- colnames(seurat_obj)
  
  # get obs
  # trackstatus: class=Seurat, feature=set_obs, status=done
  obs <- seurat_obj@meta.data
  rownames(obs) <- NULL
  
  # construct var_names
  # trackstatus: class=Seurat, feature=set_var_names, status=done
  var_names <- rownames(seurat_obj)
  
  # construct var
  # trackstatus: class=Seurat, feature=set_var, status=done
  var <- seurat_obj@assays[[seurat_obj@active.assay]]@meta.features
  rownames(var) <- NULL
  
  # use generator to create new AnnData object
  generator <- get_generator(output_class)
  ad <- generator$new(
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names,
    ...
  )
  
  # trackstatus: class=Seurat, feature=set_X, status=wip
  # trackstatus: class=Seurat, feature=set_layers, status=wip
  for (assay_name in names(seurat_obj@assays)) {
    # TODO: Maybe we shouldn't use counts but instead data
    assay_data <- SeuratObject::GetAssayData(seurat_obj, "counts", assay = assay_name)
    
    if (nrow(assay_data) != length(var_names) || !identical(rownames(assay_data), var_names)) {
      warning(
        "Skipping assay '", assay_name, "' because it has different feature names ",
        "than the active assay ('", seurat_obj@active.assay, "')."
      )
      next
    }
    if (ncol(assay_data) != length(obs_names) || !identical(colnames(assay_data), obs_names)) {
      warning(
        "Skipping assay '", assay_name, "' because it has different cell names ",
        "than the active assay ('", seurat_obj@active.assay, "')."
      )
      next
    }
    
    # remove names
    dimnames(assay_data) <- list(NULL, NULL)
    if (assay_name == seurat_obj@active.assay) {
      ad$X <- Matrix::t(assay_data)
    } else {
      ad$layers[[assay_name]] <- Matrix::t(assay_data)
    }
  }
  
  return(ad)
}
