#' @rdname Seurat
#'
#' @title Convert between AnnData and Seurat
#'
#' @description `to_Seurat()` converts an AnnData object to a Seurat object.
#'
#' @param obj An AnnData object
#' @inheritParams to_DimReduc
#' @inheritDotParams SeuratObject::CreateAssayObject
#' @returns A \link[SeuratObject]{SeuratObject}.
#'
#' @importFrom Matrix t
#'
#' @export
#' @examples
#' adata <- dummy_data("InMemoryAnnData")
#' to_Seurat(adata)
# TODO: Add parameters to choose which how X and layers are translated into 
# counts, data and scaled.data
to_Seurat <- function(obj,
                      key_map = list("X_pca" = "PCs", 
                                     "X_umap" = NULL,
                                     "X_mofa" = "LFs"),
                      ...) { # nolint
  requireNamespace("SeuratObject")
  
  stopifnot(inherits(obj, "AbstractAnnData"))
  
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
  # TODO: should x_ be passed to counts or to data?
  # TODO: creating a seurat object when th AnnData doesn't contain X or layers
  # probably doesn't make any sense
  x_ <-
    if (!is.null(obj$X)) {
      Matrix::t(obj$X)
    } else {
      mat <- Matrix::sparseMatrix(
        i = integer(0),
        p = c(0L),
        x = integer(0),
        dims = c(obj$n_var(), obj$n_obs())
      )
      attr(mat, "is_X_null") <- TRUE # nolint
      mat
    }
  dimnames(x_) <- list(var_names_, obs_names_)
  x_assay <- SeuratObject::CreateAssayObject(counts = x_, 
                                             ...)
  
  # create seurat object
  if (ncol(var_) > 0) {
    # don't add var metadata if the data frame does not contain any columns
    x_assay <- SeuratObject::AddMetaData(x_assay,
                                         metadata = var_)
  }
  seur <- SeuratObject::CreateSeuratObject(x_assay, 
                                          meta.data = obs_)
  
  ## add layers
  # trackstatus: class=Seurat, feature=get_layers, status=wip
  # TODO: should values be passed to counts or to data?
  for (key in obj$layers_keys()) {
    layer_ <- Matrix::t(obj$layers[[key]])
    dimnames(layer_) <- list(var_names_, obs_names_)
    seur[[key]] <- SeuratObject::CreateAssayObject(counts = layer_)
  } 
  ## Add DimReduc objects
  drl <- to_DimReduc(obj = obj)
  seur <- add_DimReduc(obj = seur, 
                       drl = drl)
  ## Add obsp/varp as graphs 
  seur <- anndata_to_graphs(ad = obj,
                            seur = seur) 
  ## Return 
  return(seur)
}

.toseurat_check_obsvar_names <- function(names, label) {
  if (any(grepl("_", names))) {
    ## mimic seurat behaviour
    wrn <- wrap_message(
      shQuote(label),
      " cannot have underscores ('_') when converting to Seurat;",
      " replacing with dashes ('-')."
    )
    warning(wrn)
    names <- gsub("_", "-", names)
  }
  names
}

#' @rdname Seurat
#'
#' @title Convert between AnnData and Seurat
#' 
#' @description `from_Seurat()` converts a Seurat object to an AnnData object.
#' Only one assay can be converted at a time.
#'
#' @param obj An object inheriting from \link[Seurat]{Seurat}.
#' @param output_class Name of the AnnData class. 
#' Must be one of:
#' \itemize{
#' \item{"HDF5AnnData"}{For \link[anndataR]{HDF5AnnData} output.}
#' \item{"InMemoryAnnData"}{For \link[anndataR]{InMemoryAnnData} output.}
#' }
#' @param assay Assay to be converted. If \code{NULL}, 
#' \link[SeuratObject]{DefaultAssay} is used.
#' @param X Which of 'counts', 'data', or 'scale.data' will be used for X. 
#' By default, 'counts' will be used (if it is
#'   not empty), followed by 'data', then 'scale.data'. 
#'   The remaining non-empty slots will be stored in different
#'   layers.
#' @param update Run \link[SeuratObject]{UpdateSeuratObject} on the 
#' input \code{obj} beforehand.
#' @param ... Additional arguments passed to the generator function,
#' which will be one of the following depending on the selected 
#' \code{output_class}:
#' \itemize{
#' \item{"HDF5AnnData"}{For \link[anndataR]{HDF5AnnData} output.}
#' \item{"InMemoryAnnData"}{For \link[anndataR]{InMemoryAnnData} output.}
#' }
#' @returns An \link[anndataR]{AnnData} object.
#'
#' @export
#' @examples
#' obj <- dummy_data("Seurat")
#' adata <- from_Seurat(obj)
#' adata
from_Seurat <- function(obj, 
                        output_class = c("InMemoryAnnData", "HDF5AnnData"), 
                        assay = NULL, 
                        X = "counts", 
                        update = FALSE,
                        ...) { # nolint
  # devoptera::args2vars(from_Seurat)
  
  stopifnot(inherits(obj, "Seurat"))
  
  if(isTRUE(update)) obj <- SeuratObject::UpdateSeuratObject(obj)
  output_class <- output_class[1]
  
  if (!is.null(X)) {
    if (!X %in% c("counts", "data", "scale.data")) {
      stop("X must be NULL or one of: 'counts', 'data', 'scale.data'")
    }
  }
  
  # If a specific assay is selected, use it
  if (!is.null(assay)) {
    if (!assay %in% names(obj@assays)) {
      stp <- paste(
        "'assay' must be NULL or one of:", 
        paste0("'", names(obj@assays), "'", collapse = ", ")
      )
      stop(stp)
    }
    assay_name <- assay
  } else {
    assay_name <- SeuratObject::DefaultAssay(obj)
    
    # If Seurat object contains multiple assays,
    # notify user the Default one is used.
    if (length(names(obj@assays)) > 1) {
      messager(
        "There are ", length(names(obj@assays)), 
        "assays in the Seurat object; using the default assay ",
        paste0("(",shQuote(assay_name),")."),
        "You can use the `assay` parameter to select a specific assay."
      )
    }
  }
  
  # get obs_names
  # trackstatus: class=Seurat, feature=set_obs_names, status=done
  obs_names <- colnames(obj)
  
  # get obs
  # trackstatus: class=Seurat, feature=set_obs, status=done
  obs <- obj@meta.data
  rownames(obs) <- NULL
  
  # construct var_names
  # trackstatus: class=Seurat, feature=set_var_names, status=done
  var_names <- rownames(obj@assays[[assay_name]])
  
  # construct var
  # trackstatus: class=Seurat, feature=set_var, status=done
  assay_slots <- methods::slotNames(obj@assays[[assay_name]])
  slot_opts <- c("meta.features","meta.data")
  slot_opts <- slot_opts[slot_opts %in% assay_slots]
  if(length(assay_slots)>0){
    var <- methods::slot(obj@assays[[assay_name]],
                         name = slot_opts[1])
  } 
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
  
  if (!is.null(X)) {
    # Check if the slot is not empty
    if (all(dim(SeuratObject::GetAssayData(obj, 
                                           layer = X, 
                                           assay = assay_name)) == 0)) {
      stop("The ",shQuote(X)," slot is empty.")
    } 
    assay_data <- SeuratObject::GetAssayData(obj, 
                                             layer = X, 
                                             assay = assay_name) 
    # Remove names
    dimnames(assay_data) <- list(NULL, NULL)
    ad$X <- Matrix::t(assay_data)
  } else {
    # Cannot compare other values with NULL
    X <- "none"
  }
  
  # Add the remaining non-empty "slots" (renamed "layers" in Seurat>=5.0).
  layer_opts <- names(obj@assays[[assay_name]]@layers)
  layers <- c("counts", "data", "scale.data")
  layers <- layers[layers != X & layers %in% layer_opts]
  for (layer in layers) { 
    if (!all(dim(SeuratObject::GetAssayData(obj, layer = layer)) == 0)) {
      assay_data <- SeuratObject::GetAssayData(obj,
                                               layer = layer,
                                               assay = assay_name)
      dimnames(assay_data) <- list(NULL, NULL)
      ad$layers[[layer]] <- Matrix::t(assay_data)
    }
  }
  ## Add obsm and varm
  obsm_varm <- from_DimReduc(obj)
  ad$obsm <- obsm_varm$obsm
  ad$varm <- obsm_varm$varm  
  ## Add obsp/varp from graphs 
  ad <- graphs_to_anndata(seur = obj, 
                          ad = ad) 
  ## Return 
  return(ad)
}


#' Graphs to AnnData
#' 
#' Transfer Seurat graphs to relevant AnnData slots.
#' @param ad An AnnData object.
#' @param seur A Seurat object.
#' @returns An AnnData object. 
#' @keywords internal 
graphs_to_anndata <- function(seur,
                              ad){
  graph_nms <- names(seur@graphs)
  if("obsp" %in% graph_nms){
    ad$obsp <- seur@graphs$obsp 
  }
  if("varp" %in% graph_nms){
    ad$varp <- seur@graphs$varp 
  } 
  extra_graphs <- graph_nms[!graph_nms %in% c("obsp","varp")]
  for(gn in extra_graphs){
    g <- seur@graphs[[gn]] 
    if(all(dim(g)==nrow(seur))){
      messager("Inferring graph",shQuote(gn),"as a obs x obs matrix.",
               "Adding to 'obsp' slot.")
      ad$obsp[[gn]] <- g
    } else if (all(dim(g)==ncol(seur))){
      messager("Inferring graph",shQuote(gn),"as a var x var matrix.",
               "Adding to 'varp' slot.")
      ad$varp[[gn]] <- g
    } else {
      messager("Unable to infer where graph",shQuote(gn),"should be placed.")
    }
  }
  return(ad)
}

#' Graphs to AnnData
#' 
#' Transfer Seurat graphs to relevant AnnData slots.
#' @returns A \link[SeuratObject]{SeuratObject}.
#' @inheritParams graphs_to_anndata
#' @keywords internal 
anndata_to_graphs <- function(ad,
                              seur){
  seur@graphs <- list(ad$obsp,
                      ad$varp)
  return(seur)
}