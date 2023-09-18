#' @rdname Seurat
#'
#' @title Convert Between AnnData and Seurat
#'
#' @description `to_Seurat()` converts an AnnData object to a Seurat object.
#' NOTE: Only 1 assay per object is currently. 
#' @param obj An AnnData object.
#'
#' @importFrom Matrix t
#'
#' @export
#' @examples
#' ad <- example_anndata()
#' to_Seurat(ad) 
to_Seurat <- function(obj) { # nolint
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
        dims = c(obj$n_vars(), obj$n_obs())
      )
      attr(mat, "is_X_null") <- TRUE # nolint
      mat
    }
  dimnames(x_) <- list(var_names_, obs_names_)
  x_assay <- SeuratObject::CreateAssayObject(counts = x_)
  
  # create seurat object
  if (ncol(var_) > 0) {
    # don't add var metadata if the data frame does not contain any columns
    x_assay <- SeuratObject::AddMetaData(x_assay, metadata = var_)
  }
  seurat_obj <- SeuratObject::CreateSeuratObject(x_assay, meta.data = obs_)
  
  # add layers
  # trackstatus: class=Seurat, feature=get_layers, status=wip
  # TODO: should values be passed to counts or to data?
  for (key in obj$layers_keys()) {
    layer_ <- t(obj$layers[[key]])
    dimnames(layer_) <- list(var_names_, obs_names_)
    seurat_obj[[key]] <- SeuratObject::CreateAssayObject(counts = layer_)
  }
  
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
#' @inheritParams example_data
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
