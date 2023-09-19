#' @rdname Seurat
#'
#' @title Convert Between AnnData and Seurat
#'
#' @description `to_Seurat()` converts an AnnData object to a Seurat object.
#'
#' @param obj An AnnData object
#'
#' @importFrom Matrix t
#'
#' @export
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(cell = 1:3),
#'   obs_names = letters[1:3],
#'   var = data.frame(gene = 1:5),
#'   var_names = letters[1:5]
#' )
#' to_Seurat(ad)
# TODO: Add parameters to choose which how X and layers are translated into counts, data and scaled.data
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
#' Only one assay can be converted at a time.
#'
#' @param seurat_obj An object inheriting from Seurat.
#' @param output_class Name of the AnnData class. Must be one of `"HDF5AnnData"` or `"InMemoryAnnData"`.
#' @param assay Assay to be converted. If NULL, `DefaultAssay()` is used.
#' @param X Which of 'counts', 'data', or 'scale.data' will be used for X. By default, 'counts' will be used (if it is
#'   not empty), followed by 'data', then 'scale.data'. The remaining non-empty slots will be stored in different
#'   layers.
#' @param ... Additional arguments passed to the generator function.
#'
#' @export
# TODO: add tests with Seurat objects not created by anndataR
from_Seurat <- function(seurat_obj, output_class = c("InMemoryAnnData", "HDF5AnnData"), assay = NULL, X = "counts", ...) { # nolint

  stopifnot(inherits(seurat_obj, "Seurat"))

  if (!is.null(X)) {
    if (!X %in% c("counts", "data", "scale.data")) {
      stop("X must be NULL or one of: 'counts', 'data', 'scale.data'")
    }
  }

  # If a specific assay is selected, use it
  if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) {
      stop("'assay' must be NULL or one of: ", paste0("'", names(seurat_obj@assays), "'", collapse = ", "))
    }
    assay_name <- assay
  } else {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)

    # If Seurat object contains multiple assays, notify user the Default one is used
    if (length(names(seurat_obj@assays)) > 1) {
      message(
        "There are ", length(names(seurat_obj@assays)), " assays in the Seurat object; using the default assay ('",
        assay_name, "'). You can use the `assay` parameter to select a specific assay."
      )
    }
  }

  # get obs_names
  # trackstatus: class=Seurat, feature=set_obs_names, status=done
  obs_names <- colnames(seurat_obj)

  # get obs
  # trackstatus: class=Seurat, feature=set_obs, status=done
  obs <- seurat_obj@meta.data
  rownames(obs) <- NULL

  # construct var_names
  # trackstatus: class=Seurat, feature=set_var_names, status=done
  var_names <- rownames(seurat_obj@assays[[assay_name]])

  # construct var
  # trackstatus: class=Seurat, feature=set_var, status=done
  var <- seurat_obj@assays[[assay_name]]@meta.features
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
    if (all(dim(SeuratObject::GetAssayData(seurat_obj, slot = X, assay = assay_name)) == 0)) {
      stop("The '", X, "' slot is empty.")
    }

    assay_data <- SeuratObject::GetAssayData(seurat_obj, slot = X, assay = assay_name)

    # Remove names
    dimnames(assay_data) <- list(NULL, NULL)
    ad$X <- Matrix::t(assay_data)
  } else {
    # Cannot compare other values with NULL
    X <- "none"
  }

  # Add the remaining non-empty slots as layers
  slots <- c("counts", "data", "scale.data")
  slots <- slots[slots != X]

  for (slot in slots) {
    if (!all(dim(SeuratObject::GetAssayData(seurat_obj, slot = slot)) == 0)) {
      assay_data <- SeuratObject::GetAssayData(seurat_obj, slot = slot, assay = assay_name)
      dimnames(assay_data) <- list(NULL, NULL)
      ad$layers[[slot]] <- Matrix::t(assay_data)
    }
  }
  
  # TODO: add checks
  # TODO: can anything be stored in varp?

  # Dimensionality reduction
  for (reduction in names(seurat_obj@reductions)){
    # Check if the dimreduc was calculated by the selected assay
    if (seurat_obj@reductions[[reduction]]@assay.used != assay_name){
      next
    }
    
    # cell embeddings to obsm
    ad$obsm[[reduction]] <- seurat_obj@reductions[[reduction]]@cell.embeddings
    
    # gene loadings to varm
    feature_loadings <- seurat_obj@reductions[[reduction]]@feature.loadings
    feature_loadings_projected <- seurat_obj@reductions[[reduction]]@feature.loadings.projected
    
    if (!all(dim(feature_loadings) == 0)){
      ad$varm[[reduction]] <- seurat_obj@reductions[[reduction]]@feature.loadings
    }
    
    if (!all(dim(feature_loadings_projected) == 0)){
      ad$varm[[paste0(reduction, "_projected")]] <- seurat_obj@reductions[[reduction]]@feature.loadings.projected
    }
    
    # uns
    stdev <- seurat_obj@reductions[[reduction]]@stdev
    misc <- seurat_obj@reductions[[reduction]]@misc
    if (length(stdev) > 0){
      ad$uns[[paste0(reduction, "_stdev")]] <- stdev
    }
    
    if (length(misc) > 0){
      ad$uns[[paste0(reduction, "_misc")]] <- misc
    }
    
    # TODO: seurat_obj@reductions[[reduction]]@jackstraw
  }
  
  # Graphs (inherits from dgCMatrix)
  for (graph in names(seurat_obj@graphs)){
    # Check if graphs was created from selected assay
    if (grepl(paste0("^", assay_name), graph)){
      # pairwise distances -> obsp
      ad$obsp[[graph]] <- seurat_obj@graphs[[graph]]
    }
  }
  
  # Neighbors
  for (neighbor in names(seurat_obj@neighbors)){
    # Check if neighbors was created from selected assay
    if (grepl(paste0("^", assay_name), neighbor)){
      
      # Unlike Graphs, this is not a pairwise distance
      ad$obsm[[neighbor]] <- seurat_obj@neighbors[[neighbor]]@nn.idx
      ad$obsm[[neighbor]] <- seurat_obj@neighbors[[neighbor]]@nn.dist
      
      # Other parameters
      ad$uns[[neighbor]] <- seurat_obj@neighbors[[neighbor]]@alg.info
      
      # TODO: @alg.idx is not considered right now
      # TODO: check if @cell.names can be different from all cells
      
    }
  }

  return(ad)
}
