#' @rdname Seurat
#'
#' @title Convert between AnnData and Seurat
#'
#' @description `to_Seurat()` converts an AnnData object to a Seurat object.
#'
#' @param obj An AnnData object
#' @inheritParams .to_Seurat_DimReduc
#' @inheritDotParams SeuratObject::CreateAssayObject
#' @returns A \link[SeuratObject]{SeuratObject}.
#'
#' @importFrom Matrix t
#'
#' @export
#' @examples
#' adata <- generate_dataset(format="AnnData")
#' to_Seurat(adata)
# TODO: Add parameters to choose which how X and layers are translated into
# counts, data and scaled.data
to_Seurat <- function(obj,
                      key_map = list("X_pca" = "PCs",
                                     "X_umap" = NULL,
                                     "X_mofa" = "LFs"),
                      ...) { # nolint
  # devoptera::args2vars(to_Seurat)
  requireNamespace("SeuratObject")
  stopifnot(inherits(obj, "AbstractAnnData"))

  # translate var_names
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  var_names_ <- .to_Seurat_check_obsvar_names(obj$var_names, "var_names")
  # translate obs_names
  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  obs_names_ <- .to_Seurat_check_obsvar_names(obj$obs_names, "obs_names")
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
    key <- .to_Seurat_check_layer_names(key)
    dimnames(layer_) <- list(var_names_, obs_names_)
    seur[[key]] <- SeuratObject::CreateAssayObject(counts = layer_,
                                                   key = key)
  }
  ## Add DimReduc objects
  drl <- .to_Seurat_DimReduc(obj = obj,
                             ...)
  seur <- .add_Seurat_DimReduc(obj = seur,
                               drl = drl)
  ## Return
  return(seur)
}

.to_Seurat_check_obsvar_names <- function(names, label) {
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

.to_Seurat_check_obsvar_names <- function(names, label) {
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

.to_Seurat_check_layer_names <- function(key,
                                         pattern="_|-|[.]",
                                         replacement="",
                                         suffix="_"){
  paste0(gsub(pattern,replacement,key),suffix)
}


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
#' @keywords internal
.to_Seurat_DimReduc <- function(obj = NULL,
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
    message("obsm is NULL; returning NULL.")
    return(NULL)
  }
  if(!is.list(obsm) || is.null(names(obsm))){
    message("obsm name not provided; setting name to 'obsm'.")
    obsm <- list(obsm=obsm)
  }
  if(!is.list(varm) || is.null(names(varm))){
    message("varm name not provided; setting name to 'varm'.")
    varm <- list(varm=varm)
  }

  .nms_msg <- function(X,
                       X_name,
                       X_type="loadings",
                       dim=2){
    dim_nm <- if(dim==1) "rownames" else "colnames"
    message(paste(
      "No",dim_nm,"present in",paste0(X_type,"; setting to"),
      shQuote(paste0(X_name,"[1:",ncol(X),"]")
      ))
    )
  }
  # Add embeddings
  lapply(stats::setNames(names(obsm),
                         names(obsm)),
         function(key){
           message(paste("Preparing DimReduc for:",shQuote(key)))
           ## Prepare obsm
           obsm_nm2 <- .to_Seurat_check_layer_names(gsub('X_', '', key))
           if(!key %in% names(key_map)){
             key_map[[key]] <- key
           }
           embeddings <- obsm[[key]] |> as.matrix()
           ### Precheck embedding is indeed numeric
           if(!is.numeric(embeddings[,1])){
             message(paste(
               "obsm",shQuote(key),"is not numeric; returning NULL."
             ))
             return(NULL)
           }
           if(is.null(rownames(embeddings))){
             if(is.null(obs_names)) stop("Must provide `obs_names`.")
             rownames(embeddings) <- obs_names
           }
           if(is.null(colnames(embeddings)) ||
              !all(grepl(paste0(obsm_nm2,"_"),colnames(embeddings)))){
             .nms_msg(X = embeddings,
                      X_name = obsm_nm2,
                      X_type = "embeddings")
             colnames(embeddings) <- paste(obsm_nm2,
                                           seq(ncol(embeddings)),sep="")
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
                 .nms_msg(X = embeddings,
                          X_name = obsm_nm2,
                          X_type = "embeddings")
                 colnames(loadings) <- paste(obsm_nm2,
                                             seq(ncol(loadings)),sep="")
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
             key = obsm_nm2,
             stdev = emb_stdev,
             ...
           )
         })
}

.add_Seurat_DimReduc <- function(obj,
                                 drl){
  if(length(drl)>0){
    for(nm in names(drl)){
      obj[[nm]] <- drl[[nm]]
    }
  }
  return(obj)
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
#' @returns \link[anndataR]{InMemoryAnnData} or \link[anndataR]{HDF5AnnData}
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
  generator <- get_anndata_constructor(output_class)
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
  ## Add obsm and varm
  obsm_varm <- .from_Seurat_DimReduc(seurat_obj)
  ad$obsm <- obsm_varm$obsm
  ad$varm <- obsm_varm$varm
  return(ad)
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
#' @keywords internal
.from_Seurat_DimReduc <- function(obj,
                                  rm_rownames=TRUE,
                                  rm_colnames=TRUE,
                                  drop_null=TRUE
){

  requireNamespace("SeuratObject")

  ## Handle multiple object types
  if(methods::is(obj,"DimReduc")){
    stop("Must provide DimReduc objects as a named list.")
  } else if(SeuratObject::IsNamedList(obj)){
    if(!all( unlist(lapply(obj,methods::is,"DimReduc")) ) ){
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

