#' Convert an AnnData object to Seurat
#'
#' @param obj An AnnData object
#'
#' @importFrom Matrix t
#'
#' @export
#' @examples
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(cell = 1:3),
#'   obs_names = letters[1:3],
#'   var = data.frame(gene = 1:5),
#'   var_names = letters[1:5]
#' )
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
  obs_ <- obj$obs
  rownames(obs_) <- obs_names_

  # translate X
  # trackstatus: class=Seurat, feature=get_X, status=wip
  # TODO: what if there is no X?
  # TODO: should x_ be passed to counts or to data?
  x_ <- Matrix::t(obj$X)
  dimnames(x_) <- list(var_names_, obs_names_)
  x_assay <- SeuratObject::CreateAssayObject(counts = x_)

  # create seurat object
  x_assay <- SeuratObject::AddMetaData(x_assay, metadata = var_)
  seurat_obj <- SeuratObject::CreateSeuratObject(x_assay, meta.data = obs_)

  # add layers
  # trackstatus: class=Seurat, feature=get_layers, status=wip
  # TODO: should x_ be passed to counts or to data?
  for (key in obj$layers_keys()) {
    layer_ <- t(obj$layers[[key]])
    dimnames(layer_) <- list(var_names_, obs_names_)
    seurat_obj[[key]] <- Seurat::CreateAssayObject(counts = layer_)
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

# todo: add option to determine default X assay name
# todo: add option to decide on var key resolution
# todo: add option to decide on whether to use counts or data for the different layers
# todo: add tests with Seurat objects not created by anndataR
# todo: add option to decide on which backend to use ( https://github.com/scverse/anndataR/issues/51 )
from_Seurat <- function(seurat_obj) { # nolint
  # get obs_names
  # trackstatus: class=Seurat, feature=set_obs_names, status=done
  obs_names <- colnames(seurat_obj)

  # get obs
  # trackstatus: class=Seurat, feature=set_obs, status=done
  obs <- seurat_obj@meta.data
  rownames(obs) <- NULL

  # construct var_names
  # trackstatus: class=Seurat, feature=set_var_names, status=done
  var_names <- unique(unlist(lapply(seurat_obj@assays, rownames)))

  # detect conflicting var keys
  var_keys_checker <- do.call(rbind, lapply(names(seurat_obj@assays), function(assay_name) {
    var_keys <- colnames(seurat_obj@assays[[assay_name]]@meta.features)
    if (length(var_keys) > 0) {
      data.frame(
        assay = assay_name,
        var_key = var_keys
      )
    } else {
      NULL
    }
  }))
  duplicated_var_keys <- unique(var_keys_checker$var_key[duplicated(var_keys_checker$var_key)])
  var_keys_checker$conflict <- var_keys_checker$var_key %in% duplicated_var_keys
  var_keys_checker$renamed_var_key <- ifelse(
    var_keys_checker$conflict,
    paste0(var_keys_checker$assay, "_", var_keys_checker$var_key),
    var_keys_checker$var_key
  )

  # construct var
  # NOTE: This will create duplicate data if the meta features are the same.
  # trackstatus: class=Seurat, feature=set_var, status=done
  var <- data.frame(
    row.names = var_names
  )
  for (assay_name in names(seurat_obj@assays)) {
    meta_feat <- seurat_obj@assays[[assay_name]]@meta.features
    var_keys_assay <- var_keys_checker[var_keys_checker$assay == assay_name, , drop = FALSE]

    if (nrow(var_keys_assay) > 0) {
      var[rownames(meta_feat), var_keys_assay$renamed_var_key] <- meta_feat[, var_keys_assay$var_key, drop = FALSE]
    }
  }
  rownames(var) <- NULL

  # construct X
  # TODO: Maybe we shouldn't use counts but instead data
  # TODO: check whether X has all of the obs_names and var_names
  # trackstatus: class=Seurat, feature=set_X, status=wip
  X <- Matrix::t(SeuratObject::GetAssayData(seurat_obj, "counts", assay = "RNA"))
  dimnames(X) <- list(NULL, NULL)

  # create AnnData
  ad <- InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names
  )

  for (assay_name in names(seurat_obj@assays)) {
    # TODO: Maybe we shouldn't use counts but instead data
    # TODO: check whether X has all of the obs_names and var_names
    # trackstatus: class=Seurat, feature=set_layers, status=wip
    layer <- Matrix::t(SeuratObject::GetAssayData(seurat_obj, "counts", assay = "RNA"))
    dimnames(layer) <- list(NULL, NULL)
    ad$layers[[assay_name]] <- layer
  }

  return(ad)
}
