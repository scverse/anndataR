#' Convert a Seurat object to an AnnData object
#'
#' `to_Seurat()` converts an AnnData object to a Seurat object.
#'
#' @param obj An AnnData object
#'
#' @importFrom Matrix t
#'
#' @noRd
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' to_Seurat(ad)
# nolint start: object_name_linter
to_Seurat <- function(
  adata,
  assay_name = "RNA",
  counts_layer = ".X",
  data_layer = NULL,
  reductions = list(
    list(
      name = "pca",
      key = "PC_",
      obsm_embedding = "X_pca",
      varm_loadings = "PCs"
    ),
    list(
      name = "umap",
      key = "umap_",
      obsm_embedding = "X_umap"
    )
  )
) {
  # nolint end: object_name_linter
  requireNamespace("SeuratObject")

  stopifnot(inherits(adata, "AbstractAnnData"))

  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  obs_names <- adata$obs_names[]
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  var_names <- adata$var_names[]

  obj <- SeuratObject::CreateSeuratObject(
    # trackstatus: class=Seurat, feature=get_X, status=done
    # trackstatus: class=Seurat, feature=get_layers, status=done
    counts = .to_seurat_get_matrix(adata, counts_layer),
    data = .to_seurat_get_matrix(adata, data_layer),
    # trackstatus: class=Seurat, feature=get_obs, status=done
    meta.data = adata$obs,
    assay = assay_name
  )

  # trackstatus: class=Seurat, feature=get_var, status=done
  assay <- obj@assays[[assay_name]]
  if (!is.null(var)) {
    assay <- SeuratObject::AddMetaData(assay, metadata = var)
  }

  colnames(obj) <- obs_names
  rownames(obj) <- var_names

  # trackstatus: class=Seurat, feature=get_obsm, status=done
  for (reduction in reductions) {
    if (!is.list(reduction)) {
      stop("reductions must be a list of lists")
    }
    if (!all(c("name", "key", "obsm_embedding") %in% names(reduction))) {
      stop("reduction must contain 'name', 'key', and 'obsm_embedding'")
    }
    if (!.to_seurat_is_atomic_character(reduction$name)) {
      stop("reduction$name must be a character scalar")
    }
    dr <- .to_seurat_process_reduction(
      adata = adata,
      key = reduction$key,
      obsm_embedding = reduction$obsm_embedding,
      varm_loadings = reduction$varm_loadings,
      assay_name = assay_name
    )
    if (!is.null(dr)) {
      obj[[reduction$name]] <- dr
    }
  }
  # trackstatus: class=Seurat, feature=get_varm, status=wip
  # trackstatus: class=Seurat, feature=get_obsp, status=wip
  # trackstatus: class=Seurat, feature=get_varp, status=wip
  # trackstatus: class=Seurat, feature=get_uns, status=wip

  obj
}

.to_seurat_is_atomic_character <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x)
}

.to_seurat_get_matrix <- function(adata, layer_name) {
  if (is.null(layer_name)) {
    return(NULL)
  }

  if (!.to_seurat_is_atomic_character(layer_name)) {
    stop("layer_name must be the name of one of the layers, \".X\", or NULL")
  }

  if (layer_name == ".X") {
    return(Matrix::t(adata$X))
  }

  if (! layer_name %in% names(adata$layers)) {
    stop("layer_name must be the name of one of the layers, \".X\", or NULL")
  }

  return(Matrix::t(adata$layers[[layer_name]]))
}

.to_seurat_process_reduction <- function(adata, assay_name, key, obsm_embedding, varm_loadings) {
  requireNamespace("SeuratObject")
  if (!.to_seurat_is_atomic_character(key)) {
    stop("key must be a character scalar")
  }
  if (!.to_seurat_is_atomic_character(obsm_embedding)) {
    stop("obsm_embedding must be a character scalar")
  }
  if (!is.null(varm_loadings) && !.to_seurat_is_atomic_character(varm_loadings)) {
    stop("varm_loadings must be a character scalar or NULL")
  }
  embed <- adata$obsm[[obsm_embedding]]

  if (is.null(embed)) {
    warning(paste0("The embedding ", obsm_embedding, " is not present in adata$obsm, skipping"))
    return(NULL)
  }

  rownames(embed) <- adata$obs_names

  loadings <-
    if (is.null(varm_loadings)) {
      NULL
    } else if (!varm_loadings %in% names(adata$varm)) {
      warning(paste0("The loadings ", varm_loadings, " is not present in adata$varm, skipping"))
      NULL
    } else {
      adata$varm[[varm_loadings]]
    }

  args <- list(
    embeddings = embed,
    key = key,
    assay = assay_name,
    global = TRUE
  )

  if (!is.null(loadings)) {
    rownames(loadings) <- adata$var_names
    args$loadings <- loadings
  }

  do.call(SeuratObject::CreateDimReducObject, args)
}


#' Convert a Seurat object to an AnnData object
#'
#' `from_Seurat()` converts a Seurat object to an AnnData object.
#' Only one assay can be converted at a time.
#'
#' For more information on the functionality of an AnnData object, see [anndataR-package].
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
#'
#' @seealso [anndataR-package]
# TODO: Add examples
# nolint start: object_name_linter
from_Seurat <- function(
    # nolint end: object_name_linter
    seurat_obj,
    output_class = c("InMemoryAnnData", "HDF5AnnData"),
    assay = NULL,
    X = "counts",
    ...) {
  output_class <- match.arg(output_class)

  stopifnot(inherits(seurat_obj, "Seurat"))

  # if (!is.null(X)) {
  #   if (!X %in% c("counts", "data", "scale.data")) {
  #     stop("X must be NULL or one of: 'counts', 'data', 'scale.data'")
  #   }
  # }

  # # If a specific assay is selected, use it
  # if (!is.null(assay)) {
  #   if (!assay %in% names(seurat_obj@assays)) {
  #     stop("'assay' must be NULL or one of: ", paste0("'", names(seurat_obj@assays), "'", collapse = ", "))
  #   }
  #   assay_name <- assay
  # } else {
  #   assay_name <- SeuratObject::DefaultAssay(seurat_obj)

  #   # If Seurat object contains multiple assays, notify user the Default one is used
  #   if (length(names(seurat_obj@assays)) > 1) {
  #     message(
  #       "There are ", length(names(seurat_obj@assays)), " assays in the Seurat object; using the default assay ('",
  #       assay_name, "'). You can use the `assay` parameter to select a specific assay."
  #     )
  #   }
  # }

  # # get obs
  # # trackstatus: class=Seurat, feature=set_obs_names, status=done
  # # trackstatus: class=Seurat, feature=set_obs, status=done
  # obs <- seurat_obj@meta.data
  # rownames(obs) <- colnames(seurat_obj) # TODO: this is probably not needed

  # # construct var
  # # trackstatus: class=Seurat, feature=set_var_names, status=done
  # # trackstatus: class=Seurat, feature=set_var, status=done
  # var <- seurat_obj@assays[[assay_name]]@meta.features
  # rownames(var) <- rownames(seurat_obj@assays[[assay_name]]) # TODO: this is probably not needed

  # # use generator to create new AnnData object
  # generator <- get_anndata_constructor(output_class)
  # ad <- generator$new(
  #   obs = obs,
  #   var = var,
  #   ...
  # )

  # if (!is.null(X)) {
  #   # Check if the slot is not empty
  #   if (all(dim(SeuratObject::GetAssayData(seurat_obj, slot = X, assay = assay_name)) == 0)) {
  #     stop("The '", X, "' slot is empty.")
  #   }

  #   assay_data <- SeuratObject::GetAssayData(seurat_obj, slot = X, assay = assay_name)

  #   # Remove names
  #   dimnames(assay_data) <- list(NULL, NULL)
  #   ad$X <- Matrix::t(assay_data)
  # } else {
  #   # Cannot compare other values with NULL
  #   X <- "none"
  # }

  # # Add the remaining non-empty slots as layers
  # slots <- c("counts", "data", "scale.data")
  # slots <- slots[slots != X]

  # for (slot in slots) {
  #   if (!all(dim(SeuratObject::GetAssayData(seurat_obj, slot = slot)) == 0)) {
  #     assay_data <- SeuratObject::GetAssayData(seurat_obj, slot = slot, assay = assay_name)
  #     dimnames(assay_data) <- list(NULL, NULL)
  #     ad$layers[[slot]] <- Matrix::t(assay_data)
  #   }
  # }

  # return(ad)
}
