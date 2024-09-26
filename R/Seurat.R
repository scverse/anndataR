# TODO: export this?
guess_seurat_layers <- function(adata) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  layers <- list()

  if (!is.null(adata$X)) {
    layers["data"] <- list(NULL)
  }

  for (layer_name in names(adata$layers)) {
    layers[[layer_name]] <- layer_name
  }

  layers
}
# TODO: export this?
guess_seurat_embeddings <- function(adata) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  embeddings <- list()

  for (embedding_name in names(adata$obsm)) {
    if (grepl("^X_", embedding_name)) {
      name <- gsub("^X_", "", embedding_name)
      out <-
        if (embedding_name == "X_pca") {
          list(key = "PC_", obsm = "X_pca", varm = "PCs")
        } else {
          list(key = paste0(name, "_"), obsm = embedding_name, varm = NULL)
        }

      embeddings[[name]] <- out
    }
  }

  embeddings
}

guess_seurat_graphs <- function(adata) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  graphs <- list()

  for (graph_name in names(adata$obsp)) {
    if (graph_name == "connectivities") {
      graphs[["nn"]] <- graph_name
    } else if (grepl("^connectivities_", graph_name)) {
      new_name <- gsub("^connectivities_", "", graph_name)
      graphs[[new_name]] <- graph_name
    }
  }

  ## TODO: rename 'neighbors' to 'nn'?

  graphs
}

#' Convert a Seurat object to an AnnData object
#'
#' `to_Seurat()` converts an AnnData object to a Seurat object.
#'
#' @param obj An AnnData object
#' @param assay_name Name of the assay to be created
#' @param layer_mapping A named list mapping layer names to the names of the layers in the AnnData object. If the layer
#'   name is `NULL`, the `X` slot will be used. In the named list, at least 'counts' or 'data' must be present.
#' @param embedding_mapping A named list mapping embedding names to the names of the embeddings in the AnnData object. Each
#'   embedding must be a list with keys 'key', 'obsm', and 'varm'. The 'key' is the prefix of the embedding names, 'obsm'
#'   is the name of the embedding in `obsm`, and 'varm' is the name of the loadings in `varm`. If 'varm' is `NULL`, no
#'   loadings will be added.
#' @param graph_mapping A named list mapping graph names to the names of the graphs in the AnnData object.
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
    layer_mapping = guess_seurat_layers(adata),
    embedding_mapping = guess_seurat_embeddings(adata),
    graph_mapping = guess_seurat_graphs(adata)) {
  # nolint end: object_name_linter
  requireNamespace("SeuratObject")

  stopifnot(inherits(adata, "AbstractAnnData"))

  # store obs and var names
  obs_names <- adata$obs_names[]
  var_names <- adata$var_names[]

  # check seurat layers
  if (is.null(names(layer_mapping))) {
    names(layer_mapping) <- layer_mapping
  }
  if (!"counts" %in% names(layer_mapping) && !"data" %in% names(layer_mapping)) {
    stop(paste0(
      "layer_mapping must contain at least an item named \"counts\" or \"data\". Found names: ",
      paste(names(layer_mapping), collapse = ", ")
    ))
  }

  # trackstatus: class=Seurat, feature=get_obs, status=done
  # trackstatus: class=Seurat, feature=get_X, status=done
  # trackstatus: class=Seurat, feature=get_layers, status=done
  obj <- SeuratObject::CreateSeuratObject(
    meta.data = adata$obs,
    assay = assay_name,
    counts = .to_seurat_get_matrix_by_key(adata, layer_mapping, "counts"),
    data = .to_seurat_get_matrix_by_key(adata, layer_mapping, "data")
  )

  # trackstatus: class=Seurat, feature=get_var, status=done
  assay <- obj@assays[[assay_name]]
  if (!is.null(var)) {
    assay <- SeuratObject::AddMetaData(assay, metadata = adata$var)
  }

  # make sure obs and var names are set properly
  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  colnames(obj) <- obs_names
  rownames(obj) <- var_names

  # copy other layers
  for (i in seq_along(layer_mapping)) {
    from <- layer_mapping[[i]]
    to <- names(layer_mapping)[[i]]
    if (to == "counts" || to == "data") {
      next
    }
    SeuratObject::LayerData(obj, assay = assay_name, layer = to) <- adata$layers[[from]]
  }

  # copy embeddings
  # trackstatus: class=Seurat, feature=get_obsm, status=done
  # trackstatus: class=Seurat, feature=get_varm, status=done
  if (!is.null(embedding_mapping)) {
    if (!is.list(embedding_mapping) || (length(embedding_mapping) > 0 && is.null(names(embedding_mapping)))) {
      stop("embedding_mapping must be a named list")
    }
    for (i in seq_along(embedding_mapping)) {
      embedding_name <- names(embedding_mapping)[[i]]
      embedding <- embedding_mapping[[i]]

      if (!is.list(embedding) ||
        is.null(names(embedding)) ||
        !all(names(embedding) %in% c("key", "obsm", "varm")) ||
        !all(c("key", "obsm", "varm") %in% names(embedding))) {
        stop("each embedding must be a list with keys 'key', 'obsm', and 'varm'")
      }
      dr <- .to_seurat_process_reduction(
        adata = adata,
        key = embedding$key,
        obsm_embedding = embedding$obsm,
        varm_loadings = embedding$varm,
        assay_name = assay_name
      )
      if (!is.null(dr)) {
        obj[[embedding_name]] <- dr
      }
    }
  }

  # trackstatus: class=Seurat, feature=get_obsp, status=done
  for (i in seq_along(graph_mapping)) {
    graph_name <- names(graph_mapping)[[i]]
    graph <- graph_mapping[[i]]
    if (!is.character(graph) || length(graph) != 1) {
      stop("each graph must be a character vector of length 1")
    }
    obsp <- adata$obsp[[graph]]
    if (!is.null(obsp)) {
      dimnames(obsp) <- list(obs_names, obs_names)
      obsp_gr <- Seurat::as.Graph(obsp)
      obj[[paste0(assay_name, "_", graph_name)]] <- obsp_gr
    }
  }

  # trackstatus: class=Seurat, feature=get_uns, status=wip
  # TODO: should we store everything that is not stored elsewhere (e.g. unused obsm, varm, obsp) in misc?
  obj@misc <- adata$uns

  # trackstatus: class=Seurat, feature=get_varp, status=missing
  # TODO: could store varp in misc?

  obj
}

.to_seurat_is_atomic_character <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x)
}


.to_seurat_get_matrix_by_key <- function(adata, mapping, key) {
  if (!key %in% names(mapping)) {
    return(NULL)
  }

  layer_name <- mapping[[key]]

  if (is.null(layer_name)) {
    return(Matrix::t(adata$X))
  }

  if (!.to_seurat_is_atomic_character(layer_name)) {
    stop("layer_name must be the name of one of the layers or NULL")
  }

  if (!layer_name %in% names(adata$layers)) {
    stop("layer_name must be the name of one of the layers or NULL")
  }

  return(Matrix::t(adata$layers[[layer_name]]))
}

.to_seurat_get_matrix <- function(adata, layer_name) {
  if (is.null(layer_name)) {
    return(Matrix::t(adata$X))
  }

  if (!.to_seurat_is_atomic_character(layer_name)) {
    stop("layer_name must be the name of one of the layers or NULL")
  }

  if (!layer_name %in% names(adata$layers)) {
    stop("layer_name must be the name of one of the layers or NULL")
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
