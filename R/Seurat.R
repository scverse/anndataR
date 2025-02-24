#' @title Convert a Seurat object to an AnnData object
#'
#' @description
#' `to_Seurat()` converts an AnnData object to a Seurat object. Only one assay can be converted at a time.
#' Arguments are used to configure the conversion. If `NULL`, the functions `to_Seurat_guess_*` will be used to guess
#' the mapping.
#'
#' @param adata An AnnData object to be converted
#' @param assay_name Name of the assay to be created (default: "RNA").
#' @param layers_mapping A named list to map AnnData layers to Seurat layers. See section "Layer mapping" for more
#'   details.
#' @param reduction_mapping A named list to map AnnData reductions to Seurat reductions. Each item in the list must be a
#'   named list with keys 'key', 'obsm', and 'varm'. See section "Reduction mapping" for more details.
#' @param graph_mapping A named list to map AnnData graphs to Seurat graphs. Each item in the list must be a character
#'   vector of length 1. See section "Graph mapping" for more details.
#' @param misc_mapping A named list to map miscellaneous data to the names of the data in the Seurat object. See section
#'   "Miscellaneous mapping" for more details.
#'
#' @section Layer mapping:
#'
#' A named list to map AnnData layers to Seurat layers. Each item in the list must be a character vector of length 1,
#' where the values correspond to the names of the layers in the AnnData object, and the names correspond
#' to the names of the layers in the resulting Seurat object. A value of `NULL` corresponds to the AnnData `X` slot.
#'
#' Example: `layers_mapping = list(counts = "counts", data = NULL, foo = "bar")`.
#'
#' If `NULL`, the internal function `to_Seurat_guess_layers` will be used to guess the layer mapping as follows:
#'
#' * All AnnData layers are copied to Seurat layers by name.
#'
#' @section Reduction mapping:
#'
#' A named list to map AnnData `$obsm` and `$varm` to Seurat reductions. Each item in the list must be a named list
#' with keys `'key'`, `'obsm'`, and `'varm'`.
#
#' Example: `reduction_mapping = list(pca = list(key = "PC_", obsm = "X_pca", varm = "PCs"))`.
#'
#' If `NULL`, the internal function `to_Seurat_guess_reductions` will be used to guess the reduction mapping as follows:
#'
#' * All `$obsm` items starting with `X_` are copied by name.
#'
#' @section Graph mapping:
#'
#' A named list mapping graph names to the names of the graphs in the AnnData object. Each item in the list must be a
#' character vector of length 1. The values correspond to the names of the graphs in the resulting Seurat object, while
#' the names correspond to the names of the graphs in the AnnData object.
#'
#' Example: `graph_mapping = list(nn = "connectivities")`.
#'
#' If `NULL`, the internal function `to_Seurat_guess_graphs` will be used to guess the graph mapping as follows:
#'
#' * An obsp named `connectivities` will be mapped to `nn`.
#' * Other graphs starting with `connectivities_` are stripped of the prefix and copied by name.
#'
#' @section Miscellaneous mapping:
#'
#' A named list mapping miscellaneous data to the names of the data in the AnnData object. Each item in the list must be
#' a vector with one or two elements. The first element must be one of: 'X', 'layers', 'obs', 'obsm', 'obsp', 'var',
#' 'varm', 'varp', 'uns'. The second element is the name of the data in the corresponding slot. If the second element is
#' not present, the whole slot as specified by the first element will be used.
#'
#' Example: `misc_mapping = list(uns = "uns", varp_neighbors = c("varp", "neighbors"))`.
#'
#' If `NULL`, the internal function `to_Seurat_guess_misc` will be used to guess the miscellaneous mapping as follows:
#'
#' * If `$uns` is defined, all values in `$uns` are copied to the Seurat misc.
#'
#' @return A Seurat object
#'
#' @importFrom Matrix t
#'
#' @export
#'
#' @rdname to_Seurat
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
  layers_mapping = NULL,
  reduction_mapping = NULL,
  graph_mapping = NULL,
  misc_mapping = NULL
) {
  # nolint end: object_name_linter
  check_requires("Converting AnnData to Seurat", "SeuratObject")

  stopifnot(inherits(adata, "AbstractAnnData"))

  if (is.null(layers_mapping)) {
    layers_mapping <- to_Seurat_guess_layers(adata)
  }
  if (is.null(reduction_mapping)) {
    reduction_mapping <- to_Seurat_guess_reductions(adata)
  }
  if (is.null(graph_mapping)) {
    graph_mapping <- to_Seurat_guess_graphs(adata)
  }
  if (is.null(misc_mapping)) {
    misc_mapping <- to_Seurat_guess_misc(adata)
  }

  if (length(adata$layers) == 0 && is.null(adata$X)) {
    stop("to_Seurat: adata must have an $X slot or at least one layer")
  }

  # store obs and var names
  obs_names <- adata$obs_names[]
  var_names <- adata$var_names[]

  # check seurat layers
  if (is.null(names(layers_mapping))) {
    names(layers_mapping) <- layers_mapping
  }
  if (
    !"counts" %in% names(layers_mapping) && !"data" %in% names(layers_mapping)
  ) {
    stop(paste0(
      "layers_mapping must contain at least an item named \"counts\" or \"data\". Found names: ",
      paste(names(layers_mapping), collapse = ", ")
    ))
  }

  # trackstatus: class=Seurat, feature=get_obs, status=done
  # trackstatus: class=Seurat, feature=get_X, status=done
  # trackstatus: class=Seurat, feature=get_layers, status=done
  counts <- .to_seurat_get_matrix_by_key(adata, layers_mapping, "counts")
  data <- .to_seurat_get_matrix_by_key(adata, layers_mapping, "data")
  if (!is.null(counts)) {
    dimnames(counts) <- list(adata$var_names, adata$obs_names)
  }
  if (!is.null(data)) {
    dimnames(data) <- list(adata$var_names, adata$obs_names)
  }
  obj <- SeuratObject::CreateSeuratObject(
    meta.data = adata$obs,
    assay = assay_name,
    counts = counts,
    data = data
  )

  # trackstatus: class=Seurat, feature=get_var, status=done
  if (!is.null(adata$var)) {
    obj@assays[[assay_name]] <- SeuratObject::AddMetaData(
      obj@assays[[assay_name]],
      metadata = adata$var
    )
  }

  # make sure obs and var names are set properly
  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  colnames(obj) <- obs_names
  rownames(obj) <- var_names

  # copy other layers
  for (i in seq_along(layers_mapping)) {
    from <- layers_mapping[[i]]
    to <- names(layers_mapping)[[i]]
    if (!to %in% c("counts", "data")) {
      SeuratObject::LayerData(
        obj,
        assay = assay_name,
        layer = to
      ) <- adata$layers[[from]]
    }
  }

  # copy reductions
  # trackstatus: class=Seurat, feature=get_obsm, status=wip
  # trackstatus: class=Seurat, feature=get_varm, status=wip
  if (!is.null(reduction_mapping)) {
    if (
      !is.list(reduction_mapping) ||
        (length(reduction_mapping) > 0 && is.null(names(reduction_mapping)))
    ) {
      stop("reduction_mapping must be a named list")
    }
    for (i in seq_along(reduction_mapping)) {
      reduction_name <- names(reduction_mapping)[[i]]
      reduction <- reduction_mapping[[i]]

      if (
        !is.list(reduction) ||
          is.null(names(reduction)) ||
          !all(names(reduction) %in% c("key", "obsm", "varm")) ||
          !all(c("key", "obsm") %in% names(reduction))
      ) {
        stop(
          "each reduction must be a list with keys 'key', 'obsm', and 'varm'"
        )
      }
      dr <- .to_seurat_process_reduction(
        adata = adata,
        key = reduction$key,
        obsm_embedding = reduction$obsm,
        varm_loadings = reduction$varm,
        assay_name = assay_name
      )
      if (!is.null(dr)) {
        obj[[reduction_name]] <- dr
      }
    }
  }

  # trackstatus: class=Seurat, feature=get_obsp, status=wip
  for (i in seq_along(graph_mapping)) {
    graph_name <- names(graph_mapping)[[i]]
    graph <- graph_mapping[[i]]
    if (!is.character(graph) || length(graph) != 1) {
      stop(
        "item ",
        graph_name,
        " in graph_mapping is not a character vector of length 1"
      )
    }
    obsp <- adata$obsp[[graph]]
    if (!is.null(obsp)) {
      dimnames(obsp) <- list(obs_names, obs_names)
      obsp_gr <- Seurat::as.Graph(obsp)
      obj[[paste0(assay_name, "_", graph_name)]] <- obsp_gr
    }
  }

  # trackstatus: class=Seurat, feature=get_uns, status=done
  # trackstatus: class=Seurat, feature=get_varp, status=done
  for (i in seq_along(misc_mapping)) {
    misc_name <- names(misc_mapping)[[i]]
    misc <- misc_mapping[[i]]
    if (!is.character(misc) || length(misc) <= 0 || length(misc) > 2) {
      stop("misc_mapping must be a named list with one or two elements")
    }
    misc_slot <- misc[[1]]
    misc_key <- misc[[2]]
    expected_slots <- c(
      "X",
      "layers",
      "obs",
      "obsm",
      "obsp",
      "var",
      "varm",
      "varp",
      "uns"
    )
    if (!misc_slot %in% expected_slots) {
      stop(paste0(
        "The first element of each item in misc_mapping must be one of: ",
        paste0("'", expected_slots, "'", collapse = ", ")
      ))
    }
    misc_data <- adata[[misc_slot]]
    if (length(misc) == 2) {
      if (!misc_key %in% names(misc_data)) {
        stop(paste0(
          "misc_mapping: adata$",
          misc_slot,
          "[[",
          misc_key,
          "]] does not exist"
        ))
      }
      misc_data <- misc_data[[misc_key]]
    }
    if (!is.null(misc_data)) {
      SeuratObject::Misc(obj, misc_name) <- misc_data
    }
  }

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
    stop("layer_name must be character vector of length 1")
  }

  if (!layer_name %in% names(adata$layers)) {
    stop("layer_name '", layer_name, "' is not an item in adata$layers")
  }

  Matrix::t(adata$layers[[layer_name]])
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

  Matrix::t(adata$layers[[layer_name]])
}

.to_seurat_process_reduction <- function(
  adata,
  assay_name,
  key,
  obsm_embedding,
  varm_loadings
) {
  if (!.to_seurat_is_atomic_character(key)) {
    stop("key must be a character scalar")
  }
  if (!.to_seurat_is_atomic_character(obsm_embedding)) {
    stop("obsm_embedding must be a character scalar")
  }
  if (
    !is.null(varm_loadings) && !.to_seurat_is_atomic_character(varm_loadings)
  ) {
    stop("varm_loadings must be a character scalar or NULL")
  }
  embed <- adata$obsm[[obsm_embedding]]

  if (is.null(embed)) {
    stop(paste0(
      "The reduction ",
      obsm_embedding,
      " is not present in adata$obsm"
    ))
  }

  rownames(embed) <- adata$obs_names

  loadings <-
    if (is.null(varm_loadings)) {
      new(Class = "matrix")
    } else if (!varm_loadings %in% names(adata$varm)) {
      stop(paste0(
        "The loadings ",
        varm_loadings,
        " is not present in adata$varm"
      ))
    } else {
      load <- adata$varm[[varm_loadings]]
      rownames(load) <- adata$var_names
      load
    }

  SeuratObject::CreateDimReducObject(
    embeddings = embed,
    loadings = loadings,
    key = key,
    assay = assay_name,
    global = TRUE
  )
}

# nolint start: object_name_linter object_length_linter
to_Seurat_guess_layers <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  layers <- list()

  if (!is.null(adata$X)) {
    # guess the name of the X slot
    layer_name_for_x <-
      if (!"counts" %in% names(adata$layers)) {
        "counts"
      } else {
        "data"
      }

    layers[layer_name_for_x] <- list(NULL)
  }

  for (layer_name in names(adata$layers)) {
    layers[[layer_name]] <- layer_name
  }

  layers
}

# nolint start: object_name_linter object_length_linter
to_Seurat_guess_reductions <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  reductions <- list()

  for (reduction_name in names(adata$obsm)) {
    if (grepl("^X_", reduction_name)) {
      name <- gsub("^X_", "", reduction_name)
      out <-
        if (reduction_name == "X_pca") {
          list(key = "PC_", obsm = "X_pca", varm = "PCs")
        } else {
          list(key = paste0(name, "_"), obsm = reduction_name, varm = NULL)
        }

      reductions[[name]] <- out
    }
  }

  reductions
}

# nolint start: object_name_linter object_length_linter
to_Seurat_guess_graphs <- function(adata) {
  # nolint end: object_name_linter object_length_linter
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

  graphs
}

# nolint start: object_name_linter object_length_linter
to_Seurat_guess_misc <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  misc_mapping <- list()

  if (!is.null(adata$uns)) {
    for (key in names(adata$uns)) {
      misc_mapping[[key]] <- c("uns", key)
    }
  }

  # TODO: copy obsm which were not used as embeddings?
  # TODO: copy varm which were not used as loadings?
  # Then again, the user can do this manually if needed

  misc_mapping
}

#' Convert a Seurat object to an AnnData object
#'
#' `from_Seurat()` converts a Seurat object to an AnnData object.
#' Only one assay can be converted at a time. Arguments are used to configure the conversion.
#' If `NULL`, the functions `from_Seurat_guess_*` will be used to guess the mapping.
#'
#' For more information on the functionality of an AnnData object, see [anndataR-package].
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_class Name of the AnnData class. Must be one of `"HDF5AnnData"` or `"InMemoryAnnData"`.
#' @param assay_name The name of the assay to be converted. If `NULL`, the default assay will be used
#' ([SeuratObject::DefaultAssay()]).
#' @param x_mapping A mapping of a Seurat layer to the AnnData `X` slot. If `NULL`, no data will be copied to the
#' `X` slot.
#' @param layers_mapping A named list mapping layer names to the names of the layers in the Seurat object. Each item in
#' the list must be a character vector of length 1. See section "`$layers` mapping" for more details.
#' @param obsm_mapping A named list mapping reductions to the names of the reductions in the Seurat object. Each item in
#' the list must be a vector of length 2. See section "`$obsm` mapping" for more details.
#' @param varm_mapping A named list mapping PCA loadings to the names of the PCA loadings in the Seurat object.
#' Each item in the list must be a character vector of length 1. See section "`$varm` mapping" for more details.
#' @param obsp_mapping A named list mapping graph names to the names of the graphs in the Seurat object.
#' Each item in the list must be a character vector of length 1. See section "`$obsp` mapping" for more details.
#' @param varp_mapping A named list mapping miscellaneous data to the names of the data in the Seurat object.
#' Each item in the list must be a named list with one or two elements. See section "`$varp` mapping" for more details.
#' @param uns_mapping A named list mapping miscellaneous data to the names of the data in the Seurat object.
#' Each item in the list must be a named list with one or two elements. See section "`$uns` mapping" for more details.
#' @param ... Additional arguments passed to the generator function.
#'
#' @section `$X` mapping:
#'
#' A mapping of a Seurat layer to the AnnData `X` slot. Its value must be `NULL` or a character vector of the Seurat
#' layer name to copy. If `NULL`, no data will be copied to the `X` slot.
#'
#' @section `$layers` mapping:
#'
#' A named list to map AnnData layers to Seurat layers. Each item in the list must be a character vector of length 1.
#' The `$X` key maps to the `X` slot.
#'
#' Example: `layers_mapping = list(counts = "counts", foo = "bar")`.
#'
#' If `NULL`, the internal function `from_Seurat_guess_layers` will be used to guess the layer mapping as follows:
#'
#' * All Seurat layers are copied to AnnData `layers` by name.
#' * This means that the AnnData `X` slot will be `NULL` (empty). If you want to copy data to the `X` slot,
#'   you must define the layer mapping explicitly.
#'
#' @section `$obsm` mapping:
#'
#' A named list to map Seurat reductions to AnnData `$obsm`.
#'
#' Each item in the list must be a vector of length 2,
#' where the name corresponds to the name of the resulting `$obsm` slot, and the values correspond to the
#' the location of the data in the Seurat object.
#'
#' Example: `obsm_mapping = list(pca = c("reductions", "pca"), umap = c("reductions", "umap"))`.
#'
#' If `NULL`, the internal function `from_Seurat_guess_obsms` will be used to guess the obsm mapping as follows:
#'
#' * All Seurat reductions are prefixed with `X_` and copied to AnnData `$obsm`.
#'
#' @section `$varm` mapping:
#'
#' A named list to map Seurat reduction loadings to AnnData `$varm`.
#'
#' Each item in the list must be a character vector of length 2, where the name corresponds to the name of the
#' resulting `$varm` slot, and the value corresponds to the location of the data in the Seurat object.
#'
#' Example: `varm_mapping = list(PCs = c("reductions", "pca")`.
#'
#' If `NULL`, the internal function `from_Seurat_guess_varms` will be used to guess the varm mapping as follows:
#'
#' * The name of the PCA loadings is copied by name.
#'
#' @section `$obsp` mapping:
#'
#' A named list to map Seurat graphs to AnnData `$obsp`.
#'
#' Example: `obsp_mapping = list(nn = "connectivities")`.
#'
#' If `NULL`, the internal function `from_Seurat_guess_obsps` will be used to guess the obsp mapping as follows:
#'
#' * All Seurat graphs are copied to `$obsp` by name.
#'
#' @section `$varp` mapping:
#'
#' A named list to map Seurat miscellaneous data to AnnData `$varp`. The name of each item corresponds to the
#' resulting `$varp` slot, while the value of each item must be a fector which corresponds to the location of the data
#' in the Seurat object.
#'
#' Example: `varp_mapping = list(foo = c("misc", "foo"))`.
#'
#' If `NULL`, the internal function `from_Seurat_guess_varps` will be used to guess the varp mapping as follows:
#'
#' * No data is mapped to `$varp`.
#'
#' @section `$uns` mapping:
#'
#' A named list to map Seurat miscellaneous data to AnnData `uns`. Each item in the list must be a character of
#' length 2. The first element must be `"misc"`. The second element is the name of the data in the corresponding slot.
#'
#' Example: `uns_mapping = list(foo = c("misc", "foo"))`.
#'
#' If `NULL`, the internal function `from_Seurat_guess_uns` will be used to guess the uns mapping as follows:
#'
#' * All Seurat miscellaneous data is copied to `uns` by name.
#'
#' @return An AnnData object
#'
#' @export
#'
#' @examples
#' library(Seurat)
#'
#' counts <- matrix(rbinom(20000, 1000, .001), nrow = 100)
#' obj <- CreateSeuratObject(counts = counts)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj, npcs = 10L)
#' obj <- FindNeighbors(obj)
#' obj <- RunUMAP(obj, dims = 1:10)
#' from_Seurat(obj)
# nolint start: object_name_linter
from_Seurat <- function(
  # nolint end: object_name_linter
  seurat_obj,
  output_class = c("InMemoryAnnData", "HDF5AnnData"),
  assay_name = NULL,
  x_mapping = NULL,
  layers_mapping = NULL,
  obsm_mapping = NULL,
  varm_mapping = NULL,
  obsp_mapping = NULL,
  varp_mapping = NULL,
  uns_mapping = NULL,
  ...
) {
  check_requires("Converting Seurat to AnnData", "SeuratObject")

  output_class <- match.arg(output_class)

  stopifnot(inherits(seurat_obj, "Seurat"))

  if (is.null(assay_name)) {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  }

  seurat_assay <- seurat_obj@assays[[assay_name]]

  if (is.null(seurat_assay)) {
    stop(paste0(
      "The assay '",
      assay_name,
      "' does not exist in the Seurat object"
    ))
  }

  if (!inherits(seurat_assay, "Assay5")) {
    stop(paste0(
      "Assay '",
      assay_name,
      "' is not a valid Seurat v5 assay.\n",
      "Please use `SeuratObject::UpdateSeuratObject()` to upgrade the object to Seurat v5."
    ))
  }

  if (is.null(layers_mapping)) {
    layers_mapping <- from_Seurat_guess_layers(seurat_obj, assay_name)
  }
  if (is.null(obsm_mapping)) {
    obsm_mapping <- from_Seurat_guess_obsms(seurat_obj, assay_name)
  }
  if (is.null(varm_mapping)) {
    varm_mapping <- from_Seurat_guess_varms(seurat_obj, assay_name)
  }
  if (is.null(obsp_mapping)) {
    obsp_mapping <- from_Seurat_guess_obsps(seurat_obj, assay_name)
  }
  if (is.null(varp_mapping)) {
    varp_mapping <- from_Seurat_guess_varps(seurat_obj)
  }
  if (is.null(uns_mapping)) {
    uns_mapping <- from_Seurat_guess_uns(seurat_obj)
  }

  # fetch obs
  # trackstatus: class=Seurat, feature=set_obs_names, status=done
  # trackstatus: class=Seurat, feature=set_obs, status=done
  obs <- seurat_obj@meta.data

  # fetch var
  # trackstatus: class=Seurat, feature=set_var_names, status=done
  # trackstatus: class=Seurat, feature=set_var, status=done
  var <- seurat_assay@meta.data

  # use generator to create new AnnData object
  generator <- get_anndata_constructor(output_class)

  tryCatch(
    {
      adata <- generator$new(
        obs = obs,
        var = var,
        ...
      )

      # fetch X
      # trackstatus: class=Seurat, feature=set_X, status=done
      if (!is.null(x_mapping)) {
        adata$X <- Matrix::t(seurat_assay@layers[[x_mapping]])
      }

      # fetch layers
      # trackstatus: class=Seurat, feature=set_layers, status=done
      for (i in seq_along(layers_mapping)) {
        layer <- layers_mapping[[i]]
        layer_name <- names(layers_mapping)[[i]]

        adata$layers[[layer_name]] <- Matrix::t(seurat_assay@layers[[layer]])
      }

      # fetch obsm
      # trackstatus: class=Seurat, feature=set_obsm, status=wip
      for (i in seq_along(obsm_mapping)) {
        obsm <- obsm_mapping[[i]]
        obsm_name <- names(obsm_mapping)[[i]]

        if (!is.character(obsm) || length(obsm) != 2) {
          stop("each obsm_mapping must be a character vector of length 2")
        }

        obsm_slot <- obsm[[1]]
        obsm_key <- obsm[[2]]

        if (obsm_slot == "reductions") {
          adata$obsm[[obsm_name]] <- SeuratObject::Embeddings(
            seurat_obj,
            obsm_key
          )
        } else if (obsm_slot == "misc") {
          adata$obsm[[obsm_name]] <- seurat_obj@misc[[obsm_key]]
        }
      }

      # fetch varm
      # trackstatus: class=Seurat, feature=set_varm, status=wip
      for (i in seq_along(varm_mapping)) {
        varm <- varm_mapping[[i]]
        varm_name <- names(varm_mapping)[[i]]

        if (!is.character(varm) || length(varm) < 2 || length(varm) > 3) {
          stop("each varm_mapping must be a character vector of length 2 or 3")
        }

        varm_slot <- varm[[1]]
        varm_key <- varm[[2]]

        if (varm_slot == "reductions") {
          adata$varm[[varm_name]] <- SeuratObject::Loadings(
            seurat_obj,
            varm_key
          )
        } else if (varm_slot == "misc") {
          data <- seurat_obj@misc[[varm_key]]
          if (length(varm) == 3) {
            data <- data[[varm[[3]]]]
          }
          adata$varm[[varm_name]] <- data
        }
      }

      # fetch obsp
      # trackstatus: class=Seurat, feature=set_obsp, status=wip
      for (i in seq_along(obsp_mapping)) {
        obsp <- obsp_mapping[[i]]
        obsp_name <- names(obsp_mapping)[[i]]

        if (!is.character(obsp) || length(obsp) < 2 || length(obsp) > 3) {
          stop("each obsp_mapping must be a character vector of length 2 or 3")
        }

        key1 <- obsp[[1]]
        key2 <- obsp[[2]]

        if (key1 == "graphs") {
          adata$obsp[[obsp_name]] <- as(
            seurat_obj@graphs[[key2]],
            "sparseMatrix"
          )
        } else if (key1 == "misc") {
          data <- seurat_obj@misc[[key2]]
          if (length(obsp) == 3) {
            data <- data[[obsp[[3]]]]
          }
          adata$obsp[[obsp_name]] <- data
        }
      }

      # fetch varp
      # trackstatus: class=Seurat, feature=set_varp, status=wip
      for (i in seq_along(varp_mapping)) {
        varp <- varp_mapping[[i]]
        varp_name <- names(varp_mapping)[[i]]

        if (!is.character(varp) || length(varp) < 2 || length(varp) > 3) {
          stop("each varp_mapping must be a character vector of length 2 or 3")
        }

        key1 <- varp[[1]]
        key2 <- varp[[2]]

        if (key1 == "misc") {
          data <- seurat_obj@misc[[key2]]
          if (length(varp) == 3) {
            data <- data[[varp[[3]]]]
          }
          adata$varp[[varp_name]] <- data
        }
      }

      # fetch uns
      # trackstatus: class=Seurat, feature=set_uns, status=wip
      for (i in seq_along(uns_mapping)) {
        uns <- uns_mapping[[i]]
        uns_name <- names(uns_mapping)[[i]]

        if (!is.character(uns) || length(uns) < 2 || length(uns) > 3) {
          stop("each uns_mapping must be a character vector of length 2 or 3")
        }

        key1 <- uns[[1]]
        key2 <- uns[[2]]

        if (key1 == "misc") {
          data <- seurat_obj@misc[[key2]]
          if (length(uns) == 3) {
            data <- data[[uns[[3]]]]
          }
          adata$uns[[uns_name]] <- data
        }
      }

      adata
    },
    error = function(e) {
      if (output_class == "HDF5AnnData") {
        on.exit(cleanup_HDF5AnnData(adata))
      }
      stop(e)
    }
  )
}

# nolint start: object_name_linter object_length_linter
from_Seurat_guess_layers <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The provided object must be a Seurat object")
  }

  seurat_assay <- seurat_obj@assays[[assay_name]]

  if (!inherits(seurat_assay, "Assay5")) {
    stop("The provided object must be a Seurat v5 assay")
  }

  layers_mapping <- list()

  for (layer_name in SeuratObject::Layers(seurat_assay)) {
    layers_mapping[[layer_name]] <- layer_name
  }

  layers_mapping
}

# nolint start: object_name_linter object_length_linter
from_Seurat_guess_obsms <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The provided object must be a Seurat object")
  }

  obsm_mapping <- list()

  for (reduction_name in SeuratObject::Reductions(seurat_obj)) {
    # Check if the dimreduc was calculated by the selected assay
    reduction <- seurat_obj@reductions[[reduction_name]]
    if (reduction@assay.used != assay_name) {
      next
    }

    obsm_mapping[[paste0("X_", reduction_name)]] <- c(
      "reductions",
      reduction_name
    )
  }

  obsm_mapping
}

# nolint start: object_name_linter object_length_linter
from_Seurat_guess_varms <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The provided object must be a Seurat object")
  }

  varm_mapping <- list()

  if ("pca" %in% SeuratObject::Reductions(seurat_obj)) {
    # Check if the dimreduc was calculated by the selected assay
    reduction <- seurat_obj@reductions[["pca"]]
    if (reduction@assay.used == assay_name) {
      varm_mapping[["PCs"]] <- c("reductions", "pca")
    }
  }

  varm_mapping
}

# nolint start: object_name_linter object_length_linter
from_Seurat_guess_obsps <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The provided object must be a Seurat object")
  }

  obsp_mapping <- list()

  for (graph_name in SeuratObject::Graphs(seurat_obj)) {
    graph <- seurat_obj@graphs[[graph_name]]

    if (graph@assay.used != assay_name) {
      next
    }

    dest_name <- gsub(paste0(assay_name, "_"), "", graph_name)

    if (dest_name == "nn") {
      dest_name <- "connectivities"
    }

    obsp_mapping[[dest_name]] <- c("graphs", graph_name)
  }

  obsp_mapping
}

# nolint start: object_name_linter object_length_linter
from_Seurat_guess_varps <- function(seurat_obj) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The provided object must be a Seurat object")
  }

  list()
}

# nolint start: object_name_linter object_length_linter
from_Seurat_guess_uns <- function(seurat_obj) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The provided object must be a Seurat object")
  }

  uns_mapping <- list()

  for (uns_name in names(seurat_obj@misc)) {
    uns_mapping[[uns_name]] <- c("misc", uns_name)
  }

  uns_mapping
}
