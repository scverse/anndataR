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
#' @param object_metadata_mapping A named list to map observation-level metadata to object-level metadata in the Seurat
#'   object. See section "Metadata mapping" for more details.
#' @param assay_metadata_mapping A named list to map variable-level metadata to assay-level metadata in the Seurat
#'   object.
#'   See section "Metadata mapping" for more details.
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
#' If `NULL`, the internal function `.to_Seurat_guess_layers` will be used to guess the layer mapping as follows:
#'
#' * All AnnData layers are copied to Seurat layers by name.
#'
#' @section Metadata mapping:
#'
#' A named list or vector to map observation-level and feature-level metadata to object-level and assay-level metadata
#' in the Seurat object.
#'
#' Each value in the `object_metadata_mapping` list or vector corresponds to the names of the `obs` slot in the AnnData
#' object, and each of the names correspond to the names of the metadata in the resulting Seurat object.
#'
#' Example: `object_metadata_mapping = c(cellType = "cell_type")`.
#'
#' Each value in the `assay_metadata_mapping` list or vector corresponds to the names of the `var` slot in the AnnData
#' object, and the names correspond to the names of the metadata in the resulting Seurat object.
#'
#' Example: `assay_metadata_mapping = list(geneInfo = "gene_info")`.
#'
#' By default, all metadata in the `obs` and `var` slots will be copied to the Seurat object.
#'
#' @section Reduction mapping:
#'
#' A named list to map AnnData `$obsm` and `$varm` to Seurat reductions. Each item in the list must be a named list
#' with keys `'key'`, `'obsm'`, and can contain the key `'varm'`.
#
#' Example: `reduction_mapping = list(pca = list(key = "PC_", obsm = "X_pca", varm = "PCs"))`.
#'
#' If `NULL`, the internal function `.to_Seurat_guess_reductions` will be used to guess the reduction mapping as
#' follows:
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
#' If `NULL`, the internal function `.to_Seurat_guess_graphs` will be used to guess the graph mapping as follows:
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
#' If `NULL`, the internal function `.to_Seurat_guess_misc` will be used to guess the miscellaneous mapping as follows:
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
  object_metadata_mapping = NULL,
  assay_metadata_mapping = NULL,
  reduction_mapping = NULL,
  graph_mapping = NULL,
  misc_mapping = NULL
) {
  # nolint end: object_name_linter
  check_requires("Converting AnnData to Seurat", "SeuratObject")

  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  if (is.null(object_metadata_mapping)) {
    object_metadata_mapping <- .to_Seurat_guess_object_metadata(adata)
  }
  if (is.null(layers_mapping)) {
    layers_mapping <- .to_Seurat_guess_layers(adata)
  }
  if (is.null(assay_metadata_mapping)) {
    assay_metadata_mapping <- .to_Seurat_guess_assay_metadata(adata)
  }
  if (is.null(reduction_mapping)) {
    reduction_mapping <- .to_Seurat_guess_reductions(adata)
  }
  if (is.null(graph_mapping)) {
    graph_mapping <- .to_Seurat_guess_graphs(adata)
  }
  if (is.null(misc_mapping)) {
    misc_mapping <- .to_Seurat_guess_misc(adata)
  }

  if (length(adata$layers) == 0 && is.null(adata$X)) {
    cli_abort(
      "{.arg adata} must have a valid {.field X} or at least one {.field layer}"
    )
  }

  # store obs and var names
  obs_names <- adata$obs_names[]
  var_names <- adata$var_names[]

  # check seurat layers
  if (is.null(names(layers_mapping))) {
    names(layers_mapping) <- layers_mapping
  }
  if (
    !("counts" %in% names(layers_mapping)) &&
      !("data" %in% names(layers_mapping))
  ) {
    cli_abort(c(
      "{.arg layers_mapping} must contain an item named {.val counts} and/or {.val data}",
      "i" = "Found names: {.val {names(layers_mapping)}}"
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
  object_metadata <- .to_Seurat_process_metadata(
    adata,
    object_metadata_mapping,
    "obs"
  )
  obj <- SeuratObject::CreateSeuratObject(
    meta.data = object_metadata,
    assay = assay_name,
    counts = counts,
    data = data
  )

  # trackstatus: class=Seurat, feature=get_var, status=done
  assay_metadata <- .to_Seurat_process_metadata(
    adata,
    assay_metadata_mapping,
    "var"
  )

  if (ncol(assay_metadata) != 0) {
    obj@assays[[assay_name]] <- SeuratObject::AddMetaData(
      obj@assays[[assay_name]],
      metadata = assay_metadata
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
      ) <- Matrix::t(adata$layers[[from]])
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
      cli_abort(
        "{.arg reduction_mapping} must be a named {.cls list}, got {.cls {class(reduction_mapping)}}"
      )
    }

    reduction_fmt_msg <- paste(
      "Each item in {.arg reduction_mapping} must be a {.cls list}",
      "with names {.val {keys_str}}"
    )
    # nolint start object_usage_linter
    keys_str <- cli::cli_vec(
      c("key", "obsm", "varm"),
      list("vec-last" = " and (optionally) ")
    )
    # nolint end
    for (i in seq_along(reduction_mapping)) {
      reduction_name <- names(reduction_mapping)[[i]]
      reduction <- reduction_mapping[[i]]

      if (!is.list(reduction)) {
        cli_abort(c(
          reduction_fmt_msg,
          "i" = "Item {.val {i}} has class {.cls {class(reduction)}}"
        ))
      }

      if (
        is.null(names(reduction)) ||
          !all(names(reduction) %in% c("key", "obsm", "varm")) ||
          !all(c("key", "obsm") %in% names(reduction))
      ) {
        cli_abort(c(
          reduction_fmt_msg,
          "i" = "Item {.val {i}} has names: {.val {names(reduction)}}"
        ))
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
      cli_abort(c(
        "Each item in {.arg graph_mapping} must be a {.cls character} vector of length 1",
        "i" = "{.code graph_mapping[[{i}]]} is {.obj_type_friendly {graph}}"
      ))
    }
    obsp <- adata$obsp[[graph]]
    if (!is.null(obsp)) {
      dimnames(obsp) <- list(obs_names, obs_names)
      obsp_gr <- Seurat::as.Graph(obsp)
      if (rlang::is_empty(obsp_gr@assay.used)) {
        obsp_gr@assay.used <- assay_name
      }
      obj[[paste0(assay_name, "_", graph_name)]] <- obsp_gr
    }
  }

  # trackstatus: class=Seurat, feature=get_uns, status=done
  # trackstatus: class=Seurat, feature=get_varp, status=done
  for (i in seq_along(misc_mapping)) {
    misc_name <- names(misc_mapping)[[i]]
    misc <- misc_mapping[[i]]
    if (!is.character(misc) || length(misc) <= 0 || length(misc) > 2) {
      cli_abort(c(
        paste(
          "Each item in {.arg misc_mapping} must be a {.cls character} vector",
          "with one or two elements"
        ),
        "i" = "{.code misc_mapping[[{i}]]} is {.obj_type_friendly {misc}}"
      ))
    }

    misc_slot <- misc[1]
    misc_key <- misc[2]
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
      cli_abort(c(
        paste(
          "The first element in each item of {.arg misc_mapping}",
          "must be one of: {.or {.val {expected_slots}}}"
        ),
        "i" = "{.code misc_mapping[[{i}]][1]}: {.val {misc_slot}}"
      ))
    }
    misc_data <- adata[[misc_slot]]
    if (length(misc) == 2) {
      misc_key <- misc[[2]]
      if (!misc_key %in% names(misc_data)) {
        misc_str <- cli::cli_vec(misc, list("vec-last" = ", ")) # nolint object_usage_linter
        cli_abort(paste(
          "The requested item {.code adata${misc_slot}[[{misc_key}]]}",
          "does not exist for {.code misc_mapping[[{i}]]}:",
          "{.val misc_name} = c({.val {misc_str}})"
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
    cli_abort(
      "{.arg layer_name} must be a {.cls character} vector of length 1",
      call = rlang::caller_env()
    )
  }

  if (!layer_name %in% names(adata$layers)) {
    cli_abort(
      "layer name {.val {layer_name}} is not an item in {.code adata$layers}",
      call = rlang::caller_env()
    )
  }

  Matrix::t(adata$layers[[layer_name]])
}

.to_seurat_get_matrix <- function(adata, layer_name) {
  if (is.null(layer_name)) {
    return(Matrix::t(adata$X))
  }

  if (!.to_seurat_is_atomic_character(layer_name)) {
    cli_abort(
      "{.arg layer_name} must be a {.cls character} vector of length 1",
      call = rlang::caller_env()
    )
  }

  if (!layer_name %in% names(adata$layers)) {
    cli_abort(
      "{.arg layer_name} must be the name of a layer or {.val NULL}",
    )
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
    cli_abort(
      "{.arg key} must be a {.cls character} vector of length 1",
      call = rlang::caller_env()
    )
  }
  if (!.to_seurat_is_atomic_character(obsm_embedding)) {
    cli_abort(
      "{.arg obsm_embedding} must be a {.cls character} vector of length 1",
      call = rlang::caller_env()
    )
  }
  if (
    !is.null(varm_loadings) && !.to_seurat_is_atomic_character(varm_loadings)
  ) {
    cli_abort(
      "{.arg varm_loadings} must be a {.cls character} vector of length 1 or {.val NULL}",
      call = rlang::caller_env()
    )
  }
  embed <- adata$obsm[[obsm_embedding]]

  if (is.null(embed)) {
    cli_abort(
      c(
        "The requested item {.val {obsm_embedding}} does not exist in {.code adata$obsm}",
        "i" = "{.code names(adata$obsm)}: {.val {names(adata$obsm)}}"
      ),
      call = rlang::caller_env()
    )
  }

  rownames(embed) <- adata$obs_names

  loadings <-
    if (is.null(varm_loadings)) {
      new(Class = "matrix")
    } else if (!varm_loadings %in% names(adata$varm)) {
      cli_abort(
        c(
          "The requested item {.val {varm_loadings}} does not exist in {.code adata$varm}",
          "i" = "{.code names(adata$varm)}: {.val {names(adata$varm)}}"
        ),
        call = rlang::caller_env()
      )
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
.to_Seurat_guess_layers <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
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
.to_Seurat_guess_reductions <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
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
.to_Seurat_guess_graphs <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
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
.to_Seurat_guess_misc <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  if (!inherits(adata, "AbstractAnnData")) {
    "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
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

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_object_metadata <- function(adata) {
  # nolint end: object_name_linter object_length_linter

  object_metadata_mapping <- as.list(names(adata$obs))
  names(object_metadata_mapping) <- names(adata$obs)

  object_metadata_mapping
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_assay_metadata <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  assay_metadata_mapping <- as.list(names(adata$var))
  names(assay_metadata_mapping) <- names(adata$var)

  assay_metadata_mapping
}

# nolint start: object_name_linter
.to_Seurat_process_metadata <- function(adata, mapping, slot) {
  # nolint end: object_name_linter
  mapped <- adata[[slot]][unlist(mapping)]
  names(mapped) <- names(mapping)
  mapped
}
