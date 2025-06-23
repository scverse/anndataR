#' Convert an `AnnData` to a `Seurat`
#'
#' Convert an `AnnData` object to a `Seurat` object
#'
#' @param adata The `AnnData` object to convert.
#' @param assay_name Name of the assay to be created in the new `Seurat` object
#' @param x_mapping A string specifying the name of the layer in the resulting
#'   `Seurat` object where the data in the `X` slot of `adata` will be mapped to
#' @param layers_mapping A named vector where names are names of `Layers` in the
#'   resulting `Seurat` object and values are keys of `layers` in `adata`. See
#'   below for details.
#' @param object_metadata_mapping A named vector where names are cell metadata
#'   columns in the resulting `Seurat` object and values are columns of `obs` in
#'   `adata`. See below for details.
#' @param assay_metadata_mapping A named vector where names are variable
#'   metadata columns in the assay of the resulting `Seurat` object and values
#'   are columns of `var` in `adata`. See below for details.
#' @param reduction_mapping A named vector where names are names of `Embeddings`
#'   in the resulting `Seurat` object and values are keys of `obsm` in `adata`.
#'   Alternatively, a named list where names are names of `Embeddings` in the
#'   resulting `Seurat` object and values are vectors with the items `"key"`,
#'   `"embeddings"` and (optionally) `"loadings"`. See below for details.
#' @param graph_mapping  A named vector where names are names of `Graphs` in the
#'   resulting `Seurat` object and values are keys of `obsp` in `adata`. See
#'   below for details.
#' @param misc_mapping A named vector where names are names of `Misc` in the
#'   resulting `Seurat` object and values are keys of `uns` in `adata`. See
#'   below for details.
#'
#' @details
#'
#' ## Mapping arguments
#'
#' All mapping arguments expect a named character vector where names are the
#' names of the slot in the `Seurat` object and values are the keys of the
#' corresponding slot of `adata`. If `TRUE`, the conversion function will guess
#' which items to copy as described in the conversion table below. In most
#' cases, the default is to copy all items using the same names except where the
#' correspondence between objects is unclear. The `reduction_mapping` argument
#' can also accept a more complex list format, see below for details. To avoid
#' copying anything to a slot, set the mapping argument to `FALSE`. Empt
#'  mapping arguments (`NULL`, `c()`, `list()`) will be treated as `FALSE` with
#' a warning. If an unnamed vector is provided, the values will be used as
#' names.
#'
#' ### Examples:
#'
#' - `TRUE` will guess which items to copy as described in the conversion
#'   table
#' - `c(seurat_item = "adata_item")` will copy `adata_item` from the slot in
#'   `adata` to `seurat_item` in the corresponding slot of the new `Seurat`
#'   object
#' - `FALSE` will avoid copying anything to the slot
#' - `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`
#'
#' ## Conversion table
#'
# nolint start: line_length_linter
#'
#' | **From `AnnData`** | **To `Seurat`** | **Example mapping argument** | **Default if `NULL`** |
#' |--------------------|-------------------------------|------------------------------|-----------------------|
#' | `adata$X` | `Layers(seurat)` | `x_mapping = "counts"` _OR_ `layers_mapping = c(counts = NA)` | The data in `adata$X` is copied to a layer named `X` |
#' | `adata$layers` | `Layers(seurat)` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
#' | `adata$obs` | `seurat[[]]` | `object_metadata_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#' | `adata$var` | `seurat[[assay_name]][[]]` | `assay_metadata_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#' | `adata$obsm` | `Embeddings(sce)` | `reduction_mapping = c(pca = "X_pca")` **OR** `reduction_mapping = list(pca = c(key = "PC_", obsm = "X_pca", varm = "PCs"))`  | All items that can be coerced to a numeric matrix are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
#' | `adata$obsp` | `Graphs(seurat)` | `graph_mapping = c(nn = "connectivities")` | All items are copied by name |
#' | `adata$varp` | _NA_  | _NA_ | There is no corresponding slot for `varp` |
#' | `adata$uns` | `Misc(seurat)` | `misc_mapping = c(project_metadata = "metadata")` | All items are copied by name |
#'
# nolint end: line_length_linter
#'
#' ## The `reduction_mapping` argument
#'
#' For the simpler named vector format, the names should be the names of
#' `Embeddings` in the resulting `Seurat` object, and the values
#' should be the keys of `obsm` in `adata`. A key will created from the `obsm`
#' key.
#'
#' For more advanced mapping, use the list format where each item is a vector
#' with the following names defining arguments to
#' [SeuratObject::CreateDimReducObject()]:
#'
#' - `key`: the key of the resulting [`SeuratObject::DimReduc`] object, passed
#'   to the `key` argument after processing with [SeuratObject::Key()]
#' - `embeddings`: a key of the `obsm` slot in `adata`,
#'   `adata$obsm[[embeddings]]` is passed to the `embeddings` argument
#' - `loadings`: a key of the `varm` slot in `adata` (optional),
#'   `adata$varm[[loadings]]` is passed to the `loadings` argument
#'
#' ## The `x_mapping` and `layers_mapping` arguments
#'
#' In order to specify where the data in `adata$X` will be stored in the
#' `Layers(seurat)` slot of the resulting object, you can use either the `x_mapping`
#' argument or the `layers_mapping` argument.
#' If you use `x_mapping`, it should be a string specifying the name of the layer
#' in `Layers(seurat)` where the data in `adata$X` will be stored.
#' If you use `layers_mapping`, it should be a named vector where names are names
#' of `Layers(seurat)` and values are keys of `layers` in `adata`.
#' In order to indicate the `adata$X` slot, you use `NA` as the value in the vector.
#' The name you provide for `x_mapping` may not be a name in `layers_mapping`.
#'
#' @return A `Seurat` object containing the requested data from `adata`
#' @keywords internal
#'
#' @family object converters
#'
#' @examplesIf rlang::is_installed("Seurat")
#'   ad <- AnnData(
#'     X = matrix(1:5, 3L, 5L),
#'     obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'     var = data.frame(row.names = letters[1:5], gene = 1:5)
#'   )
#'
#'   # Default usage
#'   seurat <- ad$as_Seurat(
#'     assay_name = "RNA",
#'     x_mapping = "counts",
#'     layers_mapping = TRUE,
#'     object_metadata_mapping = TRUE,
#'     assay_metadata_mapping = TRUE,
#'     reduction_mapping = TRUE,
#'     graph_mapping = TRUE,
#'     misc_mapping = TRUE
#'   )
#'
#' @importFrom Matrix t
# nolint start: object_name_linter
as_Seurat <- function(
  adata,
  assay_name = "RNA",
  x_mapping = NULL,
  layers_mapping = TRUE,
  object_metadata_mapping = TRUE,
  assay_metadata_mapping = TRUE,
  reduction_mapping = TRUE,
  graph_mapping = TRUE,
  misc_mapping = TRUE
) {
  # nolint end: object_name_linter
  check_requires("Converting AnnData to Seurat", c("Seurat", "SeuratObject"))

  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  object_metadata_mapping <- get_mapping(
    object_metadata_mapping,
    .as_Seurat_guess_object_metadata,
    adata,
    "object_metadata_mapping"
  )
  layers_mapping <- get_mapping(
    layers_mapping,
    .as_Seurat_guess_layers,
    adata,
    "layers_mapping"
  )
  assay_metadata_mapping <- get_mapping(
    assay_metadata_mapping,
    .as_Seurat_guess_assay_metadata,
    adata,
    "assay_metadata_mapping"
  )
  reduction_mapping <- get_mapping(
    reduction_mapping,
    .as_Seurat_guess_reductions,
    adata,
    "reduction_mapping"
  )
  graph_mapping <- get_mapping(
    graph_mapping,
    .as_Seurat_guess_graphs,
    adata,
    "graph_mapping"
  )
  misc_mapping <- get_mapping(
    misc_mapping,
    .as_Seurat_guess_misc,
    adata,
    "misc_mapping"
  )

  if (length(adata$layers) == 0 && is.null(adata$X)) {
    cli_abort(
      "{.arg adata} must have a valid {.field X} or at least one {.field layer}"
    )
  }

  # store obs and var names
  obs_names <- adata$obs_names
  var_names <- adata$var_names

  # check seurat layers (which includes the X mapping)
  if (!is.null(x_mapping)) {
    layers_mapping <- setNames(
      c(NA, layers_mapping),
      c(x_mapping, names(layers_mapping))
    )
  } else if (!is.null(adata$X)) {
    layers_mapping[["X"]] <- NA
  }
  if (any(duplicated(names(layers_mapping)))) {
    cli_abort(
      "{.arg layers_mapping} or {.arg x_mapping} must not contain any duplicate names",
      "i" = "Found duplicate names: {.val {names(layers_mapping)[duplicated(names(layers_mapping))]}}"
    )
  }

  object_metadata <- .as_Seurat_process_metadata(
    adata,
    object_metadata_mapping,
    "obs"
  )

  obj <- .as_Seurat_create_object_with_layers(
    adata,
    layers_mapping,
    object_metadata,
    assay_name
  )

  # trackstatus: class=Seurat, feature=get_var, status=done
  assay_metadata <- .as_Seurat_process_metadata(
    adata,
    assay_metadata_mapping,
    "var"
  )

  if (ncol(assay_metadata) != 0) {
    obj[[assay_name]] <- SeuratObject::AddMetaData(
      obj[[assay_name]],
      metadata = assay_metadata
    )
  }

  # make sure obs and var names are set properly
  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  colnames(obj) <- obs_names
  rownames(obj) <- var_names

  if (!rlang::is_empty(reduction_mapping)) {
    reductions <- .as_Seurat_process_reduction_mapping(
      adata,
      reduction_mapping,
      assay_name
    )
    for (.red in names(reductions)) {
      obj[[.red]] <- reductions[[.red]]
    }
  }

  # trackstatus: class=Seurat, feature=get_obsp, status=done
  if (!rlang::is_empty(graph_mapping)) {
    for (i in seq_along(graph_mapping)) {
      graph_name <- names(graph_mapping)[[i]]
      graph <- graph_mapping[[i]]

      if (!(graph %in% adata$obsp_keys())) {
        cli_abort(c(
          "The requested item {.val {graph}} does not exist in {.code adata$obsp}",
          "i" = "{.code adata$obsp_keys()}: {.val {adata$obsp_keys()}}"
        ))
      }

      obsp <- adata$obsp[[graph]]
      dimnames(obsp) <- list(obs_names, obs_names)
      obsp_gr <- Seurat::as.Graph(obsp)
      SeuratObject::DefaultAssay(obsp_gr) <- assay_name
      obj[[paste0(assay_name, "_", graph_name)]] <- obsp_gr
    }
  }

  # trackstatus: class=Seurat, feature=get_uns, status=done
  if (!rlang::is_empty(misc_mapping)) {
    misc_names <- names(misc_mapping)

    for (i in seq_along(misc_names)) {
      misc_name <- misc_names[[i]]
      uns_key <- misc_mapping[[i]]

      if (!(uns_key %in% adata$uns_keys())) {
        cli_abort(c(
          "The requested item {.val {uns_key}} does not exist in {.code adata$uns}",
          "i" = "{.code adata$uns_keys()}: {.val {adata$uns_keys()}}"
        ))
      }

      SeuratObject::Misc(obj, misc_name) <- adata$uns[[uns_key]]
    }
  }

  obj
}

#' Convert an AnnData object to a Seurat object
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated, use `adata$as_Seurat()` instead
#'
#' @param ... Arguments passed to [as_Seurat()]
#'
#' @return A `Seurat` object
#' @export
# nolint start: object_name_linter
to_Seurat <- function(...) {
  # nolint end: object_name_linter
  lifecycle::deprecate_warn(
    when = "0.99.0",
    what = "to_Seurat()",
    with = "adata$as_Seurat()"
  )
  as_Seurat(...)
}

# nolint start: object_name_linter
.as_Seurat_is_atomic_character <- function(x) {
  # nolint end: object_name_linter
  is.character(x) && length(x) == 1 && !is.na(x)
}

# nolint start: object_name_linter
.as_Seurat_get_matrix_by_key <- function(adata, mapping, key) {
  # nolint end: object_name_linter
  if (!key %in% names(mapping)) {
    return(NULL)
  }
  layer_name <- mapping[[key]]

  .as_Seurat_get_matrix(adata, layer_name)
}

# nolint start: object_name_linter
.as_Seurat_get_matrix <- function(adata, layer_name) {
  # nolint end: object_name_linter
  if (is.na(layer_name)) {
    return(to_R_matrix(adata$X))
  }

  if (!.as_Seurat_is_atomic_character(layer_name)) {
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

  # check if dgRMatrix and convert to dgCMatrix
  to_R_matrix(adata$layers[[layer_name]])
}

# nolint start: object_name_linter
.as_Seurat_create_object_with_layers <- function(
  # nolint end: object_name_linter
  adata,
  layers_mapping,
  object_metadata,
  assay_name
) {
  if (!"counts" %in% names(layers_mapping)) {
    # If there's no counts in the layers_mapping,
    # we will consider the first layer as "dummy counts"
    # And use it temporarily as the counts layer
    dummy_counts <- names(layers_mapping)[[1]]
  } else {
    dummy_counts <- "counts"
  }

  counts <- .as_Seurat_get_matrix_by_key(adata, layers_mapping, dummy_counts)
  dimnames(counts) <- list(adata$var_names, adata$obs_names)

  obj <- SeuratObject::CreateSeuratObject(
    meta.data = object_metadata,
    assay = assay_name,
    counts = counts
  )

  # If we have used the dummy counts layer,
  # we need to add the actual dummy counts layer to the object
  # and remove the counts layer
  if (!"counts" %in% names(layers_mapping)) {
    SeuratObject::LayerData(obj, layer = dummy_counts) <- counts
    obj[[assay_name]]$counts <- NULL
  }

  # Add all other layers to the object
  for (i in seq_along(layers_mapping)) {
    layer_name <- names(layers_mapping)[i]
    if (layer_name != dummy_counts) {
      SeuratObject::LayerData(
        obj,
        layer = layer_name,
      ) <- .as_Seurat_get_matrix_by_key(adata, layers_mapping, layer_name)
    }
  }

  obj
}

# nolint start: object_name_linter
.as_Seurat_process_reduction <- function(
  # nolint end: object_name_linter
  adata,
  assay_name,
  key,
  obsm_embedding,
  varm_loadings
) {
  if (!.as_Seurat_is_atomic_character(key)) {
    cli_abort(
      "{.arg key} must be a {.cls character} vector of length 1",
      call = rlang::caller_env()
    )
  }

  if (!.as_Seurat_is_atomic_character(obsm_embedding)) {
    cli_abort(
      "{.arg obsm_embedding} must be a {.cls character} vector of length 1",
      call = rlang::caller_env()
    )
  }

  if (
    !is.null(varm_loadings) && !.as_Seurat_is_atomic_character(varm_loadings)
  ) {
    cli_abort(
      paste(
        "{.arg varm_loadings} must be a {.cls character} vector of length 1 or",
        "{.val NULL}, but is a {.cls {class(varm_loadings)[1]}}"
      ),
      call = rlang::caller_env()
    )
  }

  if (!(obsm_embedding %in% adata$obsm_keys())) {
    cli_abort(
      c(
        "The requested item {.val {obsm_embedding}} does not exist in {.code adata$obsm}",
        "i" = "{.code adata$obsm_keys()}: {.val {adata$obsm_keys()}}"
      ),
      call = rlang::caller_env()
    )
  }

  embed <- as.matrix(adata$obsm[[obsm_embedding]])
  rownames(embed) <- adata$obs_names

  if (!is.numeric(embed)) {
    cli_abort(
      "Embedding matrix for key {.val {key}} must be {.cls numeric} but is {.cls {typeof(embed)}}"
    )
  }

  loadings <-
    if (is.null(varm_loadings)) {
      new(Class = "matrix")
    } else if (!(varm_loadings %in% adata$varm_keys())) {
      cli_abort(
        c(
          "The requested item {.val {varm_loadings}} does not exist in {.code adata$varm}",
          "i" = "{.code adata$varm_keys()}: {.val {adata$varm_keys()}}"
        ),
        call = rlang::caller_env()
      )
    } else {
      load <- as.matrix(adata$varm[[varm_loadings]])
      rownames(load) <- adata$var_names
      load
    }

  if (!grepl(SeuratObject::.KeyPattern(), key)) {
    new_key <- suppressWarnings(Seurat::Key(key))
    cli_warn(paste(
      "Key {.val {key}} does not match the expected pattern,",
      "it has been replaced with {.val {new_key}}"
    ))
    key <- new_key
  }

  allowed_colnames <- paste0(key, seq_len(ncol(embed)))
  if (!identical(colnames(embed), allowed_colnames)) {
    if (!rlang::is_empty(colnames(embed))) {
      cli_warn(paste(
        "Embedding column names do not match what is allowed for key {.val {key}},",
        "setting them to the allowed column names"
      ))
    }
    colnames(embed) <- allowed_colnames
  }

  SeuratObject::CreateDimReducObject(
    embeddings = embed,
    loadings = loadings,
    key = key,
    assay = assay_name,
    global = TRUE
  )
}

# trackstatus: class=Seurat, feature=get_obsm, status=done
# trackstatus: class=Seurat, feature=get_varm, status=done
# nolint start: object_length_linter object_name_linter
.as_Seurat_process_reduction_mapping <- function(
  adata,
  reduction_mapping,
  assay_name
) {
  # nolint end: object_length_linter object_name_linter

  # If reduction_mapping is a vector convert it to the list format
  if (is.atomic(reduction_mapping)) {
    reduction_mapping <- purrr::map(reduction_mapping, function(.obsm) {
      c(
        key = suppressWarnings(SeuratObject::Key(.obsm)),
        embeddings = .obsm
      )
    })
  }

  reductions <- list()

  # Process each reduction
  for (i in seq_along(reduction_mapping)) {
    reduction_name <- names(reduction_mapping)[[i]]
    reduction <- reduction_mapping[[i]]

    if (
      is.null(names(reduction)) ||
        !all(names(reduction) %in% c("key", "embeddings", "loadings")) ||
        !all(c("key", "embeddings") %in% names(reduction))
    ) {
      cli_abort(c(
        paste(
          "Each item in {.arg reduction_mapping} must be a {.cls character} vector",
          "with names {style_vec(c('key', 'embeddings', 'loadings'),",
          "last = ' and (optionally) ')}"
        ),
        "i" = "Item {.val {i}} has names: {.val {names(reduction)}}"
      ))
    }

    varm_loadings <- reduction["loadings"]
    if (is.na(varm_loadings)) {
      varm_loadings <- NULL
    }

    dr <- .as_Seurat_process_reduction(
      adata = adata,
      assay_name = assay_name,
      key = reduction["key"],
      obsm_embedding = reduction["embeddings"],
      varm_loadings = varm_loadings
    )

    if (!is.null(dr)) {
      reductions[[reduction_name]] <- dr
    }
  }

  reductions
}

# nolint start: object_name_linter
.as_Seurat_process_metadata <- function(adata, mapping, slot) {
  # nolint end: object_name_linter
  mapped <- adata[[slot]][unlist(mapping)]
  names(mapped) <- names(mapping)
  mapped
}

# nolint start: object_name_linter object_length_linter
.as_Seurat_guess_layers <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$layers_keys())
}


# nolint start: object_name_linter
.as_Seurat_process_metadata <- function(adata, mapping, slot) {
  # nolint end: object_name_linter
  mapped <- adata[[slot]][unlist(mapping)]
  names(mapped) <- names(mapping)
  mapped
}


# nolint start: object_name_linter object_length_linter
.as_Seurat_guess_reductions <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  purrr::map(adata$obsm_keys(), function(.obsm) {
    if (!is.numeric(as.matrix(adata$obsm[[.obsm]]))) {
      return(NULL)
    }

    mapping <- c(
      # Make sure we have valid keys here to avoid warnings later
      key = suppressWarnings(SeuratObject::Key(.obsm)),
      embeddings = .obsm
    )
    if (.obsm == "X_pca" && "PCs" %in% names(adata$varm)) {
      mapping["loadings"] <- "PCs"
    }

    mapping
  }) |>
    setNames(adata$obsm_keys()) |>
    purrr::compact()
}

# nolint start: object_name_linter object_length_linter
.as_Seurat_guess_graphs <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$obsp_keys())
}

# nolint start: object_name_linter object_length_linter
.as_Seurat_guess_misc <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$uns_keys())
}

# nolint start: object_name_linter object_length_linter
.as_Seurat_guess_object_metadata <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$obs_keys())
}

# nolint start: object_name_linter object_length_linter
.as_Seurat_guess_assay_metadata <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$var_keys())
}
