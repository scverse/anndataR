#' Convert an `AnnData` to a `Seurat`
#'
#' Convert an `AnnData` object to a `Seurat` object
#'
#' @param adata The `AnnData` object to convert.
#' @param assay_name Name of the assay to be created in the new `Seurat` object
#' @param layers_mapping A named vector where names are names of `Layers` in the resulting
#'   `Seurat` object and values are keys of `layers` in `adata`.
#'   See below for details.
#' @param object_metadata_mapping A named vector where names are cell
#'   metadata columns in the resulting `Seurat` object and values are columns of `obs` in `adata`.
#'   See below for details.
#' @param assay_metadata_mapping A named vector where names are variable
#'   metadata columns in the assay of the resulting `Seurat` object and values are columns of `var` in `adata`.
#'   See below for details.
#' @param reduction_mapping A named vector where names are names of `Embeddings` in the resulting
# `Seurat` object and values are keys of `obsm` in `adata`. Alternatively, a named list where names are names of `Embeddings`
#' in the resulting `Seurat` object and values are vectors with the items `"key"`, `"obsm"` and
#' (optionally) `"varm"`. See below for details.
#' @param graph_mapping  A named vector where names are names of `Graphs` in the resulting
#'   `Seurat` object and values are keys of `obsp` in `adata`.
#'   See below for details.
#' @param misc_mapping A named vector where names are names of `Misc` in the resulting
#'   `Seurat` object and values are keys of `uns` in `adata`.
#'   See below for details.
#'
#' @details
#'
#'   ## Mapping arguments
#'
#'   All mapping arguments expect a named character
#'   vector where names are the names of the slot in the `Seurat` object and
#'   values are the keys of the corresponding slot of `adata`. If `NULL`,
#'   the conversion function will guess which items to copy as described in the
#'   conversion tables conversion table below. In most cases, the default is to
#'   copy all items using the same names except where the correspondence between
#'   objects is unclear. The `reduction_mapping` argument can also accept a more complex list format, see below for details. To avoid copying anything to a slot, provide an empty
#'   vector. If an unnamed vector is provided, the values will be used as names
#'
#'   ### Examples:
#'
#'   - `NULL` will guess which items to copy as described in the conversion
#'     table
#'   - `c(seurat_item = "adata_item")` will copy `adata_item` from the slot in `adata` to
#'     `seurat_item` in the corresponding slot of new `Seurat` object
#'   - `c()` will avoid copying anything to the slot
#'   - `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`
#'
#'   ## Conversion table
#'
#'   | **From `AnnData`** | **To `Seurat`** | **Example mapping argument** | **Default if `NULL`** |
#'   |--------------------|-------------------------------|------------------------------|-----------------------|
#'   | `adata$layers` | `Layers(seurat)` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
#'   | `adata$obs` | `seurat[[]]` | `object_metadata_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#'   | `adata$var` | `seurat[[assay_name]][[]]` | `assay_metadata_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#'   | `adata$obsm` | `Embeddings(sce)` | `reduction_mapping = c(pca = "X_pca")` **OR** `reduction_mapping = list(pca = c(key = "PC_", obsm = "X_pca", varm = "PCs"))`  | All items that can be coerced to a numeric matrix are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
#'   | `adata$obsp` | `Graphs(seurat)` | `graph_mapping = c(nn = "connectivities")` | All items are copied by name |
#'   | `adata$varp` | _NA_  | _NA_ | There is no corresponding slot for `varp` |
#'   | `adata$uns` | `Misc(seurat)` | `misc_mapping = c(project_metadata = "metadata")` | All items are copied by name |
#'
#' ## The `reduction_mapping` argument
#'
#' For the simpler named vector format, the names should be the names of
#'   `Embeddings` in the resulting `Seurat` object, and the values
#'   should be the keys of `obsm` in `adata`.
#'
#' For more advanced mapping, use the list format where each item is a vector with the
#'   following names:
#'
#'   - `key`: the key of the resulting [`SeuratObject::DimReduc`] object
#'   - `obsm`: a key of the `obsm` slot in `adata`
#'   - `varm`: a key of the `varm` slot in `adata` (optional)
#'
#' All keys are processed by [SeuratObject::Key()] to match requirements
#'
#' @return A `Seurat` object containing the requested data from `adata`
#' @keywords internal
#'
#' @examplesIf rlang::is_installed("Seurat")
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' ad$as_Seurat()
#'
#' @importFrom Matrix t
# nolint start: object_name_linter
as_Seurat <- function(
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

  object_metadata_mapping <- self_name(object_metadata_mapping) %||%
    .to_Seurat_guess_object_metadata(adata)
  layers_mapping <- self_name(layers_mapping) %||%
    .to_Seurat_guess_layers(adata)
  assay_metadata_mapping <- self_name(assay_metadata_mapping) %||%
    .to_Seurat_guess_assay_metadata(adata)
  reduction_mapping <- self_name(reduction_mapping) %||%
    .to_Seurat_guess_reductions(adata)
  graph_mapping <- self_name(graph_mapping) %||% .to_Seurat_guess_graphs(adata)
  misc_mapping <- self_name(misc_mapping) %||% .to_Seurat_guess_misc(adata)

  if (length(adata$layers) == 0 && is.null(adata$X)) {
    cli_abort(
      "{.arg adata} must have a valid {.field X} or at least one {.field layer}"
    )
  }

  # store obs and var names
  obs_names <- adata$obs_names
  var_names <- adata$var_names

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

  # copy other layers
  for (i in seq_along(layers_mapping)) {
    from <- layers_mapping[[i]]
    to <- names(layers_mapping)[[i]]
    if (!to %in% c("counts", "data")) {
      SeuratObject::LayerData(
        obj,
        assay = assay_name,
        layer = to
      ) <- to_R_matrix(adata$layers[[from]])
    }
  }

  if (!rlang::is_empty(reduction_mapping)) {
    reductions <- .to_Seurat_process_reduction_mapping(
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

.to_seurat_is_atomic_character <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x)
}

.to_seurat_get_matrix_by_key <- function(adata, mapping, key) {
  if (!key %in% names(mapping)) {
    return(NULL)
  }

  layer_name <- mapping[[key]]

  .to_seurat_get_matrix(adata, layer_name)
}

.to_seurat_get_matrix <- function(adata, layer_name) {
  if (is.na(layer_name)) {
    return(to_R_matrix(adata$X))
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

  # check if dgRMatrix and convert to dgCMatrix
  to_R_matrix(adata$layers[[layer_name]])
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
.to_Seurat_process_reduction_mapping <- function(
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
        obsm = .obsm
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
        !all(names(reduction) %in% c("key", "obsm", "varm")) ||
        !all(c("key", "obsm") %in% names(reduction))
    ) {
      cli_abort(c(
        paste(
          "Each item in {.arg reduction_mapping} must be a {.cls character} vector",
          "with names {style_vec(c('key', 'obsm', 'varm'), last = ' and (optionally) ')}"
        ),
        "i" = "Item {.val {i}} has names: {.val {names(reduction)}}"
      ))
    }

    varm_loadings <- reduction["varm"]
    if (is.na(varm_loadings)) {
      varm_loadings <- NULL
    }

    dr <- .to_seurat_process_reduction(
      adata = adata,
      assay_name = assay_name,
      key = reduction["key"],
      obsm_embedding = reduction["obsm"],
      varm_loadings = varm_loadings
    )

    if (!is.null(dr)) {
      reductions[[reduction_name]] <- dr
    }
  }

  reductions
}

# nolint start: object_name_linter
.to_Seurat_process_metadata <- function(adata, mapping, slot) {
  # nolint end: object_name_linter
  mapped <- adata[[slot]][unlist(mapping)]
  names(mapped) <- names(mapping)
  mapped
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_layers <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  layers <- self_name(adata$layers_keys())

  if (!is.null(adata$X)) {
    # guess the name of the X slot
    layer_name_for_x <-
      if (!"counts" %in% adata$layers_keys()) {
        "counts"
      } else {
        "data"
      }

    layers[layer_name_for_x] <- NA
  }

  layers
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_reductions <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  purrr::map(adata$obsm_keys(), function(.obsm) {
    if (!is.numeric(as.matrix(adata$obsm[[.obsm]]))) {
      return(NULL)
    }

    mapping <- c(
      # Make sure we have valid keys here to avoid warnings later
      key = suppressWarnings(SeuratObject::Key(.obsm)),
      obsm = .obsm
    )
    if (.obsm == "X_pca" && "PCs" %in% names(adata$varm)) {
      mapping["varm"] <- "PCs"
    }

    mapping
  }) |>
    setNames(adata$obsm_keys()) |>
    purrr::compact()
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_graphs <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$obsp_keys())
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_misc <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$uns_keys())
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_object_metadata <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$obs_keys())
}

# nolint start: object_name_linter object_length_linter
.to_Seurat_guess_assay_metadata <- function(adata) {
  # nolint end: object_name_linter object_length_linter
  self_name(adata$var_keys())
}

# See as_AnnData() for function documentation
from_Seurat <- function(
    seurat_obj,
    assay_name = NULL,
    x_mapping = NULL,
    layers_mapping = NULL,
    obs_mapping = NULL,
    var_mapping = NULL,
    obsm_mapping = NULL,
    varm_mapping = NULL,
    obsp_mapping = NULL,
    varp_mapping = NULL,
    uns_mapping = NULL,
    output_class = c("InMemory", "HDF5AnnData"),
    ...
) {
  check_requires("Converting Seurat to AnnData", c("SeuratObject", "Seurat"))

  output_class <- match.arg(output_class)

  if (!inherits(seurat_obj, "Seurat")) {
    cli_abort(
      "{.arg seurat_obj} must be a {.cls Seurat} object but has class {.cls {class(seurat_obj)}}"
    )
  }

  if (is.null(assay_name)) {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  }

  if (!(assay_name %in% SeuratObject::Assays(seurat_obj))) {
    cli_abort(c(
      "{.arg assay_name} is not an assay in {.arg seurat_obj}",
      "i" = "{.code Assays(seurat_obj)}: {.val {SeuratObject::Assays(seurat_obj)}}"
    ))
  }

  # Set the default assay so we can easily get the dimensions etc.
  SeuratObject::DefaultAssay(seurat_obj) <- assay_name

  # For any mappings that are not set, using the guessing function
  layers_mapping <- self_name(layers_mapping) %||%
    .from_Seurat_guess_layers(seurat_obj, assay_name)
  obs_mapping <- self_name(obs_mapping) %||%
    .from_Seurat_guess_obs(seurat_obj, assay_name)
  var_mapping <- self_name(var_mapping) %||%
    .from_Seurat_guess_var(seurat_obj, assay_name)
  obsm_mapping <- self_name(obsm_mapping) %||%
    .from_Seurat_guess_obsms(seurat_obj, assay_name)
  varm_mapping <- self_name(varm_mapping) %||%
    .from_Seurat_guess_varms(seurat_obj, assay_name)
  obsp_mapping <- self_name(obsp_mapping) %||%
    .from_Seurat_guess_obsps(seurat_obj, assay_name)
  varp_mapping <- self_name(varp_mapping) %||%
    .from_Seurat_guess_varps(seurat_obj)
  uns_mapping <- self_name(uns_mapping) %||% .from_Seurat_guess_uns(seurat_obj)

  generator <- get_anndata_constructor(output_class)
  adata <- generator$new(shape = rev(dim(seurat_obj)), ...)

  # Fill in slots in the object
  .from_Seurat_process_obs(
    adata,
    seurat_obj,
    assay_name,
    obs_mapping
  )

  .from_Seurat_process_var(
    adata,
    seurat_obj,
    assay_name,
    var_mapping
  )

  # trackstatus: class=Seurat, feature=set_X, status=done
  if (!is.null(x_mapping)) {
    adata$X <- to_py_matrix(SeuratObject::LayerData(seurat_obj, x_mapping))
  }

  .from_Seurat_process_layers(
    adata,
    seurat_obj,
    assay_name,
    layers_mapping
  )

  .from_Seurat_process_obsm(
    adata,
    seurat_obj,
    assay_name,
    obsm_mapping
  )

  .from_Seurat_process_varm(
    adata,
    seurat_obj,
    assay_name,
    varm_mapping
  )

  .from_Seurat_process_obsp(
    adata,
    seurat_obj,
    assay_name,
    obsp_mapping
  )

  .from_Seurat_process_varp(
    adata,
    seurat_obj,
    assay_name,
    varp_mapping
  )

  .from_Seurat_process_uns(
    adata,
    seurat_obj,
    assay_name,
    uns_mapping
  )

  adata
}

# trackstatus: class=Seurat, feature=set_obs_names, status=done
# trackstatus: class=Seurat, feature=set_obs, status=done
# nolint start: object_name_linter
.from_Seurat_process_obs <- function(
  adata,
  seurat_obj,
  assay_name,
  obs_mapping
) {
  # nolint end: object_name_linter

  if (!rlang::is_empty(obs_mapping)) {
    if (!all(obs_mapping %in% names(seurat_obj[[]]))) {
      missing <- setdiff(obs_mapping, names(seurat_obj[[]])) # nolint object_usage_linter
      cli_abort(paste(
        "The requested obs item(s) {.val {missing}} do not exist in the Seurat",
        "object ({.code names(seurat_obj[[]])})"
      ))
    }
    adata$obs <- seurat_obj[[unlist(obs_mapping)]] |>
      setNames(names(obs_mapping))
  } else {
    # Store an empty data.frame to keep the obs names
    adata$obs <- data.frame(row.names = colnames(seurat_obj))
  }
}

# trackstatus: class=Seurat, feature=set_var_names, status=done
# trackstatus: class=Seurat, feature=set_var, status=done
# nolint start: object_name_linter
.from_Seurat_process_var <- function(
  adata,
  seurat_obj,
  assay_name,
  var_mapping
) {
  # nolint end: object_name_linter
  if (!rlang::is_empty(var_mapping)) {
    if (!all(var_mapping %in% names(seurat_obj[[assay_name]][[]]))) {
      missing <- setdiff(var_mapping, names(seurat_obj[[assay_name]][[]])) # nolint object_usage_linter
      cli_abort(paste(
        "The requested var item(s) {.val {missing}} do not exist in the Seurat",
        "object ({.code names(seurat_obj[[{.val {assay_name}}]][[]])})"
      ))
    }
    adata$var <- seurat_obj[[assay_name]][[unlist(var_mapping)]] |>
      setNames(names(var_mapping))
  } else {
    # Store an empty data.frame to keep the var names
    adata$var <- data.frame(row.names = rownames(seurat_obj))
  }
}

# trackstatus: class=Seurat, feature=set_layers, status=done
# nolint start: object_name_linter
.from_Seurat_process_layers <- function(
  # nolint end: object_name_linter
  adata,
  seurat_obj,
  assay_name,
  layers_mapping
) {
  if (rlang::is_empty(layers_mapping)) {
    return(invisible())
  }

  adata$layers <- purrr::map(layers_mapping, function(.layer) {
    to_py_matrix(
      SeuratObject::LayerData(seurat_obj, assay = assay_name, layer = .layer)
    )
  })
}

# trackstatus: class=Seurat, feature=set_obsm, status=done
# nolint start: object_name_linter
.from_Seurat_process_obsm <- function(
  adata,
  seurat_obj,
  assay_name,
  obsm_mapping
) {
  # nolint end: object_name_linter
  if (rlang::is_empty(obsm_mapping)) {
    return(invisible())
  }

  adata$obsm <- purrr::map(obsm_mapping, function(.reduction) {
    if (!(.reduction %in% SeuratObject::Reductions(seurat_obj))) {
      cli_abort(c(
        "Reduction {.val {.reduction}} not found in Seurat object.",
        "i" = "Available reductions: {.val {SeuratObject::Reductions(seurat_obj)}}"
      ))
    }
    SeuratObject::Embeddings(seurat_obj, .reduction)
  })
}

# trackstatus: class=Seurat, feature=set_varm, status=done
# nolint start: object_name_linter
.from_Seurat_process_varm <- function(
  adata,
  seurat_obj,
  assay_name,
  varm_mapping
) {
  # nolint end: object_name_linter
  if (rlang::is_empty(varm_mapping)) {
    return(invisible())
  }

  adata$varm <- purrr::map(varm_mapping, function(.reduction) {
    if (!(.reduction %in% SeuratObject::Reductions(seurat_obj))) {
      cli_abort(c(
        "Reduction {.val {.reduction}} not found in Seurat object.",
        "i" = "Available reductions: {.val {SeuratObject::Reductions(seurat_obj)}}"
      ))
    }
    SeuratObject::Loadings(seurat_obj, .reduction)
  })
}

# trackstatus: class=Seurat, feature=set_obsp, status=done
# nolint start: object_name_linter
.from_Seurat_process_obsp <- function(
  adata,
  seurat_obj,
  assay_name,
  obsp_mapping
) {
  # nolint end: object_name_linter
  if (rlang::is_empty(obsp_mapping)) {
    return(invisible())
  }

  adata$obsp <- purrr::map(obsp_mapping, function(.graph) {
    if (!(.graph %in% SeuratObject::Graphs(seurat_obj))) {
      cli_abort(c(
        "Graph {.val {graph_name}} not found in Seurat object.",
        "i" = "Available graphs: {.val {SeuratObject::Graphs(seurat_obj)}}"
      ))
    }
    as(seurat_obj[[.graph]], "sparseMatrix")
  })
}

# trackstatus: class=Seurat, feature=set_varp, status=done
# nolint start: object_name_linter
.from_Seurat_process_varp <- function(
  adata,
  seurat_obj,
  assay_name,
  varp_mapping
) {
  # nolint end: object_name_linter
  if (rlang::is_empty(varp_mapping)) {
    return(invisible())
  }

  adata$varp <- purrr::map(varp_mapping, function(.varp) {
    # Check if the misc data exists
    if (!(.varp %in% names(SeuratObject::Misc(seurat_obj)))) {
      cli_abort(c(
        "Misc data {.val {.varp}} not found in Seurat object.",
        "i" = "Available misc data: {.val {names(SeuratObject::Misc(seurat_obj))}}"
      ))
    }
    SeuratObject::Misc(seurat_obj, .varp)
  })
}

# trackstatus: class=Seurat, feature=set_uns, status=done
# nolint start: object_name_linter
.from_Seurat_process_uns <- function(
  adata,
  seurat_obj,
  assay_name,
  uns_mapping
) {
  # nolint end: object_name_linter
  if (rlang::is_empty(uns_mapping)) {
    return(invisible())
  }

  adata$uns <- purrr::map(uns_mapping, function(.misc) {
    if (!(.misc %in% names(SeuratObject::Misc(seurat_obj)))) {
      cli_abort(c(
        "Misc data {.val {.misc}} not found in Seurat object.",
        "i" = "Available misc data: {.val {names(SeuratObject::Misc(seurat_obj))}}"
      ))
    }
    SeuratObject::Misc(seurat_obj, .misc)
  })
}

# nolint start: object_name_linter object_length_linter
.from_Seurat_guess_layers <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter
  seurat_assay <- Seurat::GetAssay(seurat_obj, assay = assay_name)
  layers <- SeuratObject::Layers(seurat_assay)
  setNames(layers, layers)
}

# nolint start: object_name_linter
.from_Seurat_guess_obs <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter
  setNames(names(seurat_obj[[]]), names(seurat_obj[[]]))
}

# nolint start: object_name_linter
.from_Seurat_guess_var <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter
  assay <- seurat_obj[[assay_name]]
  setNames(names(assay[[]]), names(assay[[]]))
}

# nolint start: object_name_linter object_length_linter
.from_Seurat_guess_obsms <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter

  obsm_mapping <- c()

  for (reduction_name in SeuratObject::Reductions(seurat_obj)) {
    # Check if the dimreduc was calculated by the selected assay
    reduction <- seurat_obj[[reduction_name]]
    if (SeuratObject::DefaultAssay(reduction) != assay_name) {
      next
    }

    obsm_mapping[reduction_name] <- reduction_name
  }

  obsm_mapping
}

# nolint start: object_name_linter object_length_linter
.from_Seurat_guess_varms <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter

  varm_mapping <- c()

  for (reduction_name in SeuratObject::Reductions(seurat_obj)) {
    reduction <- seurat_obj[[reduction_name]]
    loadings <- SeuratObject::Loadings(seurat_obj, reduction_name)

    if (
      !SeuratObject::IsMatrixEmpty(loadings) &&
        nrow(loadings) == nrow(seurat_obj) &&
        SeuratObject::DefaultAssay(reduction) == assay_name
    ) {
      varm_mapping[reduction_name] <- reduction_name
    }
  }

  varm_mapping
}

# nolint start: object_name_linter object_length_linter
.from_Seurat_guess_obsps <- function(seurat_obj, assay_name) {
  # nolint end: object_name_linter object_length_linter

  obsp_mapping <- c()

  for (graph_name in SeuratObject::Graphs(seurat_obj)) {
    graph <- seurat_obj[[graph_name]]

    if (SeuratObject::DefaultAssay(graph) != assay_name) {
      next
    }

    dest_name <- gsub(paste0(assay_name, "_"), "", graph_name)

    obsp_mapping[dest_name] <- graph_name
  }

  obsp_mapping
}

# nolint start: object_name_linter object_length_linter
.from_Seurat_guess_varps <- function(seurat_obj) {
  # nolint end: object_name_linter object_length_linter
  c()
}

# nolint start: object_name_linter object_length_linter
.from_Seurat_guess_uns <- function(seurat_obj) {
  # nolint end: object_name_linter object_length_linter
  uns_names <- names(SeuratObject::Misc(seurat_obj))
  setNames(uns_names, uns_names)
}
