#' Convert a Seurat object to an AnnData object
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated, use [as_AnnData()] instead
#'
#' @param seurat_obj See [as_AnnData()]
#' @param assay_name See [as_AnnData()]
#' @param x_mapping See [as_AnnData()]
#' @param layers_mapping See [as_AnnData()]
#' @param obs_mapping See [as_AnnData()]
#' @param var_mapping See [as_AnnData()]
#' @param obsm_mapping See [as_AnnData()]
#' @param varm_mapping See [as_AnnData()]
#' @param obsp_mapping See [as_AnnData()]
#' @param varp_mapping See [as_AnnData()]
#' @param uns_mapping See [as_AnnData()]
#' @param output_class See [as_AnnData()]
#' @param ... See [as_AnnData()]
#'
#' @return An `AnnData` object
#' @export
# nolint start: object_name_linter
from_Seurat <- function(
  # nolint end: object_name_linter
  seurat_obj,
  assay_name = NULL,
  x_mapping = NULL,
  layers_mapping = TRUE,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  obsm_mapping = TRUE,
  varm_mapping = TRUE,
  obsp_mapping = TRUE,
  varp_mapping = TRUE,
  uns_mapping = TRUE,
  output_class = c("InMemory", "HDF5AnnData"),
  ...
) {
  lifecycle::deprecate_soft(
    when = "0.99.0",
    what = "from_Seurat()",
    with = "as_AnnData()"
  )

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

  layers_mapping <- get_mapping(
    layers_mapping,
    .from_Seurat_guess_layers,
    seurat_obj,
    "layers_mapping",
    assay_name = assay_name
  )
  obs_mapping <- get_mapping(
    obs_mapping,
    .from_Seurat_guess_obs,
    seurat_obj,
    "obs_mapping",
    assay_name = assay_name
  )
  var_mapping <- get_mapping(
    var_mapping,
    .from_Seurat_guess_var,
    seurat_obj,
    "var_mapping",
    assay_name = assay_name
  )
  obsm_mapping <- get_mapping(
    obsm_mapping,
    .from_Seurat_guess_obsms,
    seurat_obj,
    "obsm_mapping",
    assay_name = assay_name
  )
  varm_mapping <- get_mapping(
    varm_mapping,
    .from_Seurat_guess_varms,
    seurat_obj,
    "varm_mapping",
    assay_name = assay_name
  )
  obsp_mapping <- get_mapping(
    obsp_mapping,
    .from_Seurat_guess_obsps,
    seurat_obj,
    "obsp_mapping",
    assay_name = assay_name
  )
  varp_mapping <- get_mapping(
    varp_mapping,
    .from_Seurat_guess_varps,
    seurat_obj,
    "varp_mapping"
  )
  uns_mapping <- get_mapping(
    uns_mapping,
    .from_Seurat_guess_uns,
    seurat_obj,
    "uns_mapping"
  )

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
