#' Convert a SingleCellExperiment object to an AnnData object
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated, use [as_AnnData()] instead
#'
#' @param sce See [as_AnnData()]
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
from_SingleCellExperiment <- function(
  # nolint end: object_name_linter
  sce,
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
    what = "from_SingleCellExperiment()",
    with = "as_AnnData()"
  )

  check_requires(
    "Converting SingleCellExperiment to AnnData",
    "SingleCellExperiment",
    "Bioc"
  )

  output_class <- match.arg(output_class)

  if (!(inherits(sce, "SingleCellExperiment"))) {
    cli_abort(
      "{.arg sce} must be a {.cls SingleCellExperiment} but has class {.cls {sce}}"
    )
  }

  layers_mapping <- get_mapping(
    layers_mapping,
    .from_SCE_guess_layers,
    sce,
    "layers_mapping",
    x_mapping = x_mapping
  )
  obs_mapping <- get_mapping(
    obs_mapping,
    .from_SCE_guess_all,
    sce,
    "obs_mapping",
    slot = SingleCellExperiment::colData
  )
  var_mapping <- get_mapping(
    var_mapping,
    .from_SCE_guess_all,
    sce,
    "var_mapping",
    slot = SingleCellExperiment::rowData
  )
  obsm_mapping <- get_mapping(
    obsm_mapping,
    .from_SCE_guess_obsm,
    sce,
    "obsm_mapping"
  )
  varm_mapping <- get_mapping(
    varm_mapping,
    .from_SCE_guess_varm,
    sce,
    "varm_mapping"
  )
  obsp_mapping <- get_mapping(
    obsp_mapping,
    .from_SCE_guess_obspvarp,
    sce,
    "obsp_mapping",
    slot = SingleCellExperiment::colPairs
  )
  varp_mapping <- get_mapping(
    varp_mapping,
    .from_SCE_guess_obspvarp,
    sce,
    "varp_mapping",
    slot = SingleCellExperiment::rowPairs
  )
  uns_mapping <- get_mapping(
    uns_mapping,
    .from_SCE_guess_all,
    sce,
    "uns_mapping",
    slot = S4Vectors::metadata
  )

  generator <- get_anndata_constructor(output_class)
  adata <- generator$new(shape = rev(dim(sce)), ...)

  # Fill in slots in the object
  .from_SCE_process_obs(adata, sce, obs_mapping)

  .from_SCE_process_var(adata, sce, var_mapping)

  # trackstatus: class=SingleCellExperiment, feature=set_X, status=done
  if (!is.null(x_mapping)) {
    adata$X <- .from_SCE_convert(
      SummarizedExperiment::assay(sce, x_mapping, withDimnames = FALSE)
    )
  }

  .from_SCE_process_layers(adata, sce, layers_mapping)

  .from_SCE_process_obsm(adata, sce, obsm_mapping)

  .from_SCE_process_varm(adata, sce, varm_mapping)

  .from_SCE_process_obsp(adata, sce, obsp_mapping)

  .from_SCE_process_varp(adata, sce, varp_mapping)

  .from_SCE_process_uns(adata, sce, uns_mapping)

  adata
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_all <- function(sce, slot) {
  # nolint end: object_length_linter object_name_linter
  self_name(names(slot(sce)))
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_layers <- function(sce, x_mapping) {
  # nolint end: object_length_linter object_name_linter
  layers_mapping <- self_name(SummarizedExperiment::assayNames(sce))

  if (!is.null(x_mapping)) {
    layers_mapping <- layers_mapping[layers_mapping != x_mapping]
  }

  if (rlang::is_empty(layers_mapping)) {
    c()
  } else {
    layers_mapping
  }
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_obsm <- function(sce) {
  # nolint end: object_length_linter object_name_linter
  obsm_mapping <- self_name(SingleCellExperiment::reducedDimNames(sce))

  if (rlang::is_empty(obsm_mapping)) {
    c()
  } else {
    obsm_mapping
  }
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_varm <- function(sce) {
  # nolint end: object_length_linter object_name_linter
  varm_mapping <- c()

  for (reduction_name in names(SingleCellExperiment::reducedDims(sce))) {
    reduction <- SingleCellExperiment::reducedDim(sce, reduction_name)
    if (inherits(reduction, "LinearEmbeddingMatrix")) {
      varm_mapping[reduction_name] <- reduction_name
    }
  }

  varm_mapping
}

# If sce is a SummarizedExperiment, return an empty mapping
# nolint start: object_length_linter object_name_linter
.from_SCE_guess_obspvarp <- function(sce, slot) {
  # nolint end: object_length_linter object_name_linter
  .from_SCE_guess_all(sce, slot)
}

# Convert BioConductor-specific objects to base R objects
# Convert matrices
# Convert dgCMatrix to dgRMatrix
# nolint start: object_length_linter object_name_linter
.from_SCE_convert <- function(object, transpose = TRUE) {
  # nolint end: object_length_linter object_name_linter
  if (inherits(object, "DataFrame")) {
    as.data.frame(object)
  } else if (inherits(object, "SimpleList")) {
    as.list(object)
  } else if (inherits(object, "matrix") || inherits(object, "Matrix")) {
    if (inherits(object, "denseMatrix")) {
      object <- as.matrix(object)
    }
    if (transpose) {
      to_py_matrix(object)
    } else {
      object
    }
  } else {
    object
  }
}

# trackstatus: class=SingleCellExperiment, feature=set_obs_names, status=done
# trackstatus: class=SingleCellExperiment, feature=set_obs, status=done
# nolint start: object_name_linter
.from_SCE_process_obs <- function(adata, sce, obs_mapping) {
  # nolint end: object_name_linter

  if (is.null(colnames(sce)) && rlang::is_empty(obs_mapping)) {
    adata$obs <- data.frame(matrix(nrow = ncol(sce), ncol = 0))
  } else {
    adata$obs <- purrr::map(obs_mapping, function(.item) {
      # Check if the column exists
      if (!.item %in% names(SingleCellExperiment::colData(sce))) {
        cli_abort(c(
          "Column {.val {obs_mapping[[.item]]}} not found in SCE object.",
          "i" = "Available columns: {.val {names(SingleCellExperiment::colData(sce))}}"
        ))
      }
      SingleCellExperiment::colData(sce)[[.item]]
    }) |>
      as.data.frame(row.names = colnames(sce)) |>
      setNames(names(obs_mapping))
  }
}

# trackstatus: class=SingleCellExperiment, feature=set_var_names, status=done
# trackstatus: class=SingleCellExperiment, feature=set_var, status=done
# nolint start: object_name_linter
.from_SCE_process_var <- function(adata, sce, var_mapping) {
  # nolint end: object_name_linter

  # store an empty data.frame to keep the var names
  if (is.null(rownames(sce)) && rlang::is_empty(var_mapping)) {
    adata$var <- data.frame(matrix(nrow = nrow(sce), ncol = 0))
  } else {
    adata$var <- purrr::map(var_mapping, function(.item) {
      # Check if the column exists
      if (!.item %in% names(SingleCellExperiment::rowData(sce))) {
        cli_abort(c(
          "Column {.val {var_mapping[[.item]]}} not found in SCE object.",
          "i" = "Available columns: {.val {names(SingleCellExperiment::rowData(sce))}}"
        ))
      }
      SingleCellExperiment::rowData(sce)[[.item]]
    }) |>
      as.data.frame(row.names = rownames(sce)) |>
      setNames(names(var_mapping))
  }
}

# trackstatus: class=SingleCellExperiment, feature=set_layers, status=done
# nolint start: object_length_linter object_name_linter
.from_SCE_process_layers <- function(adata, sce, layers_mapping) {
  # nolint end: object_length_linter object_name_linter
  if (rlang::is_empty(layers_mapping)) {
    return(invisible())
  }

  adata$layers <- purrr::map(layers_mapping, function(.layer) {
    .from_SCE_convert(SummarizedExperiment::assay(sce, .layer))
  })
}

# trackstatus: class=SingleCellExperiment, feature=set_obsm, status=done
# nolint start: object_name_linter
.from_SCE_process_obsm <- function(adata, sce, obsm_mapping) {
  # nolint end: object_name_linter
  if (rlang::is_empty(obsm_mapping)) {
    return(invisible())
  }

  adata$obsm <- purrr::imap(obsm_mapping, function(reduction_name, obsm_key) {
    # Check if the reduction exists
    if (!reduction_name %in% names(SingleCellExperiment::reducedDims(sce))) {
      cli_abort(c(
        "Reduction {.val {reduction_name}} not found in SCE object.",
        "i" = "Available reductions: {.val {names(SingleCellExperiment::reducedDims(sce))}}"
      ))
    }

    reduction <- SingleCellExperiment::reducedDim(sce, reduction_name)
    if (inherits(reduction, "LinearEmbeddingMatrix")) {
      SingleCellExperiment::sampleFactors(reduction)
    } else {
      reduction
    }
  })
}

# trackstatus: class=SingleCellExperiment, feature=set_varm, status=done
# nolint start: object_name_linter
.from_SCE_process_varm <- function(adata, sce, varm_mapping) {
  # nolint end: object_name_linter
  if (rlang::is_empty(varm_mapping)) {
    return(invisible())
  }

  adata$varm <- purrr::map(varm_mapping, function(reduction_name) {
    # Check if the reduction exists
    if (!reduction_name %in% names(SingleCellExperiment::reducedDims(sce))) {
      cli_abort(c(
        "Reduction {.val {reduction_name}} not found in SCE object.",
        "i" = "Available reductions: {.val {names(SingleCellExperiment::reducedDims(sce))}}"
      ))
    }

    reduction <- SingleCellExperiment::reducedDim(sce, reduction_name)
    if (!inherits(reduction, "LinearEmbeddingMatrix")) {
      cli_abort(paste(
        "{.code reducedDim(sce, {.val {reduction_name}}} must be a {.cls LinearEmbeddingMatrix}",
        "but has class {.cls {class(reduction)[1]}}"
      ))
    }

    loadings <- SingleCellExperiment::featureLoadings(reduction)
    if (nrow(loadings) != adata$n_vars()) {
      cli_abort(paste(
        "The number of rows ({.val {nrow(loadings)}}) in the loadings for ",
        "reduction {.val {reduction_name}} does not match {.code adata$n_vars()}",
        "({.val {adata$n_vars()}})"
      ))
    }

    loadings
  })
}

# trackstatus: class=SingleCellExperiment, feature=set_obsp, status=done
# nolint start: object_length_linter object_name_linter
.from_SCE_process_obsp <- function(adata, sce, obsp_mapping) {
  # nolint end: object_length_linter object_name_linter
  if (rlang::is_empty(obsp_mapping)) {
    return(invisible())
  }

  adata$obsp <- purrr::map(obsp_mapping, function(colpair_name) {
    # Check if the colPair exists
    if (!(colpair_name %in% SingleCellExperiment::colPairNames(sce))) {
      cli_abort(c(
        "colPair {.val {colpair_name}} not found in SCE object.",
        "i" = "Available colPairs: {.val {names(SingleCellExperiment::colPairs(sce))}}"
      ))
    }

    .from_SCE_convert(
      SingleCellExperiment::colPair(sce, colpair_name, asSparse = TRUE),
      transpose = FALSE
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=set_varp, status=done
# nolint start: object_length_linter object_name_linter
.from_SCE_process_varp <- function(adata, sce, varp_mapping) {
  # nolint end: object_length_linter object_name_linter
  if (rlang::is_empty(varp_mapping)) {
    return(invisible())
  }

  adata$varp <- purrr::map(varp_mapping, function(rowpair_name) {
    # Check if the rowPair exists
    if (!(rowpair_name %in% SingleCellExperiment::rowPairNames(sce))) {
      cli_abort(c(
        "rowPair {.val {rowpair_name}} not found in SCE object.",
        "i" = "Available rowPairs: {.val {names(SingleCellExperiment::rowPairs(sce))}}"
      ))
    }

    .from_SCE_convert(
      SingleCellExperiment::rowPair(sce, rowpair_name, asSparse = TRUE),
      transpose = FALSE
    )
  })
}

# trackstatus: class=SingleCellExperiment, feature=set_uns, status=done
# nolint start: object_length_linter object_name_linter
.from_SCE_process_uns <- function(adata, sce, uns_mapping) {
  # nolint end: object_length_linter object_name_linter
  if (rlang::is_empty(uns_mapping)) {
    return(invisible())
  }

  adata$uns <- purrr::map(uns_mapping, function(meta_name) {
    # Check if the metadata exists
    if (!(meta_name %in% names(S4Vectors::metadata(sce)))) {
      cli_abort(c(
        "Metadata {.val {meta_name}} not found in SCE object.",
        "i" = "Available metadata: {.val {names(S4Vectors::metadata(sce))}}"
      ))
    }

    S4Vectors::metadata(sce)[[meta_name]]
  })
}
