#' Convert an `AnnData` to a `SingleCellExperiment`
#'
#' Convert an `AnnData` object to a `SingleCellExperiment` object
#'
#' @param adata The `AnnData` object to convert
#' @param assays_mapping A named vector where names are names of `assays` in the resulting
#'   `SingleCellExperiment` object and values are keys of `layers` in `adata`.
#'   See below for details.
#' @param colData_mapping A named vector where names are columns of `colData` in the resulting
#'   `SingleCellExperiment` object and values are columns of `obs` in `adata`.
#'   See below for details.
#' @param rowData_mapping A named vector where names are columns of `rowData` in the resulting
#'   `SingleCellExperiment` object and values are columns of `var` in `adata`.
#'   See below for details.
#' @param reduction_mapping A named vector where names are names of `reducedDims` in the resulting
#'   `SingleCellExperiment` object and values are keys of `obsm` in `adata`. Alternatively, a named list where names are names of `reducedDims` in the resulting `SingleCellExperiment` object and values are vectors with the items `"obsm"` and `"varm"` and/or `"uns"`. See below for details.
#' @param colPairs_mapping A named vector where names are names of `colPairs` in the resulting
#'   `SingleCellExperiment` object and values are keys of `obsp` in `adata`.
#'   See below for details.
#' @param rowPairs_mapping A named vector where names are names of `rowPairs` in the resulting
#'   `SingleCellExperiment` object and values are keys of `varp` in `adata`.
#'   See below for details.
#' @param metadata_mapping A named vector where names are names of `metadata` in the resulting
#'   `SingleCellExperiment` object and values are keys of `uns` in `adata`.
#'   See below for details.
#'
#' @details
#'
#'   ## Mapping arguments
#'
#'   All mapping arguments expect a named character
#'   vector where names are the names of the slot in the `SingleCellExperiment` object and
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
#'   - `c(sce_item = "adata_item")` will copy `adata_item` from the slot in `adata` to
#'     `sce_item` in the corresponding slot of new `SingleCellExperiment` object
#'   - `c()` will avoid copying anything to the slot
#'   - `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`
#'
#'   ## Conversion table
#'
#'   | **From `AnnData`** | **To `SingleCellExperiment`** | **Example mapping argument** | **Default if `NULL`** |
#'   |--------------------|-------------------------------|------------------------------|-----------------------|
#'   | `adata$layers` | `assays(sce)` | `assays_mapping = c(counts = "counts")` | All items are copied by name |
#'   | `adata$obs` | `colData(sce)` | `colData_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#'   | `adata$var` | `rowData(sce)` | `rowData_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#'   | `adata$obsm` | `reducedDims(sce)` | `reduction_mapping = c(pca = "X_pca")` **OR** `reduction_mapping = list(pca = c(obsm = "X_pca", varm = "PCs", uns = "pca_metadata"))`  | All items are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
#'   | `adata$obsp` | `colPairs(sce)` | `colPairs_mapping = c(nn = "connectivities")` | All items are copied by name |
#'   | `adata$varp` | `rowPairs(sce)` | `rowPairs_mapping = c(gene_overlaps = "similarities")` | All items are copied by name |
#'   | `adata$uns` | `metadata(sce)` | `uns_mapping = c(project_metadata = "metadata")` | All items are copied by name |
#'
#' ## The `reduction_mapping` argument
#'
#' For the simpler named vector format, the names should be the names of
#'   `reducedDims` in the resulting `SingleCellExperiment` object, and the values
#'   should be the keys of `obsm` in `adata`.
#'
#'   For more advanced mapping, use the list format where each item is a vector with the
#'   following names:
#'
#'   - `obsm`: a key of the `obsm` slot in `adata`
#'   - `varm`: a key of the `varm` slot in `adata` (optional)
#'   - `uns`: a key of the `uns` slot in `adata` (optional)
#'
#'   If `'varm'` or `'uns'` is given, a [SingleCellExperiment::LinearEmbeddingMatrix]
#'   will be created for that item with `adata$varm[[varm]]` passed to the `featureLoadings` argument and
#'   `adata$uns[[uns]]` passed as `metadata`
#'
#' @return A `SingleCellExperiment` object containing the requested data from `adata`
#' @keywords internal
#'
#' @examplesIf rlang::is_installed("SingleCellExperiment")
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#'
#' sce <- ad$as_SingleCellExperiment()
as_SingleCellExperiment <- function(
  adata,
  assays_mapping = NULL,
  colData_mapping = NULL, # nolint
  rowData_mapping = NULL, # nolint
  reduction_mapping = NULL,
  colPairs_mapping = NULL, # nolint
  rowPairs_mapping = NULL, # nolint
  metadata_mapping = NULL
) {
  check_requires(
    "Converting AnnData to SingleCellExperiment",
    "SingleCellExperiment",
    "Bioc"
  )

  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  # guess mappings if not provided
  # nolint start object_name_linter
  assays_mapping <- self_name(assays_mapping) %||% .to_SCE_guess_assays(adata)
  colData_mapping <- self_name(colData_mapping) %||%
    .to_SCE_guess_all(adata, "obs")
  rowData_mapping <- self_name(rowData_mapping) %||%
    .to_SCE_guess_all(adata, "var")
  reduction_mapping <- self_name(reduction_mapping) %||%
    .to_SCE_guess_reduction(adata)
  colPairs_mapping <- self_name(colPairs_mapping) %||%
    .to_SCE_guess_all(adata, "obsp")
  rowPairs_mapping <- self_name(rowPairs_mapping) %||%
    .to_SCE_guess_all(adata, "varp")
  metadata_mapping <- self_name(metadata_mapping) %||%
    .to_SCE_guess_all(adata, "uns")
  # nolint end object_name_linter

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  sce_assays <- vector("list", length(assays_mapping))
  names(sce_assays) <- names(assays_mapping)
  for (i in seq_along(assays_mapping)) {
    from <- assays_mapping[[i]]
    to <- names(assays_mapping)[[i]]
    if (from != "X") {
      sce_assays[[to]] <- to_R_matrix(adata$layers[[from]])
    } else {
      sce_assays[[to]] <- to_R_matrix(adata$X)
    }
  }

  # construct colData
  # FIXME: probably better way to make a dataframe from a list of vectors
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- .to_SCE_process_simple_mapping(adata, colData_mapping, "obs")

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- .to_SCE_process_simple_mapping(adata, rowData_mapping, "var")

  # construct colPairs
  # trackstatus: class=SingleCellExperiment, feature=get_obsp, status=done
  col_pairs <- .to_SCE_process_simple_mapping(adata, colPairs_mapping, "obsp")

  # construct rowPairs
  # trackstatus: class=SingleCellExperiment, feature=get_varp, status=done
  row_pairs <- .to_SCE_process_simple_mapping(adata, rowPairs_mapping, "varp")

  # construct metadata
  # trackstatus: class=SingleCellExperiment, feature=get_uns, status=done
  metadata <- .to_SCE_process_simple_mapping(adata, metadata_mapping, "uns")

  # construct reducedDims
  reduceddims <- .to_SCE_process_reduction_mapping(adata, reduction_mapping)

  arguments <- list(
    assays = sce_assays,
    colPairs = col_pairs,
    rowPairs = row_pairs,
    metadata = metadata,
    checkDimnames = TRUE
  )
  # add col_data if not empty list
  if (length(col_data) > 0) {
    arguments$colData <- as(col_data, "DataFrame")
  }
  # add row_data if not empty list
  if (length(row_data) > 0) {
    arguments$rowData <- as(row_data, "DataFrame")
  }

  # construct output object
  sce <- do.call(SingleCellExperiment::SingleCellExperiment, arguments)
  rownames(sce) <- rownames(adata$var)
  colnames(sce) <- rownames(adata$obs)

  SingleCellExperiment::reducedDims(sce) <- reduceddims # only here to ensure that the dimensions are right

  sce
}

# nolint start: object_length_linter object_name_linter
.to_SCE_guess_assays <- function(adata) {
  # nolint end: object_length_linter object_name_linter
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  layers <- list()

  if (!is.null(adata$X)) {
    layer_name_for_x <-
      if (!"counts" %in% names(adata$layers)) {
        # could expand checks, to check if integers
        "counts"
      } else {
        "data"
      }
    layers[[layer_name_for_x]] <- "X"
  }

  for (layer_name in names(adata$layers)) {
    layers[[layer_name]] <- layer_name
  }

  layers
}

# nolint start: object_length_linter object_name_linter
.to_SCE_guess_all <- function(adata, slot) {
  # nolint end: object_length_linter object_name_linter
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  self_name(names(adata[[slot]]))
}

# nolint start: object_length_linter object_name_linter
.to_SCE_guess_reduction <- function(adata) {
  # nolint end: object_length_linter object_name_linter
  .to_Seurat_guess_reductions(adata) # nolint object_usage_linter
}

# nolint start: object_length_linter object_name_linter
.to_SCE_process_simple_mapping <- function(adata, mapping, slot) {
  # nolint end: object_length_linter object_name_linter

  if (rlang::is_empty(mapping)) {
    return(list())
  }

  purrr::map(mapping, function(.item) {
    adata[[slot]][[.item]]
  })
}

# trackstatus: class=SingleCellExperiment, feature=get_obsm, status=done
# trackstatus: class=SingleCellExperiment, feature=get_varm, status=done
# nolint start: object_length_linter object_name_linter
.to_SCE_process_reduction <- function(
  # nolint end: object_length_linter object_name_linter
  adata,
  obsm_key,
  varm_key,
  uns_key
) {
  if (!(obsm_key %in% adata$obsm_keys())) {
    cli_abort(
      c(
        "{.val {obsm_key}} is not an item in {.code adata$obsm}",
        "i" = "{.code adata$obsm_keys()}: {.val {adata$obsm_keys()}}"
      )
    )
  }

  embedding <- adata$obsm[[obsm_key]]
  rownames(embedding) <- adata$obs_names

  # If there is no extra info to add, return the embedding matrix
  if (is.null(varm_key) && is.null(uns_key)) {
    return(embedding)
  }

  if (!is.null(varm_key)) {
    if (!(varm_key %in% adata$varm_keys())) {
      cli_abort(
        c(
          "{.val {varm_key}} is not an item in {.code adata$varm}",
          "i" = "{.code adata$varm_keys()}: {.val {adata$varm_keys()}}"
        )
      )
    }

    loadings <- adata$varm[[varm_key]]
    rownames(loadings) <- colnames(embedding)
  } else {
    loadings <- matrix(nrow = 0, ncol = ncol(embedding))
  }

  if (!is.null(uns_key)) {
    if (!(uns_key %in% adata$uns_keys())) {
      cli_abort(
        c(
          "{.val {uns_key}} is not an item in {.code adata$uns}",
          "i" = "{.code adata$uns_keys()}: {.val {adata$uns_keys()}}"
        )
      )
    }

    metadata <- adata$uns[[uns_key]]
  } else {
    metadata <- list()
  }

  lem <- SingleCellExperiment::LinearEmbeddingMatrix(
    sampleFactors = embedding,
    featureLoadings = loadings,
    metadata = metadata
  )
  rownames(lem) <- rownames(embedding)
  colnames(lem) <- colnames(loadings)

  lem
}

# nolint start: object_length_linter object_name_linter
.to_SCE_process_reduction_mapping <- function(adata, reduction_mapping) {
  # nolint end: object_length_linter object_name_linter

  # If the mapping is a vector expand it to the list format
  if (is.atomic(reduction_mapping)) {
    reduction_mapping <- purrr::map(reduction_mapping, function(.obsm) {
      c(obsm = .obsm)
    })
  } else {
    purrr::walk(names(reduction_mapping), function(.name) {
      mapping <- reduction_mapping[[.name]]
      if (!(is.character(mapping) && ("obsm" %in% names(mapping)))) {
        cli_abort(c(
          paste(
            "If it is a list, each item in {.arg reduction_mapping} must be a
            {.cls character} vector with the name {.val obsm} and optional",
            "names {.val varm} and {.val uns}"
          ),
          "i" = paste(
            "{.code reduction_mapping[[{.val { .name }}]]}:",
            "{style_vec(mapping, wrap = TRUE)}"
          )
        ))
      }
    })
  }

  purrr::map(reduction_mapping, function(.map) {
    varm <- if ("varm" %in% names(.map)) .map[["varm"]] else NULL
    uns <- if ("uns" %in% names(.map)) .map[["uns"]] else NULL

    .to_SCE_process_reduction(
      adata,
      obsm_key = .map[["obsm"]],
      varm_key = varm,
      uns_key = uns
    )
  })
}

# See as_AnnData() for function documentation
from_SingleCellExperiment <- function(
    sce,
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

  # For any mappings that are not set, using the guessing function
  layers_mapping <- self_name(layers_mapping) %||%
    .from_SCE_guess_layers(sce, x_mapping)
  obs_mapping <- self_name(obs_mapping) %||%
    .from_SCE_guess_all(
      sce,
      SingleCellExperiment::colData
    )
  var_mapping <- self_name(var_mapping) %||%
    .from_SCE_guess_all(
      sce,
      SingleCellExperiment::rowData
    )
  obsm_mapping <- self_name(obsm_mapping) %||% .from_SCE_guess_obsm(sce)
  varm_mapping <- self_name(varm_mapping) %||% .from_SCE_guess_varm(sce)
  obsp_mapping <- self_name(obsp_mapping) %||%
    .from_SCE_guess_obspvarp(
      sce,
      SingleCellExperiment::colPairs
    )
  varp_mapping <- self_name(varp_mapping) %||%
    .from_SCE_guess_obspvarp(
      sce,
      SingleCellExperiment::rowPairs
    )
  uns_mapping <- self_name(uns_mapping) %||%
    .from_SCE_guess_all(sce, S4Vectors::metadata)

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

  if (!rlang::is_empty(obs_mapping)) {
    adata$obs <- SummarizedExperiment::colData(sce) |>
      as.data.frame() |>
      setNames(names(obs_mapping))
  } else {
    # Store an empty data.frame to keep the obs names
    if (is.null(colnames(sce))) {
      adata$obs <- data.frame(matrix(nrow = ncol(sce), ncol = 0))
    } else {
      adata$obs <- data.frame(row.names = colnames(sce))
    }
  }
}

# trackstatus: class=SingleCellExperiment, feature=set_var_names, status=done
# trackstatus: class=SingleCellExperiment, feature=set_var, status=done
# nolint start: object_name_linter
.from_SCE_process_var <- function(adata, sce, var_mapping) {
  # nolint end: object_name_linter

  if (!rlang::is_empty(var_mapping)) {
    adata$var <- SummarizedExperiment::rowData(sce) |>
      as.data.frame() |>
      setNames(names(var_mapping))
  } else {
    # Store an empty data.frame to keep the var names
    if (is.null(rownames(sce))) {
      adata$var <- data.frame(matrix(nrow = nrow(sce), ncol = 0))
    } else {
      adata$var <- data.frame(row.names = rownames(sce))
    }
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
