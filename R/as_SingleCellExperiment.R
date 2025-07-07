#' Convert an `AnnData` to a `SingleCellExperiment`
#'
#' Convert an `AnnData` object to a `SingleCellExperiment` object
#'
#' @param adata The `AnnData` object to convert
#' @param x_mapping A string specifying the name of the assay in the resulting
#'   `SingleCellExperiment` where the data in the `X` slot of `adata` will be
#'   mapped to
#' @param assays_mapping A named vector where names are names of `assays` in the
#'   resulting `SingleCellExperiment` object and values are keys of `layers` in
#'   `adata`. See below for details.
#' @param colData_mapping A named vector where names are columns of `colData` in
#'   the resulting `SingleCellExperiment` object and values are columns of `obs`
#'   in `adata`. See below for details.
#' @param rowData_mapping A named vector where names are columns of `rowData` in
#'   the resulting `SingleCellExperiment` object and values are columns of `var`
#'   in `adata`. See below for details.
#' @param reducedDims_mapping A named vector where names are names of
#'   `reducedDims` in the resulting `SingleCellExperiment` object and values are
#'   keys of `obsm` in `adata`. Alternatively, a named list where names are
#'   names of `reducedDims` in the resulting `SingleCellExperiment` object and
#'   values are vectors with the items `"sampleFactors"` and `"featureLoadings"`
#'   and/or `"metadata"`. See below for details.
#' @param colPairs_mapping A named vector where names are names of `colPairs` in
#'   the resulting `SingleCellExperiment` object and values are keys of `obsp`
#'   in `adata`. See below for details.
#' @param rowPairs_mapping A named vector where names are names of `rowPairs` in
#'   the resulting `SingleCellExperiment` object and values are keys of `varp`
#'   in `adata`. See below for details.
#' @param metadata_mapping A named vector where names are names of `metadata` in
#'   the resulting `SingleCellExperiment` object and values are keys of `uns` in
#'   `adata`. See below for details.
#'
#' @details
#'
#' ## Mapping arguments
#'
#' All mapping arguments expect a named character vector where names are the
#' names of the slot in the `SingleCellExperiment` object and values are the
#' keys of the corresponding slot of `adata`. If `TRUE`, the conversion function
#' will guess which items to copy as described in the conversion tables below.
#' In most cases, the default is to copy all items using the same names except
#' where the correspondence between objects is unclear. The
#' `reducedDims_mapping` argument can also accept a more complex list format,
#' see below for details. To avoid copying anything to a slot, set the mapping
#' argument to `FALSE`. Empty mapping arguments (`NULL`, `c()`, `list()`) will
#' be treated as `FALSE` with a warning. If an unnamed vector is provided, the
#' values will be used as names.
#'
#' ### Examples:
#'
#' - `TRUE` will guess which items to copy as described in the conversion
#'   table
#' - `c(sce_item = "adata_item")` will copy `adata_item` from the slot in
#'   `adata` to `sce_item` in the corresponding slot of the new
#'   `SingleCellExperiment` object
#' - `FALSE` will avoid copying anything to the slot
#' - `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`
#'
#' ## Conversion table
#'
# nolint start: line_length_linter
#'
#' | **From `AnnData`** | **To `SingleCellExperiment`** | **Example mapping argument** | **Default if `NULL`** |
#' |--------------------|-------------------------------|------------------------------|-----------------------|
#' | `adata$X` | `assays(sce)` | `x_mapping = "counts"` | The data in `adata$X` is copied to the assay named `X` |
#' | `adata$layers` | `assays(sce)` | `assays_mapping = c(counts = "counts")` | All items are copied by name |
#' | `adata$obs` | `colData(sce)` | `colData_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#' | `adata$var` | `rowData(sce)` | `rowData_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#' | `adata$obsm` | `reducedDims(sce)` | `reducedDims_mapping = c(pca = "X_pca")` **OR** `reducedDims_mapping = list(pca = c(obsm = "X_pca", varm = "PCs", uns = "pca_metadata"))`  | All items are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
#' | `adata$obsp` | `colPairs(sce)` | `colPairs_mapping = c(nn = "connectivities")` | All items are copied by name |
#' | `adata$varp` | `rowPairs(sce)` | `rowPairs_mapping = c(gene_overlaps = "similarities")` | All items are copied by name |
#' | `adata$uns` | `metadata(sce)` | `uns_mapping = c(project_metadata = "metadata")` | All items are copied by name |
#'
# nolint end: line_length_linter
#'
#' ## The `reducedDims_mapping` argument
#'
#' For the simpler named vector format, the names should be the names of
#' `reducedDims` in the resulting `SingleCellExperiment` object, and the values
#' should be the keys of `obsm` in `adata`.
#'
#' For more advanced mapping, use the list format where each item is a vector
#' with the following names used to create a
#' [`SingleCellExperiment::LinearEmbeddingMatrix`] (if `featureLoadings` or
#' `metadata` is provided):
#'
#' - `sampleFactors`: a key of the `obsm` slot in `adata`,
#'   `adata$obsm[[sampleFactors]]` is passed to the `sampleFactors` argument
#' - `featureLoadings`: a key of the `varm` slot in `adata` (optional),
#'   `adata$varm[[featureLoadings]]` is passed to the `featureLoadings` argument
#' - `metadata`: a key of the `uns` slot in `adata` (optional),
#'   `adata$uns[[metadata]]` is passed to the `metadata` argument
#'
#' ## The `x_mapping` and `assays_mapping` arguments
#'
#' In order to specify where the data in `adata$X` will be stored in the
#' `assays(sce)` slot of the resulting object, you can use either the `x_mapping`
#' argument or the `assays_mapping` argument.
#' If you use `x_mapping`, it should be a string specifying the name of the layer
#' in `assays(sce)` where the data in `adata$X` will be stored.
#' If you use `assays_mapping`, it should be a named vector where names are names
#' of `assays(sce)` and values are keys of `layers` in `adata`.
#' In order to indicate the `adata$X` slot, you use `NA` as the value in the vector.
#' The name you provide for `x_mapping` may not be a name in `assays_mapping`.
#'
#' @return A `SingleCellExperiment` object containing the requested data from
#'   `adata`
#' @keywords internal
#'
#' @family object converters
#'
#' @examplesIf rlang::is_installed("SingleCellExperiment")
#'   ad <- AnnData(
#'     X = matrix(1:5, 3L, 5L),
#'     layers = list(A = matrix(5:1, 3L, 5L), B = matrix(letters[1:5], 3L, 5L)),
#'     obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'     var = data.frame(row.names = letters[1:5], gene = 1:5)
#'   )
#'
#'   # Default usage
#'   sce <- ad$as_SingleCellExperiment(
#'     assays_mapping = TRUE,
#'     colData_mapping = TRUE,
#'     rowData_mapping = TRUE,
#'     reducedDims_mapping = TRUE,
#'     colPairs_mapping = TRUE,
#'     rowPairs_mapping = TRUE,
#'     metadata_mapping = TRUE
#'  )
# nolint start: object_name_linter
as_SingleCellExperiment <- function(
  adata,
  x_mapping = NULL,
  assays_mapping = TRUE,
  colData_mapping = TRUE,
  rowData_mapping = TRUE,
  reducedDims_mapping = TRUE,
  colPairs_mapping = TRUE,
  rowPairs_mapping = TRUE,
  metadata_mapping = TRUE
) {
  # nolint end: object_name_linter
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
  assays_mapping <- get_mapping(
    assays_mapping,
    .as_SCE_guess_all,
    adata,
    "assays_mapping",
    slot = "layers"
  )
  colData_mapping <- get_mapping(
    colData_mapping,
    .as_SCE_guess_all,
    adata,
    "colData_mapping",
    slot = "obs"
  )
  rowData_mapping <- get_mapping(
    rowData_mapping,
    .as_SCE_guess_all,
    adata,
    "rowData_mapping",
    slot = "var"
  )
  reducedDims_mapping <- get_mapping(
    reducedDims_mapping,
    .as_SCE_guess_reducedDims,
    adata,
    "reducedDims_mapping"
  )
  colPairs_mapping <- get_mapping(
    colPairs_mapping,
    .as_SCE_guess_all,
    adata,
    "colPairs_mapping",
    slot = "obsp"
  )
  rowPairs_mapping <- get_mapping(
    rowPairs_mapping,
    .as_SCE_guess_all,
    adata,
    "rowPairs_mapping",
    slot = "varp"
  )
  metadata_mapping <- get_mapping(
    metadata_mapping,
    .as_SCE_guess_all,
    adata,
    "metadata_mapping",
    slot = "uns"
  )
  # nolint end object_name_linter

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done

  if (!is.null(x_mapping)) {
    assays_mapping <- setNames(
      c(NA, assays_mapping),
      c(x_mapping, names(assays_mapping))
    )
  } else if (!is.null(adata$X)) {
    assays_mapping[["X"]] <- NA
  }
  if (any(duplicated(names(assays_mapping)))) {
    cli_abort(
      "{.arg assays_mapping} or {.arg x_mapping} must not contain any duplicate names",
      "i" = "Found duplicate names: {.val {names(assays_mapping)[duplicated(names(assays_mapping))]}}"
    )
  }
  sce_assays <- vector("list", length(assays_mapping))
  names(sce_assays) <- names(assays_mapping)
  for (i in seq_along(assays_mapping)) {
    from <- assays_mapping[[i]]
    to <- names(assays_mapping)[[i]]
    sce_assays[[to]] <-
      if (is.na(from)) {
        to_R_matrix(adata$X)
      } else {
        to_R_matrix(adata$layers[[from]])
      }
  }

  # construct colData
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- .as_SCE_process_dataframe_mapping(adata, colData_mapping, "obs")

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- .as_SCE_process_dataframe_mapping(adata, rowData_mapping, "var")

  # construct colPairs
  # trackstatus: class=SingleCellExperiment, feature=get_obsp, status=done
  col_pairs <- .as_SCE_process_simple_mapping(adata, colPairs_mapping, "obsp")

  # construct rowPairs
  # trackstatus: class=SingleCellExperiment, feature=get_varp, status=done
  row_pairs <- .as_SCE_process_simple_mapping(adata, rowPairs_mapping, "varp")

  # construct metadata
  # trackstatus: class=SingleCellExperiment, feature=get_uns, status=done
  metadata <- .as_SCE_process_simple_mapping(adata, metadata_mapping, "uns")

  # construct reducedDims
  reduceddims <- .as_SCE_process_reducedDims_mapping(adata, reducedDims_mapping)

  # construct output object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = sce_assays,
    colData = col_data,
    rowData = row_data,
    colPairs = col_pairs,
    rowPairs = row_pairs,
    metadata = metadata,
    checkDimnames = TRUE
  )

  # only here to ensure that the dimensions are right
  SingleCellExperiment::reducedDims(sce) <- reduceddims

  sce
}

#' Convert an AnnData object to a SingleCellExperiment object
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated, use `adata$as_SingleCellExperiment()` instead
#'
#' @param ... Arguments passed to [as_SingleCellExperiment()]
#'
#' @return A `SingleCellExperiment` object
#' @export
# nolint start: object_name_linter
to_SingleCellExperiment <- function(...) {
  # nolint end: object_name_linter
  lifecycle::deprecate_warn(
    when = "0.99.0",
    what = "to_SingleCellExperiment()",
    with = "adata$as_SingleCellExperiment()"
  )
  as_SingleCellExperiment(...)
}

# nolint start: object_length_linter object_name_linter
.as_SCE_guess_all <- function(adata, slot) {
  # nolint end: object_length_linter object_name_linter
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  self_name(names(adata[[slot]]))
}

# nolint start: object_length_linter object_name_linter
.as_SCE_guess_reducedDims <- function(adata) {
  # nolint end: object_length_linter object_name_linter
  purrr::map(adata$obsm_keys(), function(.obsm) {
    if (!is.numeric(as.matrix(adata$obsm[[.obsm]]))) {
      return(NULL)
    }

    mapping <- c(sampleFactors = .obsm)
    if (.obsm == "X_pca" && "PCs" %in% names(adata$varm)) {
      mapping["featureLoadings"] <- "PCs"
    }

    mapping
  }) |>
    setNames(adata$obsm_keys()) |>
    purrr::compact()
}

# nolint start: object_length_linter object_name_linter
.as_SCE_process_simple_mapping <- function(adata, mapping, slot) {
  # nolint end: object_length_linter object_name_linter

  if (rlang::is_empty(mapping)) {
    return(list())
  }

  purrr::map(mapping, function(.item) {
    adata[[slot]][[.item]]
  })
}

# nolint start: object_length_linter object_name_linter
.as_SCE_process_dataframe_mapping <- function(adata, mapping, slot) {
  # nolint end: object_length_linter object_name_linter
  if (slot == "obs") {
    row_names <- adata$obs_names
  } else if (slot == "var") {
    row_names <- adata$var_names
  } else {
    cli_abort(c(
      "{.arg slot} must be either {.val obs} or {.val var}, but is {.val {slot}}"
    ))
  }

  if (rlang::is_empty(mapping)) {
    return(S4Vectors::DataFrame(row.names = row_names))
  }

  purrr::map(mapping, function(.item) {
    adata[[slot]][[.item]]
  }) |>
    S4Vectors::DataFrame(row.names = row_names)
}

# trackstatus: class=SingleCellExperiment, feature=get_obsm, status=done
# trackstatus: class=SingleCellExperiment, feature=get_varm, status=done
# nolint start: object_length_linter object_name_linter
.as_SCE_process_reducedDim <- function(
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
.as_SCE_process_reducedDims_mapping <- function(adata, reducedDims_mapping) {
  # nolint end: object_length_linter object_name_linter

  # If the mapping is a vector expand it to the list format
  if (is.atomic(reducedDims_mapping)) {
    # nolint start: object_name_linter
    reducedDims_mapping <- purrr::map(reducedDims_mapping, function(.obsm) {
      # nolint end: object_name_linter
      c(sampleFactors = .obsm)
    })
  } else {
    purrr::walk(names(reducedDims_mapping), function(.name) {
      mapping <- reducedDims_mapping[[.name]]
      if (!(is.character(mapping) && ("sampleFactors" %in% names(mapping)))) {
        cli_abort(c(
          paste(
            "If it is a list, each item in {.arg reducedDims_mapping} must be a
            {.cls character} vector with the name {.val sampleFactors} and
            optional names {.val featureLoadings} and {.val metadata}"
          ),
          "i" = paste(
            "{.code names(reducedDims_mapping[[{.val { .name }}]])}:",
            "{style_vec(names(mapping))}"
          )
        ))
      }
    })
  }

  purrr::map(reducedDims_mapping, function(.map) {
    varm <- if ("featureLoadings" %in% names(.map)) {
      .map[["featureLoadings"]]
    } else {
      NULL
    }
    uns <- if ("metadata" %in% names(.map)) .map[["metadata"]] else NULL

    .as_SCE_process_reducedDim(
      adata,
      obsm_key = .map[["sampleFactors"]],
      varm_key = varm,
      uns_key = uns
    )
  })
}
