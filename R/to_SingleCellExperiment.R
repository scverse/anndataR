# nolint start: line_length_linter
# Sources used to understand how to convert between SingleCellExperiment and AnnData
# https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#2_Creating_SingleCellExperiment_instances
# https://bioconductor.org/packages/3.20/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#column-sample-data
# https://github.com/ivirshup/sc-interchange/issues
# https://github.com/ivirshup/sc-interchange/issues/2
# https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#3_Adding_low-dimensional_representations
# nolint end: line_length_linter

#' Convert an AnnData object to a SingleCellExperiment object
#'
#' `to_SingleCellExperiment()` converts an AnnData object
#'   to a SingleCellExperiment object.
#'
#' @param adata an AnnData object, e.g., InMemoryAnnData
#'
#' @param assays_mapping A named list mapping `layers` in `adata` to
#' `assay` names in the created SingleCellExperiment object.
#' The names of the list should be the names of the `assays` in the
#' resulting SingleCellExperiment object, and the values should be the names
#' of the `layers` in `adata`, and should include the `X` matrix as well.
#' If `X` is not in the list, it will be added as `counts` or `data`.
#' @param colData_mapping a named list mapping `obs` in `adata` to
#' `colData` in the created SingleCellExperiment object.
#' The names of the list should be the names of the `colData` columns in the
#' resulting SingleCellExperiment object. The values of the list should be the
#' names of the `obs` columns in `adata`.
#' @param rowData_mapping a named list mapping `var` names in `adata` to
#' `rowData` in the created SingleCellExperiment object. The names of the list
#' should be the names of the `rowData` columns in the resulting SingleCellExperiment
#' object. The values of the list should be the names of the `var` columns in `adata`.
#' @param reduction_mapping a named list mapping reduction names in `adata` to
#' reduction names in the created SingleCellExperiment object.
#' The names of the list should be the names of the `reducedDims` in the resulting
#' SingleCellExperiment object.
#' The values of the list should also be a named list with the following keys:
#' - `obsm`: the name of the `obsm` slot in `adata`
#' - `varm`: the name of the `varm` slot in `adata`
#' - `uns`: the name of the `uns` slot in `adata`
#' @param colPairs_mapping a named list mapping obsp names in `adata` to
#' colPairs names in the created SingleCellExperiment object.
#' The names of the list should be the names of the `colPairs` in the resulting
#' SingleCellExperiment object. The values of the list should be the names of the
#' `obsp` in `adata`.
#' @param rowPairs_mapping a named list mapping varp names in `adata` to
#' rowPairs names in the created SingleCellExperiment object.
#' The names of the list should be the names of the `rowPairs` in the resulting
#' SingleCellExperiment object. The values of the list should be the names of the
#' `varp` in `adata`.
#' @param metadata_mapping a named list mapping uns names in `adata` to
#' metadata names in the created SingleCellExperiment object.
#' The names of the list should be the names of the `metadata` in the resulting
#' SingleCellExperiment object. The values of the list should be the names of the
#' `uns` in `adata`.
#'
#' @return `to_SingleCellExperiment()` returns a SingleCellExperiment
#'   representing the content of `adata`.
#'
#' @examples
#' if (interactive()) {
#'   ## useful when interacting with the SingleCellExperiment !
#'   library(SingleCellExperiment)
#' }
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
#' ## construct a SingleCellExperiment from an AnnData object
#' sce <- to_SingleCellExperiment(ad)
#' sce
#' @export
# nolint start: cyclocomp_linter
to_SingleCellExperiment <- function(
  # nolint end: cyclocomp_linter
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
  if (is.null(assays_mapping)) {
    assays_mapping <- to_SCE_guess_assays(adata)
  }
  if (is.null(colData_mapping)) {
    colData_mapping <- to_SCE_guess_all(adata, "obs") # nolint
  }
  if (is.null(rowData_mapping)) {
    rowData_mapping <- to_SCE_guess_all(adata, "var") # nolint
  }
  if (is.null(reduction_mapping)) {
    reduction_mapping <- to_SCE_guess_reduction(adata)
  }
  if (is.null(colPairs_mapping)) {
    colPairs_mapping <- to_SCE_guess_all(adata, "obsp") # nolint
  }
  if (is.null(rowPairs_mapping)) {
    rowPairs_mapping <- to_SCE_guess_all(adata, "varp") # nolint
  }
  if (is.null(metadata_mapping)) {
    metadata_mapping <- to_SCE_guess_all(adata, "uns")
  }

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  sce_assays <- vector("list", length(assays_mapping))
  names(sce_assays) <- names(assays_mapping)
  for (i in seq_along(assays_mapping)) {
    from <- assays_mapping[[i]]
    to <- names(assays_mapping)[[i]]
    if (from != "X") {
      sce_assays[[to]] <- t(adata$layers[[from]])
    } else {
      sce_assays[[to]] <- t(adata$X)
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

  # construct reducedDims
  # trackstatus: class=SingleCellExperiment, feature=get_reductions, status=wip
  reduceddims <- vector("list", length(reduction_mapping))
  names(reduceddims) <- names(reduction_mapping)
  for (i in seq_along(reduction_mapping)) {
    name <- names(reduction_mapping)[[i]]
    reduction <- reduction_mapping[[i]]

    obsm_key <- reduction$obsm
    varm_key <- reduction$varm
    uns_key <- reduction$uns

    reduceddims[[name]] <- .to_SingleCellExperiment_process_reduction(
      adata,
      name,
      obsm_key,
      varm_key,
      uns_key
    )
  }

  # construct colPairs
  # trackstatus: class=SingleCellExperiment, feature=get_obsp, status=done
  col_pairs <- .to_SCE_process_simple_mapping(adata, colPairs_mapping, "obsp")

  # construct rowPairs
  # trackstatus: class=SingleCellExperiment, feature=get_varp, status=done
  row_pairs <- .to_SCE_process_simple_mapping(adata, rowPairs_mapping, "varp")

  # construct metadata
  # trackstatus: class=SingleCellExperiment, feature=get_uns, status=done
  metadata <- .to_SCE_process_simple_mapping(adata, metadata_mapping, "uns")

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
to_SCE_guess_assays <- function(adata) {
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
to_SCE_guess_all <- function(adata, slot) {
  # nolint end: object_length_linter object_name_linter
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  mapping <- names(adata[[slot]])
  names(mapping) <- names(adata[[slot]])

  mapping
}

# nolint start: object_length_linter object_name_linter
to_SCE_guess_reduction <- function(adata) {
  # nolint end: object_length_linter object_name_linter
  .to_Seurat_guess_reductions(adata) # nolint object_usage_linter
}

# nolint start: object_length_linter object_name_linter
.to_SCE_process_simple_mapping <- function(adata, mapping, slot) {
  # nolint end: object_length_linter object_name_linter
  # check if mapping contains all columns of slot
  if (length(setdiff(names(adata[[slot]]), names(mapping))) == 0) {
    adata[[slot]]
  } else {
    mapped <- lapply(seq_along(mapping), function(i) {
      adata[[slot]][[mapping[[i]]]]
    })
    names(mapped) <- names(mapping)
    mapped
  }
}

# nolint start: object_length_linter object_name_linter
.to_SingleCellExperiment_process_reduction <- function(
  # nolint end: object_length_linter object_name_linter
  adata,
  key,
  obsm_key,
  varm_key,
  uns_key
) {
  # nolint
  embedding <- adata$obsm[[obsm_key]]

  if (is.null(embedding)) {
    cli_abort(
      c(
        "{.val {obsm_key}} is not an item in {.code adata$obsm}",
        "i" = "{.code adata$obsm_keys()}: {.val {adata$obsm_keys()}}"
      )
    )
  }

  rownames(embedding) <- adata$obs_names

  if (!is.null(varm_key) && varm_key %in% names(adata$varm)) {
    loadings <- adata$varm[[varm_key]]
    rownames(loadings) <- colnames(embedding)

    metadata <- list()
    if (!is.null(uns_key) && uns_key %in% names(adata$uns)) {
      metadata <- adata$uns[[uns_key]]
    }

    lem <- SingleCellExperiment::LinearEmbeddingMatrix(
      sampleFactors = embedding,
      featureLoadings = loadings,
      metadata = metadata
    )
    rownames(lem) <- rownames(embedding)
    colnames(lem) <- colnames(loadings)
    lem
  } else {
    embedding
  }
}
