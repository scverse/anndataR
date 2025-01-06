#' Convert an AnnData object to a SingleCellExperiment object
#'
#' `to_SingleCellExperiment()` converts an AnnData object
#'   to a SingleCellExperiment object.
#'
#' @param adata an AnnData object, e.g., InMemoryAnnData
#' 
#' @param assays_mapping A named list mapping `layers` in `adata` to
#'   `assay` names in the created SingleCellExperiment object.
#' @param col_data_mapping a named list mapping `obs` in `adata` to
#'  `colData` in the created SingleCellExperiment object
#' @param row_data_mapping a named list mapping `var` names in `adata` to
#'   `rowData` in the created SingleCellExperiment object
#' @param reduction_mapping a named list mapping reduction names in `adata` to
#' reduction names in the created SingleCellExperiment object
#' @param colPairs_mapping a named list mapping obsp names in `adata` to
#' colPairs names in the created SingleCellExperiment object
#' @param rowPairs_mapping a named list mapping varp names in `adata` to
#' rowPairs names in the created SingleCellExperiment object
#' @param metadata_mapping a named list mapping uns names in `adata` to
#' metadata names in the created SingleCellExperiment object
#'
#' @return `to_SingleCellExperiment()` returns a SingleCellExperiment
#'   representing the content of `object`.
#'
#' @noRd
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
to_SingleCellExperiment <- function(
  adata,
  assays_mapping = NULL,
  col_data_mapping = NULL,
  row_data_mapping = NULL,
  reduction_mapping = NULL,
  colPairs_mapping = NULL,
  rowPairs_mapping = NULL,
  metadata_mapping = NULL) {

  check_requires("Converting AnnData to SingleCellExperiment", "SingleCellExperiment", "Bioc")

  stopifnot(inherits(adata, "AbstractAnnData"))

  # guess mappings if not provided
  if (is.null(assays_mapping)) {
    assays_mapping <- to_SCE_guess_assays(adata)
  }
  if (is.null(col_data_mapping)) {
    col_data_mapping <- to_SCE_guess_all(adata, "obs")
  }
  if(is.null(row_data_mapping)) {
    row_data_mapping <- to_SCE_guess_all(adata, "var")
  }
  if (is.null(reduction_mapping)) {
    reduction_mapping <- to_SCE_guess_reduction(adata)
  }
  if (is.null(colPairs_mapping)) {
    colPairs_mapping <- to_SCE_guess_all(adata, "obsp")
  }
  if (is.null(rowPairs_mapping)) {
    rowPairs_mapping <- to_SCE_guess_all(adata, "varp")
  }
  if (is.null(metadata_mapping)) {
    metadata_mapping <- to_SCE_guess_all(adata, "uns")
  }

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  sce_assays <- vector("list", length(assays_mapping))
  names(sce_assays) <- assays_mapping
  for (i in seq_along(assays_mapping)){
    to <- assays_mapping[[i]]
    from <- names(assays_mapping)[[i]]
    if (from != "X"){
      sce_assays[[to]] <- t(adata$layers[[from]])
    } else {
      sce_assays[[to]] <- t(adata$X)
    }
  }

  # construct colData
  # FIXME: probably better way to make a dataframe from a list of vectors
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- .to_SCE_process_simple_mapping(adata, col_data_mapping, "obs")
  col_data <- as(col_data, "DataFrame")

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- .to_SCE_process_simple_mapping(adata, row_data_mapping, "var")
  row_data <- as(row_data, "DataFrame")

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

    reduceddims[[name]] <- .from_SingleCellExperiment_process_reduction(adata, name, obsm_key, varm_key, uns_key)
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

  # construct output object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = sce_assays,
    colData = col_data,
    rowData = row_data,
    reducedDims = reduceddims,
    colPairs = col_pairs,
    rowPairs = row_pairs,
    metadata = metadata,
    checkDimnames = TRUE
  )

  sce
}

to_SCE_guess_assays <- function(adata) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  layers <- list()

  if (!is.null(adata$X)) {
    layer_name_for_x <- 
      if (!"counts" %in% names(adata$layers)) { # could expand checks, to check if integers
        "counts"
      } else {
        "data"
      }
    layers[["X"]] <- layer_name_for_x
  }

  for (layer_name in names(adata$layers)) {
    layers[[layer_name]] <- layer_name
  }

  layers
}

to_SCE_guess_all <- function(adata, slot){
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  mapping <- names(adata[[slot]])
  names(mapping) <- names(adata[[slot]])

  mapping
}

to_SCE_guess_reduction <- function(adata){
  list()
}

.to_SCE_process_simple_mapping <- function(adata, mapping, slot){
  # check if mapping contains all columns of slot
  if(length(setdiff(names(adata[[slot]]), names(mapping))) == 0){
    adata[[slot]]
  } else {
    mapped <- lapply(seq_along(mapping), function(i) {
      adata[[slot]][[mapping[[i]]]]
    })
    names(mapped) <- names(mapping)
  }
}

.from_SingleCellExperiment_process_reduction <- function(adata, key, obsm_key, varm_key, uns_key) {
  # check arguments

  embedding <- adata$obsm[[obsm_key]]

  if (is.null(embedding)) {
    stop(paste0("No embedding found for key '", obsm_key, "' in adata$obsm"))
  }

  rownames(embedding) <- adata$obs_names

  if (varm_key %in% names(adata$varm)) {
    loadings <- adata$varm[[varm_key]]
    rownames(loadings) <- colnames(embedding)
    
    metadata <- list()
    if (uns_key %in% names(adata$uns)) {
      metadata <- adata$uns[[uns_key]]
    }

    LinearEmbeddingMatrix(
      sampleFactors = embedding,
      featureLoadings = loadings,
      metadata = metadata
    )
  } else {
    embedding
  }
}

#' Convert a SingleCellExperiment object to an AnnData object
#'
#' `from_SingleCellExperiment()` converts a
#'   SingleCellExperiment to an AnnData object.
#'
#' @param sce An object inheriting from SingleCellExperiment.
#'
#' @param output_class Name of the AnnData class. Must be one of `"HDF5AnnData"`
#' or `"InMemoryAnnData"`.
#'
#' @param ... Additional arguments passed to the generator function.
#' See the "Details" section for more information on which parameters
#'
#' @return `from_SingleCellExperiment()` returns an AnnData object
#'   (e.g., InMemoryAnnData) representing the content of `sce`.
#'
#' @examples
#' ## construct an AnnData object from a SingleCellExperiment
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(
#'   assays = list(counts = matrix(1:5, 5L, 3L)),
#'   colData = DataFrame(cell = 1:3),
#'   rowData = DataFrame(gene = 1:5)
#' )
#' from_SingleCellExperiment(sce, "InMemory")
#'
#' @export
# nolint start: object_name_linter
from_SingleCellExperiment <- function(
    # nolint end: object_name_linter
    sce,
    output_class = c("InMemory", "HDF5AnnData"),
    x_mapping = NULL,
    layers_mapping = NULL,
    obsm_mapping = NULL,
    varm_mapping = NULL,
    obsp_mapping = NULL,
    varp_mapping = NULL,
    uns_mapping = NULL,
    ...) {

  check_requires("Converting SingleCellExperiment to AnnData", "SingleCellExperiment", "Bioc")

  output_class <- match.arg(output_class)

  stopifnot(inherits(sce, "SingleCellExperiment"))

  # TODO: guess mappings if not provided

  # trackstatus: class=SingleCellExperiment, feature=set_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=set_obs_names, status=done
  obs <- as.data.frame(
    SummarizedExperiment::colData(sce)
  )

  # trackstatus: class=SingleCellExperiment, feature=set_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=set_var_names, status=done
  var <- as.data.frame(
    SummarizedExperiment::rowData(sce)
  )

  # fetch generator
  generator <- get_anndata_constructor(output_class)

  tryCatch(
    {
      adata <- generator$new(
        obs = obs,
        var = var,
        ...
      )

      # fetch X
      # trackstatus: class=SingleCellExperiment, feature=set_X, status=wip
      x <- NULL
      if(!is.null(x_mapping)){
        x <- .to_anndata_matrix(SummarizedExperiment::assay(sce, withDimnames = FALSE))
      }

      # fetch layers
      # trackstatus: class=SingleCellExperiment, feature=set_layers, status=wip
      for (i in seq_along(layers_mapping)) {
        layer <- layers_mapping[[i]]
        layer_name <- names(layers_mapping)[[i]]

        adata$uns[[layer_name]] <- .to_anndata_matrix(assay(sce, layer))
      }

      # trackstatus: class=SingleCellExperiment, feature=set_obsm, status=wip
      for (i in seq_along(obsm_mapping)) {
        obsm <- obsm_mapping[[i]]
        obsm_name <- names(obsm_mapping)[[i]]

        if (!is.character(obsm) || length(obsm) != 2){
          stop("each obsm_mapping must be a character vector of length 2")
        }

        obsm_slot <- obsm[[1]]
        obsm_key <- obsm[[2]]

        if (obsm_slot == "reducedDims") {
          adata$obsm[[obsm_name]] <- reducedDim(sce, obsm_key)
        } else if (obsm_slot == "metadata") {
          adata$obsm[[obsm_name]] <- metadata(sce)[[obsm_key]]
        } 
      }

      # have a look at this for loadings: https://github.com/ivirshup/sc-interchange/issues/2
      # could be in rowData or metadata, or in LinearEmbeddingMatrix (which is in a reducedDims slot)?
      # if in rowData --> need a prefix? or a subslot?
      # if in metadata --> need a subslot

      # fetch varm
      # trackstatus: class=SingleCellExperiment, feature=set_varm, status=todo
      for (i in seq_along(varm_mapping)) {
        varm <- varm_mapping[[i]]
        varm_name <- names(varm_mapping)[[i]]

      }

      # fetch obsp
      # trackstatus: class=SingleCellExperiment, feature=set_obsp, status=wip
      for (i in seq_along(obsp_mapping)) {
        obsp <- obsp_mapping[[i]]
        obsp_name <- names(obsp_mapping)[[i]]

        adata$obsp[[obsp_name]] <- .to_anndata_matrix(colPair(sce, obsp))
      }

      # fetch varp
      # trackstatus: class=SingleCellExperiment, feature=set_varp, status=wip
      for (i in seq_along(varp_mapping)) {
        varp <- varp_mapping[[i]]
        varp_name <- names(varp_mapping)[[i]]

        adata$varp[[varp_name]] <- .to_anndata_matrix(rowPair(sce, varp))
      }

      # fetch uns
      # trackstatus: class=SingleCellExperiment, feature=set_uns, status=wip
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

    },
    error = function(e) {
      if (output_class == "HDF5AnnData") {
        on.exit(cleanup_HDF5AnnData(adata))
      }
      stop(e)
    }
  )

}


.to_anndata_matrix <- function(dge) {
  m <- t(mat)
  if (inherits(m, "denseMatrix")) {
    m <- as.matrix(m)
  }
  m
}
