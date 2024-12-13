#' Convert an AnnData object to a SingleCellExperiment object
#'
#' `to_SingleCellExperiment()` converts an AnnData object
#'   to a SingleCellExperiment.
#'
#' @param object an AnnData object, e.g., InMemoryAnnData
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
to_SingleCellExperiment <- function(object) { # nolint
  check_requires("Converting AnnData to SingleCellExperiment", "SingleCellExperiment", "Bioc")

  stopifnot(
    inherits(object, "AbstractAnnData")
  )

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  ## FIXME: name of 'X' from metadata[["X_name"]]
  x_name <- "X"
  assay_names <- as.character(c(
    if (!is.null(object$X)) x_name,
    object$layers_keys()
  ))

  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  assays <- vector("list", length(assay_names))
  names(assays) <- assay_names
  for (assay in assay_names) {
    value <- if (identical(assay, x_name)) {
      object$X
    } else {
      object$layers[[assay]]
    }
    ## FIXME: is transposition robust & efficient here?
    assays[[assay]] <- t(value)
  }

  # construct colData
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- S4Vectors::DataFrame(
    object$obs,
    row.names = object$obs_names
  )

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- S4Vectors::DataFrame(
    object$var,
    row.names = object$var_names
  )

  # construct output object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = col_data,
    rowData = row_data,
    metadata = list(),
    ## FIXME: assign object$uns to metadata
    checkDimnames = TRUE
  )

  ## reducedDims

  ## rowPairs

  sce
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

        adata$uns[[uns_name]] <- metadata(sce)[[uns]]
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
