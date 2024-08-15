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
    ...) {
  stopifnot(
    inherits(sce, "SingleCellExperiment")
  )

  output_class <- match.arg(output_class)

  # fetch generator
  generator <- get_anndata_constructor(output_class)

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

  # trackstatus: class=SingleCellExperiment, feature=set_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=set_layers, status=done
  x_and_layers <- lapply(
    SummarizedExperiment::assays(sce, withDimnames = FALSE),
    function(mat) {
      m <- t(mat)
      # nolint start
      # WORKAROUND: convert denseMatrix to matrix, because otherwise:
      # - Could not write element '/layers/integer_dense' of type 'dgeMatrix':
      #   no applicable method for 'h5writeDataset' applied to an object of class "c('dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'compMa
      # - Could not write element '/layers/integer_dense_with_nas' of type 'dgeMatrix':
      #   no applicable method for 'h5writeDataset' applied to an object of class "c('dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'compMatrix', 'Matrix', 'replValueSp')"
      # nolint end
      if (inherits(m, "denseMatrix")) {
        m <- as.matrix(m)
      }
      m
    }
  )
  if (length(x_and_layers) == 0L) {
    x <- NULL
    layers <- list()
    names(layers) <- character()
  } else {
    x <- x_and_layers[[1]]
    layers <- x_and_layers[-1]
  }

  generator$new(
    X = x,
    obs = obs,
    var = var,
    layers = layers,
    ...
  )
}
