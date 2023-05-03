#' @rdname SingleCellExperiment
#'
#' @title Convert AnnData to SingleCellExperiment
#'
#' @description `to_SingleCellExperiment()` converts an AnnData object
#'     to a SingleCellExperiment.
#'
#' @param object an AnnData object, e.g., InMemoryAnnData
#'
#' @return `to_SingleCellExperiment()` returns a SingleCellExperiment
#'     representing the content of `object`.
#'
#' @examples
#' if (interactive()) {
#'   ## useful when interacting with the SingleCellExperiment !
#'   library(SingleCellExperiment)
#' }
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5),
#'   obs_names = LETTERS[1:3],
#'   var_names = letters[1:5]
#' )
#' to_SingleCellExperiment(ad)
#'
#' @export
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
    assays[[assay]] <- Matrix::t(value)
  }

  # construct colData
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  colData <- S4Vectors::DataFrame(
    object$obs, row.names = object$obs_names
  )

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  rowData <- S4Vectors::DataFrame(
    object$var, row.names = object$var_names
  )

  # construct output object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    
    colData = colData,
    rowData = rowData,
    metadata = list(),
    ## FIXME: assign object$uns to metadata
    checkDimnames = TRUE
  )

  ## reducedDims

  ## rowPairs

  sce
}
