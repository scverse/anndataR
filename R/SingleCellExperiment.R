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
# TODO: fix snake_case + CamelCase
to_SingleCellExperiment <- function(object) { # nolint
  stopifnot(
    inherits(object, "AbstractAnnData")
  )

  ## mostly following zellkonverter:::.native_reader
  assay <- object$layer
  x <- object$X
  if (!is.null(x)) {
    ## FIXME: name of 'X' from metadata[["X_name"]]
    assay <- c(list(X = x), assay)
  }
  ## FIXME: better transposition -- if sparse, then always dgCMatrix
  assay <- lapply(assay, t)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assay,
    colData = S4Vectors::DataFrame(
      object$obs,
      row.names = object$obs_names
    ),
    rowData = S4Vectors::DataFrame(
      object$var,
      row.names = object$var_names
    ),
    metadata = list(),
    ## FIXME: assign object$uns to metadata
    checkDimnames = TRUE
  )

  ## reducedDims

  ## rowPairs

  sce
}
