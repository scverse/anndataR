#' @rdname SingleCellExperiment
#'
#' @title Convert AnnData to SingleCellExperiment
#'
#' @description `AnnData_to_SingleCellExperiment()` converts an
#'     AnnData object to a SingleCellExperiment.
#'
#' @param object an AnnData object, e.g., InMemoryAnnData
#'
#' @return `AnnData_to_SingleCellExperiment()` returns a
#'     SingleCellExperiment representing the content of `object`.
#'
#' @examples
#' ad <- InMemoryAnnData$new(
#'     X = matrix(1:5, 3L, 5L),
#'     obs = data.frame(cell = 1:3, row.names = LETTERS[1:3]),
#'     var = data.frame(gene = 1:5, row.names = letters[1:5])
#' )
#' to_SingleCellExperiment(ad)
#'
#' @export
to_SingleCellExperiment <- function(object) {
    stopifnot(
        inherits(object, "AbstractAnnData")
    )

    ## mostly following zellkonverter:::.native_reader
    assay <- object$layer
    X <- object$X
    if (!is.null(X))
        ## FIXME: name of 'X' from metadata[["X_name"]]
        assay <- c(list(X = X), assay)
    ## FIXME: better transposition -- if sparse, then always dgCMatrix
    assay <- lapply(assay, t)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = assay,
        metadata = list(),
        ## FIXME: metadata = object$uns
        checkDimnames = TRUE
    )

    rowData <- object$var
    if (!is.null(rowData)) {
        SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(rowData)
        rownames(sce) <- rownames(rowData)
    }

    colData <- object$obs
    if (!is.null(colData)) {
        SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(colData)
        colnames(sce) <- rownames(colData)
    }

    ## reducedDims

    ## rowPairs

    sce
}
