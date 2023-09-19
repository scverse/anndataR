#' Convert: \code{list} ==> \code{DimReducObject}
#'
#' Create a \pkg{Seurat} \code{DimReducObject}
#' from embeddings and loadings. 
#' @param obsm Sample embeddings.
#' @param varm Feature loadings.
#' @param stdev Standard deviation (if applicable) for the
#'  dimensional reduction.
#' @param assay Assay to use.
#' @param key Key to use (name of embedding).
#' @inheritDotParams SeuratObject::CreateDimReducObject
#'
#' @export
#' @importFrom SeuratObject CreateDimReducObject
to_DimReduc <- function(obsm,
                        varm,
                        stdev=NULL,
                        assay,
                        key,
                        ...) {
  
  Seurat::CreateDimReducObject(
    embeddings = obsm,
    loadings = varm,
    stdev = if (is.null(stdev)) {
      numeric()
    } else {
      as.numeric(stdev)
    },
    assay = assay,
    key = key,
    ...)
}
