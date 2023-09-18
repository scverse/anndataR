#' Read H5AD
#'
#' Read data from a H5AD file
#'
#' @param path Path to the H5AD file to read
#' @param to The type of object to return. Must be one of: "SingleCellExperiment",
#'   "Seurat", "HDF5AnnData", "InMemoryAnnData"
#'
#' @return The object specified by `to`
#' @export
#'
#' @examples
#' h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
#' # Read the H5AD as a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- read_h5ad(h5ad_file, to = "SingleCellExperiment")
#' }
#' # Read the H5AD as a Seurat object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seurat <- read_h5ad(h5ad_file, to = "Seurat")
#' }
read_h5ad <- function(path, to = c("SingleCellExperiment", "Seurat", "HDF5AnnData", "InMemoryAnnData")) {
  to <- match.arg(to)

  adata <- HDF5AnnData$new(path)

  switch(to,
    "SingleCellExperiment" = to_SingleCellExperiment(adata),
    "Seurat" = to_Seurat(adata),
    "HDF5AnnData" = adata,
    "InMemoryAnnData" = to_InMemoryAnnData(adata)
  )
}