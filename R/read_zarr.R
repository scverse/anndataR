#' Read Zarr
#'
#' Read data from a Zarr store
#'
#' @param path Path to the H5AD file to read
#' @param to The type of object to return. Must be one of: "InMemoryAnnData",
#'   "HDF5AnnData", "SingleCellExperiment", "Seurat"
#' @param ... Extra arguments provided to [to_SingleCellExperiment()] or
#'   [to_Seurat()]
#'
#' @return The object specified by `to`
#' @export
#'
#' @examples
#' h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
#'
#' # Read the H5AD as a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- read_zarr(h5ad_file, to = "SingleCellExperiment")
#' }
#'
#' # Read the H5AD as a Seurat object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seurat <- read_zarr(h5ad_file, to = "Seurat")
#' }
read_zarr <- function(
    path,
    to = c("InMemoryAnnData", "ZarrAnnData", "SingleCellExperiment", "Seurat"),
    ...) {
  to <- match.arg(to)

  adata <- ZarrAnnData$new(path)

  fun <- switch(to,
    "SingleCellExperiment" = to_SingleCellExperiment,
    "Seurat" = to_Seurat,
    "InMemoryAnnData" = to_InMemoryAnnData,
    "ZarrAnnData" = NULL
  )

  if (!is.null(fun)) {
    fun(adata, ...)
  } else {
    adata
  }
}
