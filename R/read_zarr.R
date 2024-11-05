#' Read Zarr
#'
#' Read data from a Zarr store
#'
#' @param path Path to the Zarr store to read
#' @param to The type of object to return. Must be one of: "InMemoryAnnData",
#'   "HDF5AnnData", "SingleCellExperiment", "Seurat"
#' @param ... Extra arguments provided to [to_SingleCellExperiment()] or
#'   [to_Seurat()]
#'
#' @return The object specified by `to`
#' @export
#'
#' @examples
#' file <- system.file("extdata", "example.zarr", package = "anndataR")
#' store <- pizzarr::DirectoryStore$new(file)
#'
#' # Read the Zarr store as a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- read_zarr(store, to = "SingleCellExperiment")
#' }
#'
#' # Read the Zarr store as a Seurat object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seurat <- read_zarr(store, to = "Seurat")
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
