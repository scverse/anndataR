#' Read Zarr
#'
#' Read data from a Zarr store
#'
#' @param path Path to the Zarr store to read
#' @param to The type of object to return. Must be one of: "InMemoryAnnData",
#'   "HDF5AnnData", "SingleCellExperiment", "Seurat"
#' @param ... Extra arguments provided to `adata$to_SingleCellExperiment()` or
#'   `adata$to_Seurat()`. See [AnnData()] for more information on the arguments of
#'   these functions. Note: update this documentation when
#'   [`r-lib/roxygen2#955`](https://github.com/r-lib/roxygen2/issues/955) is resolved.
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
    to = c("InMemoryAnnData", "HDF5AnnData", "SingleCellExperiment", "Seurat", "ZarrAnnData"),
    ...) {
  to <- match.arg(to)

  adata <- ZarrAnnData$new(path)

  fun <- switch(to,
    "SingleCellExperiment" = to_SingleCellExperiment,
    "Seurat" = to_Seurat,
    "InMemoryAnnData" = to_InMemoryAnnData,
    "HDF5AnnData" = to_HDF5AnnData,
    "ZarrAnnData" = NULL
  )

  if (!is.null(fun)) {
    fun(adata, ...)
  } else {
    adata
  }
}
