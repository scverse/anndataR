#' Read H5AD
#'
#' Read data from a H5AD file
#'
#' @param path Path to the H5AD file to read
#' @param to The type of object to return. Must be one of: "InMemoryAnnData",
#'   "HDF5AnnData", "SingleCellExperiment", "Seurat"
#' @param mode The mode to open the HDF5 file.
#'
#'   * `a` creates a new file or opens an existing one for read/write.
#'   * `r` opens an existing file for reading.
#'   * `r+` opens an existing file for read/write.
#'   * `w` creates a file, truncating any existing ones.
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists.
#'
#' @param ... Extra arguments provided to `adata$to_SingleCellExperiment()` or
#'   `adata$to_Seurat()`. See [AnnData()] for more information on the arguments of
#'   these functions. Note: update this documentation when
#'   [`r-lib/roxygen2#955`](https://github.com/r-lib/roxygen2/issues/955) is resolved.
#'
#' @return The object specified by `to`
#' @export
#'
#' @examples
#' h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
#'
#' # Read the H5AD as a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- read_h5ad(h5ad_file, to = "SingleCellExperiment")
#' }
#'
#' # Read the H5AD as a Seurat object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seurat <- read_h5ad(h5ad_file, to = "Seurat")
#' }
read_h5ad <- function(
    path,
    to = c("InMemoryAnnData", "HDF5AnnData", "SingleCellExperiment", "Seurat"),
    mode = c("r", "r+", "a", "w", "w-", "x"),
    ...) {
  to <- match.arg(to)
  mode <- match.arg(mode)

  adata <- HDF5AnnData$new(path, mode = mode)

  fun <- switch(to,
    "SingleCellExperiment" = to_SingleCellExperiment,
    "Seurat" = to_Seurat,
    "InMemoryAnnData" = to_InMemoryAnnData,
    "HDF5AnnData" = NULL
  )

  if (!is.null(fun)) {
    fun(adata, ...)
  } else {
    adata
  }
}
