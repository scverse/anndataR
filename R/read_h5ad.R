#' Read H5AD
#'
#' Read data from a H5AD file
#'
#' @param path Path to the H5AD file to read
#' @param as The type of object to return. One of:
#'
#'   * `"InMemoryAnnData"`: Read the H5AD file into memory as an
#'     [`InMemoryAnnData`] object
#'   * `"HDF5AnnData"`: Read the H5AD file as an [`HDF5AnnData`] object
#'   * `"SingleCellExperiment"`: Read the H5AD file as a
#'     [`SingleCellExperiment::SingleCellExperiment`] object
#'   * `"Seurat"`: Read the H5AD file as a
#'     [`SeuratObject::Seurat`] object
#' @param mode The mode to open the HDF5 file.
#'
#'   * `a` creates a new file or opens an existing one for read/write.
#'   * `r` opens an existing file for reading.
#'   * `r+` opens an existing file for read/write.
#'   * `w` creates a file, truncating any existing ones.
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists.
#' @param backed Whether to read the H5AD file in backed mode, returning an
#'  object containing [DelayedArray::DelayedMatrix] matrices. Which slots are
#'  backed depends on the value of `as`.
#' @param ... Extra arguments provided to the `as_*` conversion function for the
#'   object specified by `as`
#'
#' @return The object specified by `as`
#' @export
#'
#' @family AnnData creators
#'
#' @examples
#' h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
#'
#' # Read the H5AD as a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment")
#' }
#'
#' # Read the H5AD as a Seurat object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seurat <- read_h5ad(h5ad_file, as = "Seurat")
#' }
read_h5ad <- function(
  path,
  as = c("InMemoryAnnData", "HDF5AnnData", "SingleCellExperiment", "Seurat"),
  mode = c("r", "r+", "a", "w", "w-", "x"),
  backed = FALSE,
  ...
) {
  as <- match.arg(as)
  mode <- match.arg(mode)

  hdf5_adata <- HDF5AnnData$new(path, mode = mode, backed = backed)

  if (as == "HDF5AnnData") {
    return(hdf5_adata)
  }

  hdf5_adata$open(readonly = TRUE)
  obj <- switch(
    as,
    "SingleCellExperiment" = hdf5_adata$as_SingleCellExperiment(...),
    "Seurat" = hdf5_adata$as_Seurat(...),
    "InMemoryAnnData" = hdf5_adata$as_InMemoryAnnData(...)
  )
  hdf5_adata$close()

  obj
}
