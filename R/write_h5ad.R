#' Write H5AD
#'
#' Write an H5AD file
#'
#' @param object The object to write, either a "SingleCellExperiment" or a
#' "Seurat" object
#' @param path Path of the file to write to
#'
#' @return `path` invisibly
#' @export
#'
#' @examples
#' # Write a InMemoryAnnData as a H5AD
#' adata <- dummy_data("InMemoryAnnData")
#' write_h5ad(adata)
#'
#' # Write a SingleCellExperiment as an H5AD
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- dummy_data("SingleCellExperiment")
#'   write_h5ad(sce)
#' }
#'
#' # Write a Seurat as a H5AD
#' if (requireNamespace("SeuratObject", quietly = TRUE)) { 
#'   seur <- dummy_data("Seurat")
#'   write_h5ad(seur)
#' }
write_h5ad <- function(object, 
                       path=tempfile(fileext = ".h5ad")) {
  if (inherits(object, "SingleCellExperiment")) {
    from_SingleCellExperiment(
      object,
      output_class = "HDF5AnnData",
      file = path
    )
  } else if (inherits(object, "Seurat")) {
    from_Seurat(
      object,
      output_class = "HDF5AnnData",
      file = path
    )
  } else if (inherits(object, "AbstractAnnData")) {
    to_HDF5AnnData(object, path)
  } else {
    stop("Unable to write object of class: ", class(object))
  }

  invisible(path)
}
