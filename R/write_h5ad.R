#' Write H5AD
#'
#' Write an H5AD file
#'
#' @param object The object to write, either a "SingleCellExperiment" or a
#' "Seurat" object
#' @param path Path of the file to write to
#' @param compression The compression algorithm to use when writing the
#'  HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
#' `"none"`.
#' @param ... Additional arguments passed to conversion functions. See
#'   [SingleCellExperiment-Conversion] and [Seurat-Conversion].
#'
#' @return `path` invisibly
#' @export
#'
#' @examples
#' # Write a SingleCellExperiment as an H5AD
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- generate_dataset(example = TRUE, format = "SingleCellExperiment")
#'   h5ad_file <- tempfile(fileext = ".h5ad")
#'   write_h5ad(sce, h5ad_file)
#' }
#'
#' # Write a Seurat as a H5AD
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   # TODO: uncomment this code when the seurat converter is fixed
#'   # seurat <- generate_dataset(example = TRUE, format = "Seurat")
#'   # h5ad_file <- tempfile(fileext = ".h5ad")
#'   # write_h5ad(seurat, h5ad_file)
#' }
write_h5ad <- function(object, path, compression = c("none", "gzip", "lzf"), ...) {
  if (inherits(object, "SingleCellExperiment")) {
    from_SingleCellExperiment(
      object,
      output_class = "HDF5AnnData",
      file = path,
      compression = compression,
      ...
    )
  } else if (inherits(object, "Seurat")) {
    from_Seurat(
      object,
      output_class = "HDF5AnnData",
      file = path,
      compression = compression,
      ...
    )
  } else if (inherits(object, "AbstractAnnData")) {
    to_HDF5AnnData(object, path, compression = compression)
  } else {
    stop("Unable to write object of class: ", class(object))
  }

  invisible(path)
}
