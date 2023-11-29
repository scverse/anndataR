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
#' adata <- AnnData(
#'   X = matrix(1:15, 3L, 5L),
#'   layers = list(
#'     A = matrix(15:1, 3L, 5L),
#'     B = matrix(letters[1:15], 3L, 5L)
#'   ),
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5),
#'   obs_names = LETTERS[1:3],
#'   var_names = letters[1:5]
#' )
#' h5ad_file <- tempfile(fileext = ".h5ad")
#' write_h5ad(adata, h5ad_file)
#'
#' # Write a SingleCellExperiment as an H5AD
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   ncells <- 100
#'   counts <- matrix(rpois(20000, 5), ncol = ncells)
#'   logcounts <- log2(counts + 1)
#'
#'   pca <- matrix(runif(ncells * 5), ncells)
#'   tsne <- matrix(rnorm(ncells * 2), ncells)
#'
#'   sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts, logcounts = logcounts),
#'     reducedDims = list(PCA = pca, tSNE = tsne)
#'   )
#'   
#'   h5ad_file <- tempfile(fileext = ".h5ad")
#'   write_h5ad(sce, h5ad_file)
#' }
#'
#' # Write a Seurat as a H5AD
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   counts <- matrix(1:15, 3L, 5L)
#'   dimnames(counts) <- list(
#'     letters[1:3],
#'     LETTERS[1:5]
#'   )
#'   gene.metadata <- data.frame(
#'     row.names = LETTERS[1:5],
#'     gene = 1:5
#'   )
#'   obj <- SeuratObject::CreateSeuratObject(counts, meta.data = gene.metadata)
#'   cell.metadata <- data.frame(
#'     row.names = letters[1:3],
#'     cell = 1:3
#'   )
#'   obj <- SeuratObject::AddMetaData(obj, cell.metadata)
#'
#'   h5ad_file <- tempfile(fileext = ".h5ad")
#'   write_h5ad(obj, h5ad_file)
#' }
write_h5ad <- function(object, path) {
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
