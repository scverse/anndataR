#' Write Zarr
#'
#' Write a Zarr store
#'
#' @param object The object to write, either a "SingleCellExperiment" or a
#' "Seurat" object
#' @param path Path of the file to write to
#' @param compression The compression algorithm to use when writing
#'  Zarr arrays. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
#' `"none"`.
#' @param mode The mode to open the Zarr store.
#'
#'   * `a` creates a new file or opens an existing one for read/write.
#'   * `r+` opens an existing file for read/write.
#'   * `w` creates a file, truncating any existing ones
#'   * `w-`/`x` are synonyms creating a file and failing if it already exists.
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
#' store <- pizzarr::MemoryStore$new()
#' write_zarr(adata, store)
#'
#' # Write a SingleCellExperiment as a Zarr store
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
#'   store <- pizzarr::MemoryStore$new()
#'   write_zarr(sce, store)
#' }
#'
#' # Write a Seurat as a Zarr store
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   # TODO: uncomment this code when the seurat converter is fixed
#'   # counts <- matrix(1:15, 3L, 5L)
#'   # dimnames(counts) <- list(
#'   #   letters[1:3],
#'   #   LETTERS[1:5]
#'   # )
#'   # gene.metadata <- data.frame(
#'   #   row.names = LETTERS[1:5],
#'   #   gene = 1:5
#'   # )
#'   # obj <- SeuratObject::CreateSeuratObject(counts, meta.data = gene.metadata)
#'   # cell.metadata <- data.frame(
#'   #   row.names = letters[1:3],
#'   #   cell = 1:3
#'   # )
#'   # obj <- SeuratObject::AddMetaData(obj, cell.metadata)
#'   #
#'   # store <- pizzarr::MemoryStore$new()
#'   # write_zarr(obj, store)
#' }
write_zarr <- function(object, store, compression = c("none", "gzip", "lzf"),
                       mode = c("w-", "r", "r+", "a", "w", "x")) {
  mode <- match.arg(mode)
  if (inherits(object, "SingleCellExperiment")) {
    from_SingleCellExperiment(
      object,
      output_class = "ZarrAnnData",
      store = store,
      compression = compression,
      mode = mode
    )
  } else if (inherits(object, "Seurat")) {
    from_Seurat(
      object,
      output_class = "ZarrAnnData",
      store = store,
      compression = compression,
      mode = mode
    )
  } else if (inherits(object, "AbstractAnnData")) {
    to_ZarrAnnData(object, store, compression = compression, mode = mode)
  } else {
    stop("Unable to write object of class: ", class(object))
  }
}
