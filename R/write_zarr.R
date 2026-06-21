#' Write Zarr
#'
#' Write a Zarr file
#'
#' @param object The object to write, either a
#'   [`SingleCellExperiment::SingleCellExperiment`] or a
#'   [`SeuratObject::Seurat`] object
#' @param path Path of the file to write to
#' @param compression The compression algorithm to use when writing the Zarr
#'   file. Can be one of `"none"`, `"gzip"`, `"blosc"`, `"zstd"`,
#'   `"lzma"`, `"bz2"`, `"zlib"`, `"lz4"`. Defaults to `"none"`.
#'   See `help("compressors", package = "Rarr")`.
#' @param mode The mode to open the Zarr file.
#'
#'   * `a` creates a new file or opens an existing one for read/write
#'   * `r+` opens an existing file for read/write
#'   * `w` creates a file, truncating any existing ones
#'   * `w-`/`x` are synonyms creating a file and failing if it already exists
#' @param ... Additional arguments passed to [as_AnnData()]
#'
#' @details
#'
#' ## `NULL` values
#'
#' For compatibility with changes in Python **anndata** 0.12.0, `NULL` values
#' in `uns` are written to Zarr files as a `NULL` dataset (instead of not being
#' written at all). To disable this behaviour, set
#' `option(anndataR.write_null = FALSE)`. This may be required to allow the file
#' to be read by older versions of Python **anndata**.
#'
#' @return `path` invisibly
#' @export
#'
#' @examples
#' adata <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' zarr_store <- tempfile(fileext = ".zarr")
#' adata$write_zarr(zarr_store)
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
#'   adata <- as_AnnData(sce)
#'   zarr_store <- tempfile(fileext = ".zarr")
#'   adata$write_zarr(zarr_store)
#' }
#'
#' # Write a Seurat as a Zarr
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   library(Seurat)
#'
#'   counts <- matrix(1:15, 5L, 3L)
#'   dimnames(counts) <- list(
#'     LETTERS[1:5],
#'     letters[1:3]
#'   )
#'   cell.metadata <- data.frame(
#'     row.names = letters[1:3],
#'     cell = 1:3
#'   )
#'   obj <- CreateSeuratObject(counts, meta.data = cell.metadata)
#'   gene.metadata <- data.frame(
#'     row.names = LETTERS[1:5],
#'     gene = 1:5
#'   )
#'   obj[["RNA"]] <- AddMetaData(GetAssay(obj), gene.metadata)
#'
#'   adata <- as_AnnData(obj)
#'   zarr_store <- tempfile(fileext = ".zarr")
#'   adata$write_zarr(zarr_store)
#' }
write_zarr <- function(
  object,
  path,
  compression = c(
    "none",
    "gzip",
    "blosc",
    "zstd",
    "lzma",
    "bz2",
    "zlib",
    "lz4"
  ),
  mode = c("w-", "r", "r+", "a", "w", "x"),
  ...
) {
  mode <- match.arg(mode)
  adata <- if (inherits(object, "AbstractAnnData")) {
    object$as_ZarrAnnData(
      path,
      compression = compression,
      mode = mode
    )
  } else {
    as_AnnData(
      object,
      output_class = "ZarrAnnData",
      file = path,
      compression = compression,
      mode = mode,
      ...
    )
  }

  rm(adata)
  gc()

  invisible(path)
}
