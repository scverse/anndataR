#' Read Zarr
#'
#' Read data from a Zarr store
#'
#' @param path Path to the Zarr store to read
#' @param as The type of object to return. One of:
#'
#'   * `"InMemoryAnnData"`: Read the Zarr store into memory as an
#'     [`InMemoryAnnData`] object
#'   * `"ZarrAnnData"`: Read the Zarr store as an [`ZarrAnnData`] object
#'   * `"SingleCellExperiment"`: Read the Zarr store as a
#'     [`SingleCellExperiment::SingleCellExperiment`] object
#'   * `"Seurat"`: Read the Zarr store as a
#'     [`SeuratObject::Seurat`] object
#' @param mode The mode to open the Zarr file.
#'
#'   * `a` creates a new file or opens an existing one for read/write.
#'   * `r` opens an existing file for reading.
#'   * `r+` opens an existing file for read/write.
#'   * `w` creates a file, truncating any existing ones.
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists.
#' @param backed Whether to read the Zarr store in backed mode, returning an
#'  object containing [DelayedArray::DelayedMatrix] matrices. Which slots are
#'  backed depends on the value of `as`.
#' @param ... Extra arguments provided to the `as_*` conversion function for the
#'   object specified by `as`
#'
#' @details
#'
#' ## Large sparse matrices
#'
#' A `dgCMatrix` or `dgRMatrix` can hold at most 2^31 - 1 non-zero values. Sparse
#' matrices with more non-zero values are read as an `SVT_SparseMatrix`, which
#' requires the **SparseArray** package. To avoid reading these matrices into
#' memory, use `backed = TRUE` with any `as` other than `"InMemoryAnnData"`, which
#' always holds its matrices in memory.
#' See the
#' [large matrices article](https://anndataR.scverse.org/articles/large_matrices.html)
#' for what works where, and how much
#' memory it takes.
#'
#' @return The object specified by `as`
#' @export
#'
#' @family AnnData creators
#'
#' @examples
#' # Please use "example_v3.zarr.zip" for AnnData stored as Zarr version 3
#' zarr_dir <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
#' td <- tempdir(check = TRUE)
#' unzip(zarr_dir, exdir = td)
#' zarr_store <- file.path(td, "example_v2.zarr")
#'
#' # Read the Zarr as a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   sce <- read_zarr(zarr_store, as = "SingleCellExperiment")
#' }
#'
#' # Read the Zarr as a Seurat object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   seurat <- read_zarr(zarr_store, as = "Seurat")
#' }
read_zarr <- function(
  path,
  as = c("InMemoryAnnData", "ZarrAnnData", "SingleCellExperiment", "Seurat"),
  mode = c("r", "r+", "a", "w", "w-", "x"),
  backed = FALSE,
  ...
) {
  as <- match.arg(as)
  mode <- match.arg(mode)

  zarr_adata <- ZarrAnnData$new(path, mode = mode, backed = backed)

  if (as == "ZarrAnnData") {
    return(zarr_adata)
  }

  adata <- switch(
    as,
    "SingleCellExperiment" = zarr_adata$as_SingleCellExperiment(...),
    "Seurat" = zarr_adata$as_Seurat(...),
    "InMemoryAnnData" = zarr_adata$as_InMemoryAnnData(...)
  )

  adata
}
