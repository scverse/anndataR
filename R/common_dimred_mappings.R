#' Get dimensionality reduction mappings keyed by a specific framework
#'
#' @param from Character string specifying which framework to use as keys.
#'   One of: "sce", "seurat", "anndata_obsm", "anndata_varm", "anndata_uns"
#' @return A named list where names are the values from the specified framework
#'   and values are the complete mapping information for each dimred
#' @noRd
.get_dimred_mapping <- function(
  from = c("sce", "seurat", "anndata_obsm", "anndata_varm", "anndata_uns")
) {
  from <- match.arg(from)

  all_mappings <- .get_common_dimred_mappings()
  result <- list()

  for (mapping in all_mappings) {
    key <- mapping[[from]]
    if (!is.null(key)) {
      result[[key]] <- mapping
    }
  }

  result
}

#' Common dimensionality reduction mappings between AnnData, SingleCellExperiment, and Seurat
#'
#' This table defines the standard mappings between scanpy/AnnData naming
#' conventions, Bioconductor/SingleCellExperiment naming conventions, and
#' Seurat naming conventions for common dimensionality reduction techniques.
#'
#' Based on the conventions documented in dr.md:
#' - AnnData uses "X_pca", "X_tsne", "X_umap" in obsm and "PCs" in varm
#' - Seurat uses lowercase "pca", "tsne", "umap"
#' - SingleCellExperiment uses uppercase "PCA", "TSNE", "UMAP"
#'
#' @return A list of dimensionality reduction mappings, where each element contains
#'   the naming conventions across frameworks and associated metadata
#' @noRd
.get_common_dimred_mappings <- function() {
  list(
    # Core dimensionality reductions with documented naming conventions
    list(
      sce = "PCA",
      seurat = "pca",
      anndata_obsm = "X_pca",
      anndata_varm = "PCs",
      anndata_uns = NULL
    ),
    list(
      sce = "tSNE",
      seurat = "tsne",
      anndata_obsm = "X_tsne",
      anndata_varm = NULL,
      anndata_uns = NULL
    ),
    list(
      sce = "UMAP",
      seurat = "umap",
      anndata_obsm = "X_umap",
      anndata_varm = NULL,
      anndata_uns = NULL
    )
  )
}
