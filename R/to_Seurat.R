#' Convert an AnnData object to Seurat
#' 
#' @param obj An AnnData object
#' 
#' @importFrom methods is
#' @importFrom Seurat CreateAssayObject AddMetaData CreateSeuratObject
#' 
#' @export
#' @examples

#' @examples
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(cell = 1:3),
#'   obs_names = letters[1:3],
#'   var = data.frame(gene = 1:5),
#'   var_names = letters[1:5]
#' )
#' ad$to_Seurat()
to_Seurat <- function(obj) {
  if (!is(obj, "AbstractAnnData")) stop("Argument 'obj' should be an AnnData object.")

  # create assay object
  X_assay_metadata <- obj$var
  rownames(X_assay_metadata) <- obj$var_names
  X_assay <- Seurat::CreateAssayObject(counts = obj$X)
  X_assay <- Seurat::AddMetaData(X_assay, metadata = X_assay_metadata)

  seurat_obj_metadata <- obj$obs
  rownames(seurat_obj_metadata) <- obj$obs_names
  seurat_obj <- Seurat::CreateSeuratObject(X_assay, meta.data = seurat_obj_metadata)
  
  seurat_obj
}
