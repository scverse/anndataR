#' Convert an AnnData object to Seurat
#' 
#' @param obj An AnnData object
#' 
#' @importFrom methods is
#' @importFrom Seurat CreateAssayObject AddMetaData CreateSeuratObject
to_seurat <- function(obj) {
  if (!is(obj, "AbstractAnnData")) stop("Argument 'obj' should be an AnnData object.")

  # create assay object
  X_assay <- Seurat::CreateAssayObject(counts = obj$X)
  X_assay <- Seurat::AddMetaData(X_assay, metadata = obj$var)

  seurat_obj <- Seurat::CreateSeuratObject(X_assay, meta.data = obj$obs)
  
  seurat_obj
}
