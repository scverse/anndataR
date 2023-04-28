#' Convert an AnnData object to Seurat
#' 
#' @param obj An AnnData object
#' 
#' @export
#' @examples
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(cell = 1:3),
#'   obs_names = letters[1:3],
#'   var = data.frame(gene = 1:5),
#'   var_names = letters[1:5]
#' )
#' to_Seurat(obj)
to_Seurat <- function(obj) {
  requireNamespace("Seurat")

  stopifnot(inherits(obj, "AbstractAnnData"))

  # translate var
  var_ <- obj$var
  rownames(var_) <- obj$var_names

  # translate obs
  obs_ <- obj$obs
  rownames(obs_) <- obj$obs_names

  # create assay object
  X_assay <- Seurat::CreateAssayObject(counts = obj$X)
  X_assay <- Seurat::AddMetaData(X_assay, metadata = var_)

  # create seurat object
  seurat_obj <- Seurat::CreateSeuratObject(X_assay, meta.data = obs_)
  
  seurat_obj
}
