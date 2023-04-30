#' Convert an AnnData object to Seurat
#'
#' @param obj An AnnData object
#'
#' @importFrom Matrix t
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
#' to_Seurat(ad)
# TODO: fix snake_case + CamelCase
to_Seurat <- function(obj) { # nolint
  requireNamespace("SeuratObject")

  stopifnot(inherits(obj, "AbstractAnnData"))

  # translate var_names
  var_names_ <- .toseurat_check_obsvar_names(obj$var_names, "var_names")
  obs_names_ <- .toseurat_check_obsvar_names(obj$obs_names, "obs_names")

  # translate var
  var_ <- obj$var
  rownames(var_) <- var_names_

  # translate obs
  obs_ <- obj$obs
  rownames(obs_) <- obs_names_

  # translate X
  x_ <- Matrix::t(obj$X)
  dimnames(x_) <- list(var_names_, obs_names_)

  # create assay object
  x_assay <- SeuratObject::CreateAssayObject(counts = x_)
  x_assay <- SeuratObject::AddMetaData(x_assay, metadata = var_)

  # create seurat object
  seurat_obj <- SeuratObject::CreateSeuratObject(x_assay, meta.data = obs_)

  seurat_obj
}

.toseurat_check_obsvar_names <- function(names, label) {
  if (any(grepl("_", names))) {
    # mimic seurat behaviour
    warning(wrap_message(
      "'", label, "' ",
      "cannot have underscores ('_') when converting to Seurat, ",
      "replacing with dashes ('-')"
    ))
    names <- gsub("_", "-", names)
  }

  names
}
