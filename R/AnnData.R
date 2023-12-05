#' An in-memory AnnData object
#'
#' @description
#' This class is used to represent an AnnData object in memory.
#' AnnData stores a data matrix `X` together with annotations of
#' observations `obs` (`obsm`, `obsp`), variables `var` (`varm`, `varp`), and
#' unstructured annotations `uns`.
#'
#' To read an AnnData file from disk, use [read_h5ad()] instead.
#'
#' @param obj Object to create an AnnData from. If `obj = NULL` the remaining
#'   arguments are used instead.
#' @param obs_names A vector of unique identifiers
#'   used to identify each row of `obs` and to act as an index into the
#'   observation dimension of the AnnData object. The length of `obs_names`
#'   defines the observation dimension of the AnnData object.
#' @param var_names A vector of unique identifiers used to identify each row
#'   of `var` and to act as an index into the variable dimension of the
#'   AnnData object. The length of `var_names` defines the variable
#'   dimension of the AnnData object.
#' @param X Either `NULL` or a observation × variable matrix with
#'   dimensions consistent with `obs` and `var`.
#' @param layers Either `NULL` or a named list, where each element is an
#'   observation × variable matrix with dimensions consistent with `obs` and
#'   `var`.
#' @param obs Either `NULL` or a `data.frame` with columns containing
#'   information about observations. If `NULL`, an `n_obs`×0 data frame will
#'   automatically be generated.
#' @param var Either `NULL` or a `data.frame` with columns containing
#'   information about variables. If `NULL`, an `n_vars`×0 data frame will
#'   automatically be generated.
#' @param obsm The obsm slot is used to store multi-dimensional annotation
#'   arrays. It must be either `NULL` or a named list, where each element is a
#'   matrix with `n_obs` rows and an arbitrary number of columns.
#' @param varm The varm slot is used to store multi-dimensional annotation
#'   arrays. It must be either `NULL` or a named list, where each element is a
#'   matrix with `n_vars` rows and an arbitrary number of columns.
#' @param obsp The obsp slot is used to store sparse multi-dimensional
#'   annotation arrays. It must be either `NULL` or a named list, where each
#'   element is a sparse matrix where each dimension has length `n_obs`.
#' @param varp The varp slot is used to store sparse multi-dimensional
#'   annotation arrays. It must be either `NULL` or a named list, where each
#'   element is a sparse matrix where each dimension has length `n_vars`.
#' @param uns The uns slot is used to store unstructured annotation. It must
#'   be either `NULL` or a named list.
#' @param ... Additional arguments passed to conversion functions. See
#'   [SingleCellExperiment-Conversion] and [Seurat-Conversion].
#'
#' @export
#'
#' @examples
#' adata <- AnnData(
#'   obs_names = paste0("obs", 1:3),
#'   var_names = paste0("var", 1:4),
#'   X = matrix(1:12, nrow = 3, ncol = 4),
#'   obs = data.frame(
#'     n_counts = c(1, 2, 3),
#'     n_cells = c(1, 2, 3)
#'   ),
#'   var = data.frame(
#'     n_cells = c(1, 2, 3, 4)
#'   )
#' )
#'
#' adata
AnnData <- function(
    obj = NULL,
    obs_names = NULL,
    var_names = NULL,
    X = NULL,
    obs = NULL,
    var = NULL,
    layers = NULL,
    obsm = NULL,
    varm = NULL,
    obsp = NULL,
    varp = NULL,
    uns = NULL,
    ...) {

  if (is.null(obj)) {
    InMemoryAnnData$new(
      obs_names = obs_names,
      var_names = var_names,
      X = X,
      obs = obs,
      var = var,
      layers = layers,
      obsm = obsm,
      varm = varm,
      obsp = obsp,
      varp = varp,
      uns = uns
    )
  } else if (inherits(obj, "SingleCellExperiment")) {
    from_SingleCellExperiment(obj, output_class = "InMemoryAnnData", ...)
  } else if (inherits(obj, "Seurat")) {
    from_Seurat(obj, output_class = "InMemoryAnnData", ...)
  } else {
    stop("No implementation for converting object of class ", class(obj))
  }
}
