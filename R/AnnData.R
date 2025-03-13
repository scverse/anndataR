#' Create an in-memory AnnData object.
#'
#' For more information on the functionality of an AnnData object, see [AnnData-usage].
#'
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
#' @param shape Shape tuple (#observations, #variables). Can be provided
#'   if `X` or `obs` and `var` are not provided.
#'
#' @return An in-memory AnnData object.
#'
#' @export
#'
#' @examples
#' adata <- AnnData(
#'   X = matrix(1:12, nrow = 3, ncol = 4),
#'   obs = data.frame(
#'     row.names = paste0("obs", 1:3),
#'     n_counts = c(1, 2, 3),
#'     n_cells = c(1, 2, 3)
#'   ),
#'   var = data.frame(
#'     row.names = paste0("var", 1:4),
#'     n_cells = c(1, 2, 3, 4)
#'   )
#' )
#'
#' adata
#'
#' @seealso [AnnData-usage]
AnnData <- function(
  X = NULL,
  obs = NULL,
  var = NULL,
  layers = NULL,
  obsm = NULL,
  varm = NULL,
  obsp = NULL,
  varp = NULL,
  uns = NULL,
  shape = NULL
) {
  InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    layers = layers,
    obsm = obsm,
    varm = varm,
    obsp = obsp,
    varp = varp,
    uns = uns,
    shape = shape
  )
}
