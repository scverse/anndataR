#' Create an in-memory AnnData object.
#'
#' For more information on the functionality of an AnnData object, see [AnnData-usage]
#'
#' @param X See the `X` slot in [AnnData-usage]
#' @param layers See the `layers` slot in [AnnData-usage]
#' @param obs See the `obs` slot in [AnnData-usage]
#' @param var See the `var` slot in [AnnData-usage]
#' @param obsm See the `obsm` slot in [AnnData-usage]
#' @param varm See the `varm` slot in [AnnData-usage]
#' @param obsp See the `obsp` slot in [AnnData-usage]
#' @param varp See the `varp` slot in [AnnData-usage]
#' @param uns See the `uns` slot in [AnnData-usage]
#' @param shape Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if both
#'   `X` or `obs` and `var` are not provided.
#'
#' @return An [InMemoryAnnData] object
#' @export
#'
#' @seealso [AnnData-usage] for details of `AnnData` structure and usage
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
