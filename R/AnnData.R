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
  obs_names = NULL,
  var_names = NULL,
  X = NULL,
  obs = NULL,
  var = NULL,
  layers = NULL
) {
  InMemoryAnnData$new(
    obs_names = obs_names,
    var_names = var_names,
    X = X,
    obs = obs,
    var = var,
    layers = layers
  )
}