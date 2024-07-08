#' An AnnData object
#'
#' @description An AnnData object. This class can either be an in-memory
#' AnnData (InMemoryAnnData) or an HDF5-backed AnnData (HDF5AnnData). The
#' AnnData object stores a data matrix `X` together with annotations of
#' observations `obs` (`obsm`, `obsp`) and variables `var` (`varm`, `varp`).
#' Additional layers of data can be stored in `layers` and unstructured
#' annotations in `uns`.
#'
#' @section Functions that can be used to create AnnData objects:
#'
#'   * [AnnData()]: Create an in-memory AnnData object.
#'   * [read_h5ad()]: Read an HDF5-backed AnnData file from disk.
#'   * [from_SingleCellExperiment()]: Convert a SingleCellExperiment object to an AnnData object.
#'   * [from_Seurat()]: Convert a Seurat object to an AnnData object.
#'
#' @section Slots:
#'
#' * `X`: A matrix of observations by variables.
#' * `obs`: A data frame of observations.
#' * `var`: A data frame of variables.
#' * `layers`: A named list of matrices with the same dimensions as `X`.
#' * `obsm`: A named list of matrices with the same number of rows as `obs`.
#' * `varm`: A named list of matrices with the same number of rows as `var`.
#' * `obsp`: A named list of sparse matrices with the same number of rows and columns as the number of observations.
#' * `varp`: A named list of sparse matrices with the same number of rows and columns as the number of variables.
#' * `uns`: A named list of unstructured annotations.
#'
#' @section Methods:
#'
#' * `print()`: Print a summary of the AnnData object.
#' * `shape()`: Dimensions (observations x variables) of the AnnData object.
#' * `n_obs()`: Number of observations in the AnnData object.
#' * `n_vars()`: Number of variables in the AnnData object.
#' * `obs_keys()`: Column names of `obs`.
#' * `var_keys()`: Column names of `var`.
#' * `layers_keys()`: Element names of `layers`.
#' * `obsm_keys()`: Element names of `obsm`.
#' * `varm_keys()`: Element names of `varm`.
#' * `obsp_keys()`: Element names of `obsp`.
#' * `varp_keys()`: Element names of `varp`.
#'
#' @section Conversion methods:
#'
#' * `to_SingleCellExperiment()`: Convert to SingleCellExperiment.
#' * `to_Seurat()`: Convert to Seurat.
#' * `to_InMemoryAnnData()`: Convert to an InMemory AnnData.
#' * `to_HDF5AnnData()`: Convert to an HDF5 Backed AnnData.
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
#' @param obs_names Names of observations (alias for `rownames(obs)`).
#' @param var_names Names of variables (alias for `rownames(var)`).
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
#' @return An [AbstractAnnData] object.
#'
#' @seealso [AbstractAnnData]
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
    shape = shape) {
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
