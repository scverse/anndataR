#' AnnData structure and usage
#'
#' @description
#' The `AnnData` object stores a data matrix `X` together with annotations of
#' observations `obs` (`obsm`, `obsp`) and variables `var` (`varm`, `varp`).
#' Additional layers of data can be stored in `layers` and unstructured
#' annotations in `uns`.
#'
#' ## Back ends
#'
#' There are different back ends for `AnnData` objects that inherit from the
#' abstract [AbstractAnnData] class. For example, the [InMemoryAnnData] stores
#' data in memory or the [HDF5AnnData] is backed by a H5AD file.
#'
#' ## Usage
#'
#' The items listed as **"Slots"** are elements of the `AnnData` object that
#' contain data and can be accessed or set. **"Fields"** return information
#' about the `AnnData` object but cannot be set directly. Both, as well as
#' methods, can be accessed using the `$` operator
#'
#' For example:
#'
#' - `adata$X` returns the `X` matrix
#' - `adata$X <- x` sets the `X` matrix
#' - `adata$method()` calls a method
#'
#' @slot X The main data matrix. Either `NULL` or an observation x variable
#'   matrix (without dimnames) with dimensions consistent with `n_obs` and
#'   `n_vars`.
#' @slot layers Additional data layers. Must be `NULL` or a named list of
#'   matrices having dimensions consistent with `n_obs` and `n_vars`.
#' @slot obs Observation annotations. A `data.frame` with columns containing
#'   information about observations. The number of rows of `obs` defines the
#'   observation dimension of the `AnnData` object (`n_obs`). If `NULL`, an
#'   `n_obs` × 0 `data.frame` will automatically be generated.
#' @slot var Variable annotations. A `data.frame` with columns containing
#'   information about variables. The number of rows of `var` defines the
#'   variable dimension of the `AnnData` object (`n_vars`). If `NULL`, an
#'   `n_vars` × 0 `data.frame` will automatically be generated.
#' @slot obs_names Observation names. Either `NULL` or a vector of unique
#'   identifiers used to identify each row of `obs` and to act as an index into
#'   the observation dimension of the `AnnData` object. For compatibility with
#'   _R_ representations, `obs_names` should be a unique character vector.
#' @slot var_names Variable names. Either `NULL` or a vector of unique
#'   identifiers used to identify each row of `var` and to act as an index into
#'   the variable dimension of the `AnnData` object. For compatibility with _R_
#'   representations, `var_names` should be a unique character vector.
#' @slot obsm Multi-dimensional observation annotation. Must be `NULL` or a
#'   named list of array-like elements with number of rows equal to `n_obs`.
#' @slot varm Multi-dimensional variable annotations. Must be `NULL` or a named
#'   list of array-like elements with number of rows equal to `n_vars`.
#' @slot obsp Observation pairs. Must be `NULL` or a named list of array-like
#'   elements with number of rows and columns equal to `n_obs`.
#' @slot varp Variable pairs. Must be `NULL` or a named list of array-like
#'   elements with number of rows and columns equal to `n_vars`.
#' @slot uns Unstructured annotations. Must be `NULL` or a named list.
#'
#' @field shape Dimensions (observations x variables) of the `AnnData` object
#' @field n_obs Number of observations
#' @field n_vars Number of variables
#' @field obs_keys Keys (column names) of `obs`
#' @field var_keys Keys (column names) of `var`
#' @field layers_keys Keys (element names) of `layers`
#' @field obsm_keys Keys (element names) of `obsm`
#' @field varm_keys Keys (element names) of `varm`
#' @field obsp_keys Keys (element names) of `obsp`
#' @field varp_keys Keys (element names) of `varp`
#' @field uns_keys Keys (element names) of `uns`
#'
#' @section Methods:
#'
#' ## Conversion methods:
#'
#' \describe{
#'   \item{
#'     `as_SingleCellExperiment()`
#'   }{
#'     Convert to `SingleCellExperiment`, see [as_SingleCellExperiment()]
#'   }
#'   \item{`to_Seurat()`}{Convert to `Seurat`, see [to_Seurat()]}
#'   \item{`to_InMemoryAnnData()`}{Convert to [InMemoryAnnData]}
#'   \item{`to_HDF5AnnData()`}{Convert to [HDF5AnnData], see [to_HDF5AnnData()]}
#' }
#'
#' ## Output methods:
#'
#' \describe{
#'   \item{
#'     `write_h5ad()`
#'   }{
#'     Write the `AnnData` object to an HDF5 file, see [write_h5ad()]
#'   }
#' }
#'
#' ## General methods:
#'
#' \describe{
#'   \item{`print()`}{Print a summary of the `AnnData` object}
#' }
#'
#' @section Functions that can be used to create AnnData objects:
#'
#' \describe{
#'   \item{[AnnData()]}{Create an [InMemoryAnnData] object}
#'   \item{[read_h5ad()]}{Read an `AnnData` from a H5AD file}
#'   \item{[as_AnnData()]}{Convert other objects to an `AnnData` object}
#' }
#'
#' @seealso The documentation for the Python `anndata` package
#'   <https://anndata.readthedocs.io/en/stable/>
#' @seealso [AbstractAnnData] for the abstract class that all `AnnData` objects
#'   inherit from
#' @seealso [InMemoryAnnData] for the in-memory implementation of `AnnData`
#' @seealso [HDF5AnnData] for the HDF5-backed implementation of `AnnData`
#'
#' @name AnnData-usage
NULL
