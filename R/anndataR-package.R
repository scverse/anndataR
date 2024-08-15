#' @keywords internal
#' @importFrom methods as
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
#' @section An AnnData object has the following slots:
#'
#' Access them by using the `$` operator. Example: `adata$obs`.
#'
#' * `X`: A matrix of observations by variables.
#' * `obs`: A data frame of observations.
#' * `var`: A data frame of variables.
#' * `obs_names`: Names of observations (alias for `rownames(obs)`).
#' * `var_names`: Names of variables (alias for `rownames(var)`).
#' * `layers`: A named list of matrices with the same dimensions as `X`.
#' * `obsm`: A named list of matrices with the same number of rows as `obs`.
#' * `varm`: A named list of matrices with the same number of rows as `var`.
#' * `obsp`: A named list of sparse matrices with the same number of rows and columns as the number of observations.
#' * `varp`: A named list of sparse matrices with the same number of rows and columns as the number of variables.
#' * `uns`: A named list of unstructured annotations.
#'
#' @section An AnnData object has the following methods:
#'
#' Access them by using the `$` operator. Example: `adata$write_h5ad()`.
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
#' * `uns_keys()`: Element names of `uns`.
#' * `write_h5ad()`: Write the AnnData object to an HDF5 file.
#'
#' @section Conversion methods:
#'
#' Access them by using the `$` operator. Example: `adata$to_Seurat()`.
#'
#' * `to_SingleCellExperiment()`: Convert to SingleCellExperiment.
#' * `to_Seurat()`: Convert to Seurat.
#' * `to_InMemoryAnnData()`: Convert to an InMemory AnnData.
#' * `to_HDF5AnnData()`: Convert to an HDF5 Backed AnnData.
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
