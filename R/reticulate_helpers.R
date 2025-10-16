#' Reticulate Helper Functions for AnnData Conversion
#'
#' This file contains helper functions that enable seamless conversion between
#' R and Python AnnData objects using the reticulate package. These functions
#' provide automatic S3 method dispatch for converting AnnData objects across
#' the R-Python boundary.
#'
#' The main conversion functions include:
#'
#' * `py_to_r.anndata._core.anndata.AnnData`: Converts Python AnnData
#'   objects to R [ReticulateAnnData] objects
#' * `r_to_py.AbstractAnnData`: Converts R [AbstractAnnData] objects
#'   to Python AnnData objects
#' * `py_to_r.collections.abc.Mapping`: Converts Python mapping
#'   objects to R lists
#'
#' These functions are automatically registered as S3 methods and are called
#' when using `reticulate::py_to_r()` and `reticulate::r_to_py()` on compatible
#' objects.
#'
#' @name reticulate-helpers
#' @family object converters
#'
#' @examples
#' \donttest{
#' # Requires Python anndata to be installed
#' if (requireNamespace("reticulate", quietly = TRUE) &&
#'       reticulate::py_module_available("anndata")) {
#'
#'   library(reticulate)
#'
#'   # Create Python AnnData object
#'   ad_py <- import("anndata", convert = FALSE)
#'   py_adata <- ad_py$AnnData(X = r_to_py(matrix(1:12, 3, 4)))
#'
#'   # Automatic conversion to R (uses py_to_r.anndata._core.anndata.AnnData)
#'   r_adata <- py_to_r(py_adata)
#'
#'   # Automatic conversion back to Python (uses r_to_py.AbstractAnnData)
#'   py_adata2 <- r_to_py(r_adata)
#' }
#' }
NULL


#' @rdname reticulate-helpers
#' @importFrom reticulate py_to_r
#' @method py_to_r collections.abc.Mapping
#' @export
py_to_r.collections.abc.Mapping <- function(x) {
  out <- list()
  bi <- reticulate::import_builtins()
  keys <- bi$list(x$keys())
  for (name in keys) {
    out[[name]] <- reticulate::py_to_r(x[[name]])
  }
  out
}

#' @rdname reticulate-helpers
#' @param x A Python AnnData object
#' @return A [`ReticulateAnnData`] object wrapping the Python object
#' @export
#' @importFrom reticulate py_to_r
#' @method py_to_r anndata._core.anndata.AnnData
py_to_r.anndata._core.anndata.AnnData <- function(x) {
  ReticulateAnnData$new(py_anndata = x)
}

#' @rdname reticulate-helpers
#' @param x An AbstractAnnData object (any anndataR implementation)
#' @param convert Whether to convert the result (passed to reticulate)
#' @return A Python AnnData object
#' @export
#' @importFrom reticulate r_to_py
#' @method r_to_py AbstractAnnData
r_to_py.AbstractAnnData <- function(x, convert = TRUE) {
  if (inherits(x, "ReticulateAnnData")) {
    # If it's already a ReticulateAnnData, return the underlying Python object
    return(x$py_anndata())
  } else {
    # Convert other AnnData types to ReticulateAnnData first, then extract Python object
    ret_adata <- as_ReticulateAnnData(x)
    return(ret_adata$py_anndata())
  }
}
