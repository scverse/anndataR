#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom lifecycle deprecated
#' @importFrom methods as new
#' @importFrom purrr map_lgl map_dfr
#' @importFrom R6 R6Class
#' @importFrom rlang `%||%`
#' @importFrom stats setNames
#' @importFrom utils tail
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  # Check if the R anndata package is loaded and warn about conflicts
  if ("anndata" %in% loadedNamespaces()) {
    cli::cli_warn(c(
      "The R {.pkg anndata} package is loaded and may interfere with {.pkg anndataR}'s Python integration.",
      "i" = "This may cause automatic conversion of Python AnnData objects to R6 classes.",
      "i" = "Consider restarting R and loading {.pkg anndataR} before {.pkg anndata} to avoid conflicts.",
      "i" = "Use {.code reticulate::import('anndata', convert = FALSE)} to prevent automatic conversion."
    ))
  }
}
