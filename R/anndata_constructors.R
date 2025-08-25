#' List the available AnnData generators.
#'
#' @noRd
anndata_constructors <- function() {
  list(
    "HDF5AnnData" = HDF5AnnData,
    "InMemoryAnnData" = InMemoryAnnData
  )
}

#' Retrieve the AnnData constructor for a given class.
#'
#' @param class Name of the AnnData class. Must be one of `"HDF5AnnData"`
#' or `"InMemoryAnnData"`.
#'
#' @noRd
get_anndata_constructor <- function(
  class = c("HDF5AnnData", "InMemoryAnnData")
) {
  # TODO: also support directly passing the correct class?
  class <- match.arg(class)
  anndata_constructors()[[class]]
}
