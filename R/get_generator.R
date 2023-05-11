#' List the available AnnData generators.
list_generators <- function() {
  list(
    "HDF5AnnData" = HDF5AnnData,
    "InMemoryAnnData" = InMemoryAnnData
  )
}

#' Fetch an AnnData generator.
#'
#' @param class Name of the AnnData class. Must be one of `"HDF5AnnData"`
#' or `"InMemoryAnnData"`.
get_generator <- function(class = c("HDF5AnnData", "InMemoryAnnData")) {
  # TODO: also support directly passing the correct class?
  class <- match.arg(class)
  list_generators()[[class]]
}
