#' List the available AnnData generators.
#' @returns A named list of generators currently available via 
#' \link[anndataR]{get_generator}.
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
#' @returns A requested AnnData generator.
get_generator <- function(class = c("HDF5AnnData", "InMemoryAnnData")) {
  # TODO: also support directly passing the correct class?
  class <- match.arg(class)
  list_generators()[[class]]
}
