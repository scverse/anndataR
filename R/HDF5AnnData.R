#' @title HDF5AnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
HDF5AnnData <- R6::R6Class("HDF5AnnData",
  inherit = AbstractAnnData,
  private = list(
    .h5obj = NULL
  ),
  active = list(
    #' @description The X slot
    X = function(value) {
      if (missing(value)) {
        # return X
      } else {
        # set X
      }
    },
    #' @description The obs slot
    obs = function(value) {
      if (missing(value)) {
        # return obs
      } else {
        # set obs
      }
    },
    #' @description The var slot
    var = function(value) {
      if (missing(value)) {
        # return var
      } else {
        # set var
      }
    }
  ),
  public = list(
    #' @description HDF5AnnData constructor
    #' 
    #' @param h5obj The rhdf5 object
    initialize = function(h5obj) {
      private$.h5obj <- h5obj
    }
  )
)
