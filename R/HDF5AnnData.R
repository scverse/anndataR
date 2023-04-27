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
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        # return X
      } else {
        # set X
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        # return obs
      } else {
        # set obs
      }
    },
    #' @field var The var slot
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
      attrs <- rhdf5::h5readAttributes(h5obj, "/")
      
      if (!("encoding-type") %in% names(attrs) ||
          !("encoding-version" %in% names(attrs))) {
        stop(
          "H5AD files without encodings are not supported ",
          "(this file may have been created with Python anndata prior to v0.8.0)"
        )
      }
      
      private$.h5obj <- h5obj
    }
  )
)
