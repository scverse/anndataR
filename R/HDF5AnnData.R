#' @title HDF5AnnData
#'
#' @description
#' Representation of an anndata object, disk-backed.
HDF5AnnData <- R6Class("HDF5AnnData",
  inherit = AbstractAnnData,
  private = list(
    .h5object = NULL
  ),
  public = list(
    initialize = function(h5object){
      private$.h5object <- h5object
    },
    
    X = function() {
      
    },
    
    obs = function() {
      
    },
    
    var = function() {
      
    }
  )
)