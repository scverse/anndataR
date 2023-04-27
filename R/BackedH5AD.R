#' @title BackedH5AD
#'
#' @description
#' Representation of an anndata object, disk-backed.
#' @importFrom R6 R6Class
BackedH5AD <- R6Class("BackedH5AD",
    public = list(
      inherit = AbstractAnnData,
      h5object = NULL,

      initialize = function(h5object = NULL){
        self$h5object = h5object
      },
      
      X = function(){
        
      },
      
      obs = function(){
        
      },
      
      var = function(){
        
      }
      
    )
  )