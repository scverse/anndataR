#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#' @importFrom R6 R6Class
InMemoryAnnData <- R6Class("InMemoryAnnData", 
  public = list(
    inherit = AbstractAnnData,
    X = NULL,
    obs = NULL,
    var = NULL,
  
    #' @description Creates a new instance of an in memory AnnData object. Inherits from AbstractAnnData
    #' @param X The count matrix
    #' @param obs The associated observation metadata
    #' @param var The associated variables metadata
    initialize = function(X = NULL, obs = NULL, var = NULL) {
     self$X = X
     self$obs = obs
     self$var = var
    }
  )
)