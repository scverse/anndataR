#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
InMemoryAnnData <- R6Class("InMemoryAnnData",
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .obs = NULL,
    .var = NULL
  ),
  public = list(
    X = function() {
      self$.X
    },
    obs = function() {
      self$.obs
    },
    var = function() {
      self$.var
    },
  
    #' @description Creates a new instance of an in memory AnnData object.
    #' Inherits from AbstractAnnData
    #' 
    #' @param X The count matrix
    #' @param obs The associated observation metadata
    #' @param var The associated variables metadata
    initialize = function(X = NULL, obs = NULL, var = NULL) {
     private$.X <- X
     private$.obs <- obs
     private$.var <- var
    }
  )
)