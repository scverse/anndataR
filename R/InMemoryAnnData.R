#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#' 
#' @importFrom Matrix as.matrix
InMemoryAnnData <- R6::R6Class("InMemoryAnnData",
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .obs = NULL,
    .var = NULL
  ),
  active = list(

    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- value
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- value
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        private$.var <- value
      }
    }
  ),
  public = list(
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
