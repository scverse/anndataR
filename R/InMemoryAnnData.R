#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
InMemoryAnnData <- R6::R6Class("InMemoryAnnData",
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .obs = NULL,
    .var = NULL
  ),
  active = list(

    #' @description The X slot
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- value
      }
    },
    #' @description The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- value
      }
    },
    #' @description The var slot
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
