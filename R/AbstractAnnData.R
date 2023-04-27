.abstract_function <- function() {
  stop("This function is not implemented yet.")
}

#' @title AbstractAnnData
#'
#' @description
#' Abstract class representing an anndata object. Defines the interface.
#' @importFrom R6 R6Class
AbstractAnnData <- R6::R6Class("AbstractAnnData",
  active = list(
    #' @field X The X slot
    X = function(value) {
      .abstract_function()
    },
    #' @field obs The obs slot
    obs = function(value) {
      .abstract_function()
    },
    #' @field var The var slot
    var = function(value) {
      .abstract_function()
    }
  )
)
