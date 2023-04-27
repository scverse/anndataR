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
    },
    #' @field obs_names The obs_names slot
    obs_names = function(value) {
      .abstract_function()
    },
    #' @field var_names The var_names slot
    var_names = function(value) {
      .abstract_function()
    }
  ),
  public = list(
    #' @description Print AnnData object
    #' @param ... optional arguments to print method.
    print = function(...) {
      X_info <- if (!is.null(self$X)) {
        class(self$X)[[1]]
      } else {
        NULL
      }
      cat(
        "class: ", class(self)[[1]], "\n",
        "dim: ", self$shape()[[1]], " x ", self$shape()[[2]], "\n",
        "X: ", X_info, "\n",
        pretty_print("obs", self$obs_names), "\n",
        pretty_print("var", self$var_names), "\n",
        sep = ""
      )
    },
    #' @description Dimensions of the AnnData object
    shape = function() {
      c(
        nrow(self$obs),
        nrow(self$var)
      )
    },
    #' @description Return a new AnnData object with all objects loaded into memory.
    to_inmemory = function() {
      # should probably be stored in a separate file
      InMemoryAnnData$new(
        X = self$X,
        obs = self$obs,
        var = self$var
      )
    }
  )
)
