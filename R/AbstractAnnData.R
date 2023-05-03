.abstract_function <- function() {
  stop("This function is not implemented yet.")
}

#' @title AbstractAnnData
#'
#' @description
#'   Abstract [R6][R6::R6Class] class representing an AnnData
#'   object. Defines the interface.
#' @importFrom R6 R6Class
AbstractAnnData <- R6::R6Class("AbstractAnnData",
  active = list(
    #' @field X NULL or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and `var`.
    X = function(value) {
      .abstract_function()
    },
    #' @field layers The layers slot. Must be NULL or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      .abstract_function()
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    obs = function(value) {
      .abstract_function()
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      .abstract_function()
    },
    #' @field obs_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into
    #'   the observation dimension of the AnnData object. For
    #'   compatibility with *R* representations, `obs_names` should be a
    #'   character vector.
    obs_names = function(value) {
      .abstract_function()
    },
    #' @field var_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the AnnData object.. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    var_names = function(value) {
      .abstract_function()
    }
  ),
  public = list(
    #' @description Print a summary of the AnnData object. `print()`
    #'   methods should be implemented so that they are not
    #'   computationally expensive.
    #' @param ... Optional arguments to print method.
    print = function(...) {
      x_info <- if (!is.null(self$X)) {
        class(self$X)[[1]]
      } else {
        NULL
      }
      cat(
        "class: ", class(self)[[1]], "\n",
        "dim: ", self$n_obs(), " obs x ", self$n_vars(), " var\n",
        "X: ", x_info, "\n",
        pretty_print("layers", self$layers_keys()), "\n",
        pretty_print("obs", self$obs_keys()), "\n",
        pretty_print("var", self$var_keys()), "\n",
        sep = ""
      )
    },

    #' @description Dimensions (observations x variables) of the AnnData object.

    shape = function() {
      c(
        self$n_obs(),
        self$n_vars()
      )
    },
    #' @description Number of observations in the AnnData object.
    n_obs = function() {
      nrow(self$obs)
    },
    #' @description Number of variables in the AnnData object.
    n_vars = function() {
      nrow(self$var)
    },
    #' @description Keys ('column names') of `obs`.
    obs_keys = function() {
      names(self$obs)
    },
    #' @description Keys ('column names') of `var`.
    var_keys = function() {
      names(self$var)
    },
    #' @description Keys (element names) of `layers`.
    layers_keys = function() {
      names(self$layers)
    }
  )
)
