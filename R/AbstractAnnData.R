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
      if (missing(value)) {
        names <- names(self$obs)
        if (is.null(names)) {
          names <- character()
        }
        names
      } else {
        .abstract_function()
      }
    },
    #' @field var_names The var_names slot
    var_names = function(value) {
      if (missing(value)) {
        names <- names(self$var)
        if (is.null(names)) {
          names <- character()
        }
        names
      } else {
        .abstract_function()
      }
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
        "dim: ", self$dim()[1], " x ", self$dim()[2], "\n",
        "X: ", X_info, "\n",
        pretty_print("obs", self$obs_names), "\n",
        pretty_print("var", self$var_names), "\n",
        sep = ""
      )
    },
    #' @description Dimensions of the AnnData object
    dim = function() {
      n_row <- 0L
      if (!is.null(self$obs)) {
        n_row <- NROW(self$obs)
      } else if (!is.null(self$X)) {
        n_row <- NROW(self$X)
      }

      n_col <- 0L
      if (!is.null(self$var)) {
        n_col <- NROW(self$var)
      } else if (!is.null(self$X)) {
        n_col <- NCOL(self$X)
      }

      c(n_row, n_col)
    }
  )
)
