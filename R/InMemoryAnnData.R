#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#' 
#' @importFrom Matrix as.matrix
#'
#' @examples
#' ad <- InMemoryAnnData$new(
#'     X = matrix(1:5, 3L, 5L),
#'     obs = data.frame(cell = 1:3, row.names = LETTERS[1:3]),
#'     var = data.frame(gene = 1:5, row.names = letters[1:5])
#' )
#' ad
#'
#' @export
InMemoryAnnData <- R6::R6Class("InMemoryAnnData",
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .layers = NULL,
    .obs = NULL,
    .var = NULL,
    .obs_names = NULL,
    .var_names = NULL,

    #' @description validate a matrix (.X or .layers[...])
    .validate_matrix = function(mat, obj_name) {
      if (!is.null(mat)) {
        if (nrow(mat) != nrow(self$obs))
          stop("nrow(", obj_name, ") should be the same as nrow(obs)")
        if (ncol(mat) != nrow(self$var))
          stop("ncol(", obj_name, ") should be the same as nrow(var)")
      
        if (!is.null(rownames(mat))) {
          warning("rownames(", obj_name, ") should be NULL, removing them from the matrix")
          rownames(mat) <- NULL
        }
      
        if (!is.null(colnames(mat))) {
          warning("colnames(", obj_name, ") should be NULL, removing them from the matrix")
          colnames(mat) <- NULL
        }
      }

      mat
    },
    #' @description validate layers
    .validate_layers = function(layers) {
      if (is.null(layers)) return(layers)

      ## layers and names
      layer_names <- names(layers)
      if (!is.list(layers) || is.null(layer_names)) {
        stop("'layers' must must be a named list")
      }
      if (any(!nzchar(layer_names))) {
        stop("all 'layers' elements must have non-trivial names")
      }

      ## layer elements
      for (layer in layer_names) {
        layer_name <- paste0("layers[[", layer, "]]")
        private$.validate_matrix(layers[[layer]], layer_name)
      }

      layers
    },
    #' @description validate an obs data frame
    .validate_obs = function(obs) {
      if (is.null(obs)) stop("obs should be a data frame")
      if (.row_names_info(obs) > 0) {
        warning("obs should not have any dimnames, removing them from the matrix")
        rownames(obs) <- NULL
      }
      obs
    },

    #' @description validate a var data frame
    .validate_var = function(var) {
      if (is.null(var)) stop("var should be a data frame")
      if (.row_names_info(var) > 0) {
        warning("var should not have any rownames, removing them from the matrix")
        rownames(var) <- NULL
      }
      var
    },

    #' @description validate an obs_names vector
    .validate_obs_names = function(obs_names) {
      if (!is.null(obs_names)) {
        if (length(obs_names) != nrow(self$obs)) stop("length(obs_names) should be the same as nrow(obs)")
      }
      obs_names
    },

    #' @description validate a var_names vector
    .validate_var_names = function(var_names) {
      if (!is.null(var_names)) {
        if (length(var_names) != nrow(self$var)) stop("length(var_names) should be the same as nrow(var)")
      }
      var_names
    }
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- self$.validate_matrix(value, "X")
      }
    },
    #' @field layers The layers slot. Must be NULL or a named list with with all elements having the dimensions consistent with `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        private$.layers
      } else {
        private$.layers <- self$.validate_layers(value)
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- self$.validate_obs(value)
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        private$.var <- self$.validate_var(value)
      }
    },
    #' @field obs_names The obs_names slot
    obs_names = function(value) {
      if (missing(value)) {
        private$.obs_names
      } else {
        private$.obs_names <- self$.validate_obs_names(value)
      }
    },
    #' @field var_names The var_names slot
    var_names = function(value) {
      if (missing(value)) {
        private$.var_names
      } else {
        private$.var_names <- self$.validate_var_names(value)
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #' Inherits from AbstractAnnData
    #' 
    #' @param X The X slot
    #' @param layers The layers slot. Either NULL  or a named list, where all elements have dimensions consistent with `obs` and `var`.
    #' @param obs The obs slot
    #' @param var The var slot
    #' @param obs_names The obs_names slot
    #' @param var_names The var_names slot
    initialize = function(X = NULL, obs, var, obs_names = NULL, var_names = NULL, layers = NULL) {
      # check obs and var first
      obs <- private$.validate_obs(obs)
      var <- private$.validate_var(var)
      private$.obs <- obs
      private$.var <- var

      # check inputs
      X <- private$.validate_matrix(X, "X")
      layers <- private$.validate_layers(layers)
      obs_names <- private$.validate_obs_names(obs_names)
      var_names <- private$.validate_var_names(var_names)

      # store results
      private$.X <- X
      private$.layers <- layers
      private$.obs_names <- obs_names
      private$.var_names <- var_names
    }
  )
)
