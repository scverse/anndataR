#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#'
#' @importFrom Matrix as.matrix
#'
#' @examples
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5),
#'   obs_names = LETTERS[1:3],
#'   var_names = letters[1:5]
#' )
#' ad
#'
#' @export
InMemoryAnnData <- R6::R6Class("InMemoryAnnData", # nolint
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .layers = NULL,
    .obs = NULL,
    .var = NULL,
    .obs_names = NULL,
    .var_names = NULL,

    # @description validate a matrix (.X or .layers[...])
    # @param mat A matrix to validate
    # @param label Must be `"X"` or `"layer[[...]]"` where `...` is the name
    #   of a layer.
    .validate_matrix = function(mat, label) {
      if (!is.null(mat)) {
        if (nrow(mat) != nrow(self$obs)) {
          stop("nrow(", label, ") should be the same as nrow(obs)")
        }
        if (ncol(mat) != nrow(self$var)) {
          stop("ncol(", label, ") should be the same as nrow(var)")
        }

        if (!is.null(rownames(mat))) {
          warning(
            "rownames(", label, ") should be NULL, removing them from ",
            "the matrix"
          )
          rownames(mat) <- NULL
        }

        if (!is.null(colnames(mat))) {
          warning(
            "colnames(", label, ") should be NULL, removing them from ",
            "the matrix"
          )
          colnames(mat) <- NULL
        }
      }

      mat
    },
    # @description validate layers
    .validate_layers = function(layers) {
      if (is.null(layers)) {
        return(layers)
      }

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
    # @description validate an obs or a var data frame
    # @param df A data frame to validate. Should be an obs or a var.
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_dataframe = function(df, label) {
      if (is.null(df)) stop(label, " should be a data frame")
      if (.row_names_info(df) > 0) {
        warning(
          label, " should not have any dimnames, removing them from ",
          "the matrix"
        )
        rownames(df) <- NULL
      }
      df
    },

    # @description validate an obs_names or a var_names vector
    # @param names A vector to validate
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_names = function(names, label) {
      if (!is.null(names)) {
        stopifnot(length(names) == nrow(self[[label]]))
      }
      names
    }
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- private$.validate_matrix(value, "X")
        self
      }
    },
    #' @field layers The layers slot. Must be NULL or a named list with with
    #' all elements having the dimensions consistent with `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        private$.layers
      } else {
        private$.layers <- private$.validate_layers(value)
        self
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- private$.validate_obsvar_dataframe(value, "obs")
        self
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        private$.var <- private$.validate_obsvar_dataframe(value, "var")
        self
      }
    },
    #' @field obs_names The obs_names slot
    obs_names = function(value) {
      if (missing(value)) {
        private$.obs_names
      } else {
        private$.obs_names <- private$.validate_obsvar_names(value, "obs")
        self
      }
    },
    #' @field var_names The var_names slot
    var_names = function(value) {
      if (missing(value)) {
        private$.var_names
      } else {
        private$.var_names <- private$.validate_obsvar_names(value, "var")
        self
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #' Inherits from AbstractAnnData
    #'
    #' @param X The X slot
    #' @param layers The layers slot. Either NULL  or a named list, where all
    #'   elements have dimensions consistent with `obs` and `var`.
    #' @param obs The obs slot
    #' @param var The var slot
    #' @param obs_names The obs_names slot
    #' @param var_names The var_names slot
    initialize = function(X = NULL, # nolint
                          obs,
                          var,
                          obs_names = NULL,
                          var_names = NULL,
                          layers = NULL) {
      # check obs and var first, because these objects are used by
      # other validators
      private$.obs <- private$.validate_obsvar_dataframe(obs, "obs")
      private$.var <- private$.validate_obsvar_dataframe(var, "var")

      # then check other values
      private$.obs_names <- private$.validate_obsvar_names(obs_names, "obs")
      private$.var_names <- private$.validate_obsvar_names(var_names, "var")
      private$.X <- private$.validate_matrix(X, "X")
      private$.layers <- private$.validate_layers(layers)
    }
  )
)
