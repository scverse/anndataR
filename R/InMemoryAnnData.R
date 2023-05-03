#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#'
#' @importFrom Matrix as.matrix
#'
#' @examples
#' ## complete example
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5),
#'   obs_names = LETTERS[1:3],
#'   var_names = letters[1:5]
#' )
#' ad
#'
#' ## minimal example -- no observations or variables
#' ad <- InMemoryAnnData$new(
#'   obs = data.frame(),
#'   var = data.frame()
#' )
#' ad
#'
#' ## number of  observations or variables determined by `obs` and `var`; no
#' ## `X` or `layers`; `obs_names` and `var_names` determined automatically
#' ad <- InMemoryAnnData$new(
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5)
#' )
#' ad
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

    # @description `.validate_matrix()` checks that dimensions are
    #   consistent with `obs` and `var`, and removes dimnames if
    #   present.
    # @param mat A matrix to validate
    # @param label Must be `"X"` or `"layer[[...]]"` where `...` is
    #   the name of a layer.
    .validate_matrix = function(mat, label) {
      if (!is.null(mat)) {
        if (nrow(mat) != nrow(self$obs)) {
          stop("nrow(", label, ") should be the same as nrow(obs)")
        }
        if (ncol(mat) != nrow(self$var)) {
          stop("ncol(", label, ") should be the same as nrow(var)")
        }

        if (!is.null(rownames(mat))) {
          warning(wrap_message(
            "rownames(", label, ") should be NULL, removing them from the ",
            "matrix"
          ))
          rownames(mat) <- NULL
        }

        if (!is.null(colnames(mat))) {
          warning(wrap_message(
            "colnames(", label, ") should be NULL, removing them from the ",
            "matrix"
          ))
          colnames(mat) <- NULL
        }
      }

      mat
    },

    # @description `.validate_layers()` checks for named lists and
    #   correct dimensions on elements.
    # @param layers A named list of 0 or more matrix elements with
    #   dimensions consistent with `obs` and `var`.
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

    # @description `.validate_obsvar_dataframe()` checks that the
    #   object is a data.frame and removes explicit dimnames.
    # @param df A data frame to validate. Should be an obs or a var.
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_dataframe = function(df, label) {
      if (is.null(df)) stop(label, " should be a data frame")
      if (.row_names_info(df) > 0) {
        warning(wrap_message(
          "'", label, "' should not have any dimnames, removing them from ",
          "the matrix"
        ))
        rownames(df) <- NULL
      }
      df
    },

    # @description `.validate_obsvar_names()` checks that `*_names()`
    #   are NULL or consistent with the dimensions of `obs` or `var`.
    # @param names A vector to validate
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_names = function(names, label) {
      if (!is.null(names) && length(names) != nrow(self[[label]])) {
        stop(wrap_message(
          "length(", label, "_names) should be the same as ",
          "nrow(", label, ")"
        ))
      }
      names
    }
  ),
  active = list(
    #' @field X NULL or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and
    #"   `var`.
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- private$.validate_matrix(value, "X")
        self
      }
    },
    #' @field layers NULL or a named list with all elements having the
    #'   dimensions consistent with `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        private$.layers
      } else {
        private$.layers <- private$.validate_layers(value)
        self
      }
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- private$.validate_obsvar_dataframe(value, "obs")
        self
      }
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        private$.var <- private$.validate_obsvar_dataframe(value, "var")
        self
      }
    },
    #' @field obs_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into
    #'   the observation dimension of the AnnData object. For
    #'   compatibility with *R* representations, `obs_names` should be a
    #'   character vector.
    obs_names = function(value) {
      if (missing(value)) {
        private$.obs_names
      } else {
        private$.obs_names <- private$.validate_obsvar_names(value, "obs")
        self
      }
    },
    #' @field var_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the AnnData object.. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
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
    #'   Inherits from AbstractAnnData.
    #' @param X Either NULL or a observation x variable matrix with
    #'   dimensions consistent with `obs` and `var`.
    #' @param layers Either NULL or a named list, where each element
    #'   is an observation x variable matrix with dimensions consistent
    #'   with `obs` and `var`.
    #' @param obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    #' @param var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    #' @param obs_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into
    #'   the observation dimension of the AnnData object. For
    #'   compatibility with *R* representations, `obs_names` should be a
    #'   character vector.
    #' @param var_names Either NULL or a vector of unique identifers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the AnnData object.. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    initialize = function(
      X = NULL, obs, var, obs_names = NULL, var_names = NULL, layers = NULL
    ) {
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

#' Convert an AnnData object to an InMemoryAnnData object
#'
#' This function takes an AnnData object and converts it to an InMemoryAnnData
#' object, loading all fields into memory.
#'
#' @param adata An AnnData object to be converted to InMemoryAnnData.
#'
#' @return An InMemoryAnnData object with the same data as the input AnnData
#'   object.
#'
#' @export
#'
#' @examples
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5),
#'   obs_names = LETTERS[1:3],
#'   var_names = letters[1:5]
#' )
#' to_InMemory(ad)
to_InMemory <- function(adata) { # nolint
  InMemoryAnnData$new(
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    obs_names = adata$obs_names,
    var_names = adata$var_names,
    layers = adata$layers
  )
}
