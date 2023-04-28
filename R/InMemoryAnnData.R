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
InMemoryAnnData <- R6::R6Class("InMemoryAnnData",
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .obs = NULL,
    .var = NULL,
    .obs_names = NULL,
    .var_names = NULL,

    #' @description validate a matrix (.X or .layers[...])
    #' @param mat A matrix to validate
    #' @param label Must be `"X"` or `"layer[[...]]"` where `...` is the name of a layer.
    .validate_matrix = function(mat, label) {
      if (!is.null(mat)) {
        if (nrow(mat) != nrow(self$obs))
          stop("nrow(", label, ") should be the same as nrow(obs)")
        if (ncol(mat) != nrow(self$var))
          stop("ncol(", label, ") should be the same as nrow(var)")
      
        if (!is.null(rownames(mat))) {
          warning("rownames(", label, ") should be NULL, removing them from the matrix")
          rownames(mat) <- NULL
        }
      
        if (!is.null(colnames(mat))) {
          warning("colnames(", label, ") should be NULL, removing them from the matrix")
          colnames(mat) <- NULL
        }
      }

      mat
    },

    #' @description validate an obs or a var data frame
    #' @param df A data frame to validate. Should be an obs or a var.
    #' @param label Must be `"obs"` or `"var"`
    .validate_obsvar_dataframe = function(df, label) {
      if (is.null(df)) stop(label, " should be a data frame")
      if (.row_names_info(df) > 0) {
        warning(label, " should not have any dimnames, removing them from the matrix")
        rownames(df) <- NULL
      }
      df
    },

    #' @description validate an obs_names or a var_names vector
    #' @param names A vector to validate
    #' @param label Must be `"obs"` or `"var"`
    .validate_obsvar_names = function(names, label) {
      if (!is.null(names)) {
        if (length(names) != nrow(self[[label]])) stop("length(", label, "_names) should be the same as nrow(", label, ")")
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
        private$.X <- self$.validate_matrix(value, "X")
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- self$.validate_obsvar_dataframe(value, "obs")
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        private$.var <- self$.validate_obsvar_dataframe(value, "var")
      }
    },
    #' @field obs_names The obs_names slot
    obs_names = function(value) {
      if (missing(value)) {
        private$.obs_names
      } else {
        private$.obs_names <- self$.validate_obsvar_names(value, "obs")
      }
    },
    #' @field var_names The var_names slot
    var_names = function(value) {
      if (missing(value)) {
        private$.var_names
      } else {
        private$.var_names <- self$.validate_obsvar_names(value, "var")
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #' Inherits from AbstractAnnData
    #' 
    #' @param X The X slot
    #' @param obs The obs slot
    #' @param var The var slot
    #' @param obs_names The obs_names slot
    #' @param var_names The var_names slot
    initialize = function(X = NULL, obs, var, obs_names = NULL, var_names = NULL) {
      # check obs and var first, because these objects are used by other validators
      private$.obs <- private$.validate_obsvar_dataframe(obs, "obs")
      private$.var <- private$.validate_obsvar_dataframe(var, "var")

      # then check other values
      private$.obs_names <- private$.validate_obsvar_names(obs_names, "obs")
      private$.var_names <- private$.validate_obsvar_names(var_names, "var")
      private$.X <- private$.validate_matrix(X, "X")
    }
  )
)
