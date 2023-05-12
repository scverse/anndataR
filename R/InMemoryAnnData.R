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
#'   X = matrix(1:15, 3L, 5L),
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
#' ## minimum example
#' ad <- InMemoryAnnData$new(
#'   obs_names = letters[1:10],
#'   var_names = LETTERS[1:5]
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
    .var_names = NULL
  ),
  active = list(
    #' @field X NULL or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and
    #'   `var`.
    X = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_X, status=done
        private$.X
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_X, status=done
        private$.X <- private$.validate_matrix(value, "X")
        self
      }
    },
    #' @field layers NULL or a named list with all elements having the
    #'   dimensions consistent with `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_layers, status=done
        private$.layers
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_layers, status=done
        private$.layers <- private$.validate_layers(value)
        self
      }
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_obs, status=done
        private$.obs
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_obs, status=done
        private$.obs <- private$.validate_obsvar_dataframe(value, "obs")
        self
      }
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_var, status=done
        private$.var
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_var, status=done
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
        # trackstatus: class=InMemoryAnnData, feature=get_obs_names, status=done
        private$.obs_names
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_obs_names, status=done
        private$.obs_names <- private$.validate_obsvar_names(value, "obs")
        self
      }
    },
    #' @field var_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the AnnData object. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    var_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_var_names, status=done
        private$.var_names
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_var_names, status=done
        private$.var_names <- private$.validate_obsvar_names(value, "var")
        self
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #'   Inherits from AbstractAnnData.
    #' @param X Either `NULL` or a observation × variable matrix with
    #'   dimensions consistent with `obs` and `var`.
    #' @param layers Either `NULL` or a named list, where each element
    #'   is an observation × variable matrix with dimensions consistent
    #'   with `obs` and `var`.
    #' @param obs Either `NULL` or a `data.frame` with columns containing information
    #'   about observations. If `NULL`, an `n_obs`×0 data frame will automatically
    #'   be generated.
    #' @param var Either `NULL` or a `data.frame` with columns containing information
    #'   about variables. If `NULL`, an `n_vars`×0 data frame will automatically
    #'   be generated.
    #' @param obs_names A vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into
    #'   the observation dimension of the AnnData object. The length of
    #'   the `obs_names` defines the observation dimension of the AnnData
    #'   object.
    #' @param var_names A vector of unique identifers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the AnnData object. The length of
    #'   the `var_names` defines the variable dimension of the AnnData
    #'   object.
    initialize = function(X = NULL, obs = NULL, var = NULL, obs_names, var_names, layers = NULL) {
      # check obs_names and var_names first, because these objects are used by
      # other validators
      private$.obs_names <- private$.validate_obsvar_names(obs_names, "obs", check_length = FALSE)
      private$.var_names <- private$.validate_obsvar_names(var_names, "var", check_length = FALSE)

      # then check other values
      private$.obs <- private$.validate_obsvar_dataframe(obs, "obs")
      private$.var <- private$.validate_obsvar_dataframe(var, "var")
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
