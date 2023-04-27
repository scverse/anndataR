#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#' 
#' @importFrom Matrix as.matrix
InMemoryAnnData <- R6::R6Class("InMemoryAnnData",
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .obs = NULL,
    .var = NULL,
    .obs_names = NULL,
    .var_names = NULL
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- .inmemoryanndata_check_matrix(value, self$obs, self$var, "X")
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        private$.obs <- .inmemoryanndata_check_obs(value)
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        private$.var <- .inmemoryanndata_check_var(value)
      }
    },
    #' @field obs_names The obs_names slot
    obs_names = function(value) {
      if (missing(value)) {
        private$.obs_names
      } else {
        private$.obs_names <- .inmemoryanndata_check_obs_names(value)
      }
    },
    #' @field var_names The var_names slot
    var_names = function(value) {
      if (missing(value)) {
        private$.var_names
      } else {
        private$.var_names <- .inmemoryanndata_check_var_names(value)
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
      # check inputs
      X <- .inmemoryanndata_check_matrix(X, obs, var, "X")
      obs <- .inmemoryanndata_check_obs(obs)
      var <- .inmemoryanndata_check_var(var)
      obs_names <- .inmemoryanndata_check_obs_names(obs_names, obs)
      var_names <- .inmemoryanndata_check_var_names(var_names, var)

      # store results
      private$.X <- X
      private$.obs <- obs
      private$.var <- var
      private$.obs_names <- obs_names
      private$.var_names <- var_names
    }
  )
)

.inmemoryanndata_check_matrix <- function(mat, obs, var, obj_name) {
  if (!is.null(mat)) {
    if (nrow(mat) != nrow(obs)) stop("nrow(", obj_name, ") should be the same as nrow(obs)")
    if (ncol(mat) != nrow(var)) stop("ncol(", obj_name, ") should be the same as nrow(var)")
  
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
}

.inmemoryanndata_check_obs <- function(obs) {
  if (is.null(obs)) stop("obs should be a data frame")
  if (.row_names_info(obs) >= 0) {
    warning("obs should not have any dimnames, removing them from the matrix")
    rownames(obs) <- NULL
  }
  obs
}

.inmemoryanndata_check_var <- function(var) {
  if (is.null(var)) stop("var should be a data frame")
  if (.row_names_info(var) >= 0) {
    warning("var should not have any rownames, removing them from the matrix")
    rownames(var) <- NULL
  }
  var
}

.inmemoryanndata_check_obs_names <- function(obs_names, obs) {
  if (!is.null(obs_names)) {
    if (length(obs_names) != nrow(obs)) stop("length(obs_names) should be the same as nrow(obs)")
  }
  obs_names
}

.inmemoryanndata_check_var_names <- function(var_names, var) {
  if (!is.null(var_names)) {
    if (length(var_names) != nrow(var)) stop("length(var_names) should be the same as nrow(var)")
  }
  var_names
}