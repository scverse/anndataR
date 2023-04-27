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
    .var = NULL
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        private$.X <- .inmemoryanndata_check_matrix(value, self$obs, self$var)
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        private$.obs
      } else {
        # TODO: add checks
        private$.obs <- value
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        private$.var
      } else {
        # TODO: add checks
        private$.var <- value
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #' Inherits from AbstractAnnData
    #' 
    #' @param X A #observations × #variables data matrix. A view of the data is used if the data type matches, otherwise, a copy is made.
    #' @param obs Key-indexed one-dimensional observations annotation of length #observations.
    #' @param var Key-indexed one-dimensional variables annotation of length #variables.
    #' @param shape Shape list (#observations, #variables). Can only be provided if `X` is `NULL`.
    initialize = function(X = NULL, obs = NULL, var = NULL, shape = NULL) {
      
      # check nrow size
      nrow <- nrow(X)
      if (is.null(nrow)) nrow <- nrow(obs)
      if (is.null(nrow) && !is.null(shape)) nrow <- shape[[1]]
      if (is.null(nrow)) stop("If $X, $obs and $var are NULL, shape should be set to the dimensions of the AnnData.")

      # check ncol size
      ncol <- ncol(X)
      if (is.null(ncol)) ncol <- nrow(var)
      if (is.null(ncol) && !is.null(shape)) ncol <- shape[[2]]
      if (is.null(ncol)) stop("If $X, $obs and $var are NULL, shape should be set to the dimensions of the AnnData.")

      # check for obs names
      obs_names <- rownames(X)
      if (is.null(obs_names)) obs_names <- rownames(obs)
      if (is.null(obs_names)) obs_names <- as.character(seq_len(nrow))

      # check for var names
      var_names <- colnames(X)
      if (is.null(var_names)) var_names <- rownames(var)
      if (is.null(var_names)) var_names <- as.character(seq_len(ncol))

      # construct obs if it doesn't exist
      if (is.null(obs)) {
        obs <- data.frame(row.names = obs_names)
      }
      if (!is.data.frame(obs)) stop("$obs should be a data frame.")
      if (is.null(rownames(obs))) {
        rownames(obs) <- obs_names
      }

      # construct var if it doesn't exist
      if (is.null(var)) {
        var <- data.frame(row.names = var_names)
      }
      if (!is.data.frame(var)) stop("$var should be a data frame.")
      if (is.null(rownames(var))) {
        rownames(var) <- var_names
      }

      # check matrix
      X <- .inmemoryanndata_check_matrix(X, obs, var)

      private$.X <- X
      private$.obs <- obs
      private$.var <- var
    }
  )
)

.inmemoryanndata_check_matrix <- function(mat, obs, var) {
  if (!is.null(mat)) {
    if (nrow(mat) != nrow(obs)) stop("$X should have the same number of rows as $obs has rows.")
    if (ncol(mat) != nrow(var)) stop("$X should have the same number of columns as $var has rows.")

    rownames(mat) <- rownames(obs)
    colnames(mat) <- rownames(var)
  }

  mat
}