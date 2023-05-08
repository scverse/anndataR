#' @title HDF5AnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
HDF5AnnData <- R6::R6Class("HDF5AnnData", # nolint
  inherit = AbstractAnnData,
  private = list(
    .h5obj = NULL,
    .n_obs = NULL,
    .n_vars = NULL,
    .obs_names = NULL,
    .var_names = NULL,

    #' @description validate a value matches the observations dimension
    .validate_n_obs = function(value) {
      if (nrow(value) != self$n_obs()) {
        stop("Dimensions of value does not match the number of observations")
      }
    },

    #' @description validate a value matches the variables dimension
    .validate_n_vars = function(value) {
      if (nrow(value) != self$n_vars()) {
        stop("Dimensions of value does not match the number of variables")
      }
    },

    #' @description validate a value matches the AnnData shape
    .validate_shape = function(value) {
      if (!identical(dim(value), self$shape())) {
        stop("Dimensions of value does not match the object shape")
      }
    }
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        read_h5ad_element(private$.h5obj, "/X")
      } else {
        private$.validate_shape(value)
        write_h5ad_element(value, private$.h5obj, "/X")
      }
    },
    #' @field layers The layers slot. Must be NULL or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        read_h5ad_element(private$.h5obj, "layers")
      } else {
        write_h5ad_element(value, private$.h5obj, "layers")
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        read_h5ad_element(private$.h5obj, "/obs")
      } else {
        write_h5ad_element(value, private$.h5obj, "/obs")
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        read_h5ad_element(private$.h5obj, "/var")
      } else {
        write_h5ad_element(value, private$.h5obj, "/var")
      }
    },
    #' @field obs_names Names of observations
    obs_names = function(value) {
      if (missing(value)) {
        # obs names are cached to avoid reading all of obs whenever they are
        # accessed
        if (is.null(private$.obs_names)) {
          private$.obs_names <- rownames(self$obs)
        }
        private$.obs_names
      } else {
        obs <- self$obs
        rownames(obs) <- value
        self$obs <- obs
        private$.obs_names <- value
      }
    },
    #' @field var_names Names of variables
    var_names = function(value) {
      if (missing(value)) {
        # var names are cached to avoid reading all of var whenever they are
        # accessed
        if (is.null(private$.var_names)) {
          private$.var_names <- rownames(self$var)
        }
        private$.var_names
      } else {
        var <- self$var
        rownames(var) <- value
        self$var <- var
        private$.var_names <- value
      }
    }
  ),
  public = list(
    #' @description HDF5AnnData constructor
    #'
    #' @param h5obj The rhdf5 object
    initialize = function(h5obj) {
      attrs <- rhdf5::h5readAttributes(h5obj, "/")

      if (is.character(h5obj)) {
        h5obj <- path.expand(h5obj)
        if (!file.exists(h5obj)) {
          stop("Path to H5AD not found: ", h5obj)
        }
      }

      if (!("encoding-type") %in% names(attrs) ||
        !("encoding-version" %in% names(attrs))) {
        stop(
          "H5AD files without encodings are not supported ",
          "(this file may have been created with Python anndata prior to v0.8.0)"
        )
      }

      private$.h5obj <- h5obj
    },

    #' @description Number of observations in the AnnData object
    n_obs = function() {
      if (is.null(private$.n_obs)) {
        private$.n_obs <- nrow(self$obs)
      }
      private$.n_obs
    },

    #' @description Number of variables in the AnnData object
    n_vars = function() {
      if (is.null(private$.n_vars)) {
        private$.n_vars <- nrow(self$var)
      }
      private$.n_vars
    }
  )
)
