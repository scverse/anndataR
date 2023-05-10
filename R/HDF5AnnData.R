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
    .var_names = NULL
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_X, status=wip
        read_h5ad_element(private$.h5obj, "/X")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_X, status=wip
        value <- private$.validate_matrix(value, "X")
        write_h5ad_element(value, private$.h5obj, "/X")
      }
    },
    #' @field layers The layers slot. Must be NULL or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_layers, status=wip
        read_h5ad_element(private$.h5obj, "layers")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_layers, status=wip
        value <- private$.validate_layers(value)
        write_h5ad_element(value, private$.h5obj, "layers")
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs, status=wip
        read_h5ad_element(private$.h5obj, "/obs")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs, status=wip
        value <- private$.validate_obsvar_dataframe(value, "obs",
                                                    check_nrow = TRUE)
        write_h5ad_element(value, private$.h5obj, "/obs")
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var, status=wip
        read_h5ad_element(private$.h5obj, "/var")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var, status=wip
        value <- private$.validate_obsvar_dataframe(value, "var", 
                                                    check_nrow = TRUE)
        write_h5ad_element(value, private$.h5obj, "/var")
      }
    },
    #' @field obs_names Names of observations
    obs_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs_names, status=wip
        # obs names are cached to avoid reading all of obs whenever they are
        # accessed
        if (is.null(private$.obs_names)) {
          private$.obs_names <- rownames(self$obs)
        }
        private$.obs_names
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs_names, status=wip
        value <- private$.validate_obsvar_names(value, "obs")
        obs <- self$obs
        rownames(obs) <- value
        self$obs <- obs
        private$.obs_names <- value
      }
    },
    #' @field var_names Names of variables
    var_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var_names, status=wip
        # var names are cached to avoid reading all of var whenever they are
        # accessed
        if (is.null(private$.var_names)) {
          private$.var_names <- rownames(self$var)
        }
        private$.var_names
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var_names, status=wip
        value <- private$.validate_obsvar_names(value, "var")
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
    #' @param file The filename (character) of the `.h5ad` file. Alternatively,
    #' you can also provide an object created by `[rhdf5::H5Fopen()]`. If this
    #' file does not exist yet, you must provide values for `X`, `obs`, `var` to
    #' create a new AnnData with.
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
    #'   the variable dimension of the AnnData object. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    initialize = function(file, X, obs, var, obs_names, var_names, layers) {
      if (!inherits(file, "H5IdComponent") && !is.character(file)) {
        stop(
          "Argument 'file' should be a character path or an ",
          "H5IdComponent created by `rhdf5::H5Fopen()`."
        )
      }
      # check if file already exists
      create_new <- is.character(file) && !file.exists(file)

      if (create_new) {
        # create a new h5ad from scratch
        if (missing(X)) X <- NULL
        if (missing(obs)) stop("When creating a new .h5ad file, `obs` must be defined.")
        if (missing(var)) stop("When creating a new .h5ad file, `var` must be defined.")
        if (missing(obs_names)) obs_names <- NULL
        if (missing(var_names)) var_names <- NULL
        if (missing(layers)) layers <- NULL

        # create new H5AD from scratch
        private$.h5obj <- rhdf5::H5Fcreate(file)

        # create encoding attributes
        write_h5ad_encoding(private$.h5obj, "/", "anndata", "0.1.0")

        # write obs and var first, because these are used by other validators
        self$obs <- obs
        self$var <- var

        # write other slots later
        self$obs_names <- obs_names
        self$var_names <- var_names
        self$X <- X
        self$layers <- layers

      } else {
        # check if other arguments are defined
        slots_missing <- missing(X) && missing(obs) && missing(var) &&
          missing(obs_names) && missing(var_names) && missing(layers)
        if (!slots_missing) {
          stop(
            "Failed to create HDF5AnnData object. ",
            "Reason: If the file provided already exists, ",
            "then the other arguments (X, obs, var, ...) ",
            "should not be passed any value."
          )
        }

        # open h5 file if this isn't the already
        if (is.character(file)) {
          file <- rhdf5::H5Fopen(file)
        }

        # check for attributes
        attrs <- rhdf5::h5readAttributes(file, "/")

        if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
          stop(
            "H5AD encoding information is missing. ",
            "This file may have been created with Python anndata<0.8.0."
          )
        }

        # store file
        private$.h5obj <- file
      }
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
