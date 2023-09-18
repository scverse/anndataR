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
    .obsm = NULL,
    .varm = NULL
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
        write_h5ad_element(value, private$.h5obj, "/layers")
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs, status=wip
        read_h5ad_element(private$.h5obj, "/obs", include_index = FALSE)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs, status=wip
        value <- private$.validate_obsvar_dataframe(value, "obs")
        write_h5ad_element(
          value,
          private$.h5obj,
          "/obs",
          index = self$obs_names
        )
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var, status=wip
        read_h5ad_element(private$.h5obj, "/var", include_index = FALSE)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var, status=wip
        value <- private$.validate_obsvar_dataframe(value, "var")
        write_h5ad_element(
          value,
          private$.h5obj,
          "/var",
          index = self$var_names
        )
      }
    },
    #' @field obs_names Names of observations
    obs_names = function(value) {
      # TODO: directly write to and read from /obs/_index
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs_names, status=wip
        # obs names are cached to avoid reading all of obs whenever they are
        # accessed
        if (is.null(private$.obs_names)) {
          private$.obs_names <- read_h5ad_data_frame_index(private$.h5obj, "obs")
        }
        private$.obs_names
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs_names, status=wip
        value <- private$.validate_obsvar_names(value, "obs")
        write_h5ad_data_frame_index(value, private$.h5obj, "obs", "_index")
        private$.obs_names <- value
      }
    },
    #' @field var_names Names of variables
    var_names = function(value) {
      # TODO: directly write to and read from /var/_index
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var_names, status=wip
        # var names are cached to avoid reading all of var whenever they are
        # accessed
        if (is.null(private$.var_names)) {
          private$.var_names <- read_h5ad_data_frame_index(private$.h5obj, "var")
        }
        private$.var_names
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var_names, status=wip
        value <- private$.validate_obsvar_names(value, "var")
        write_h5ad_data_frame_index(value, private$.h5obj, "var", "_index")
        private$.var_names <- value
      }
    },
    #' @field obsm
    obsm = function(value) {
      if (missing(value)) {
        if (is.null(private$.obsm)) {
          private$.obsm <- read_h5ad_element(private$.h5obj, "obsm")
        }
      } else {
        # TODO: validate obsm
        write_h5ad_element(value, private$.h5obj, "/obsm")
      }
    },
    #' @field varm
    varm = function(value) {
      if (missing(value)) {
        if (is.null(private$.varm)) {
          private$.varm <- read_h5ad_element(private$.h5obj, "varm")
        }
      } else {
        # TODO: validate varm
        write_h5ad_element(value, private$.h5obj, "/varm")
      }
    }
  ),
  public = list(
    #' @description HDF5AnnData constructor
    #'
    #' @param file The filename (character) of the `.h5ad` file. If this
    #'   file does not exist yet, `obs_names` and `var_names` must be provided.
    #' @param obs_names A vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into the
    #'   observation dimension of the AnnData object. The length of `obs_names`
    #'   defines the observation dimension of the AnnData object.
    #' @param var_names A vector of unique identifiers used to identify each row
    #'   of `var` and to act as an index into the variable dimension of the
    #'   AnnData object. The length of `var_names` defines the variable
    #'   dimension of the AnnData object.
    #' @param X Either `NULL` or a observation × variable matrix with
    #'   dimensions consistent with `obs` and `var`.
    #' @param layers Either `NULL` or a named list, where each element is an
    #'   observation × variable matrix with dimensions consistent with `obs` and
    #'   `var`.
    #' @param obs Either `NULL` or a `data.frame` with columns containing
    #'   information about observations. If `NULL`, an `n_obs`×0 data frame will
    #'   automatically be generated.
    #' @param var Either `NULL` or a `data.frame` with columns containing
    #'   information about variables. If `NULL`, an `n_vars`×0 data frame will
    #'   automatically be generated.
    #'
    #' @details
    #' The constructor creates a new HDF5 AnnData interface object. This can
    #' either be used to either connect to an existing `.h5ad` file or to
    #' create a new one. To create a new file both `obs_names` and `var_names`
    #' must be specified. In both cases, any additional slots provided will be
    #' set on the created object. This will cause data to be overwritten if the
    #' file already exists.
    initialize = function(file, obs_names = NULL, var_names = NULL, X = NULL,
                          obs = NULL, var = NULL, layers = NULL) {
      if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("The HDF5 interface requires the 'rhdf5' package to be installed")
      }

      if (!file.exists(file)) {
        # Check obs_names and var_names have been provided
        if (is.null(obs_names)) {
          stop("When creating a new .h5ad file, `obs_names` must be defined.")
        }
        if (is.null(var_names)) {
          stop("When creating a new .h5ad file, `var_names` must be defined.")
        }

        # Create an empty H5AD using the provided obs/var names
        write_empty_h5ad(file, obs_names, var_names)

        # Set private object slots
        private$.h5obj <- file
        private$.n_obs <- length(obs_names)
        private$.n_vars <- length(var_names)
        private$.obs_names <- obs_names
        private$.var_names <- var_names
      } else {
        # Check the file is a valid H5AD
        attrs <- rhdf5::h5readAttributes(file, "/")

        if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
          stop(
            "H5AD encoding information is missing. ",
            "This file may have been created with Python anndata<0.8.0."
          )
        }

        # Set the file path
        private$.h5obj <- file

        # If obs or var names have been provided update those
        if (!is.null(obs_names)) {
          self$obs_names <- obs_names
        }

        if (!is.null(var_names)) {
          self$var_names <- var_names
        }
      }

      # Update remaining slots
      if (!is.null(X)) {
        self$X <- X
      }

      if (!is.null(obs)) {
        self$obs <- obs
      }

      if (!is.null(var)) {
        self$var <- var
      }

      if (!is.null(layers)) {
        self$layers <- layers
      }
    },

    #' @description Number of observations in the AnnData object
    n_obs = function() {
      if (is.null(private$.n_obs)) {
        private$.n_obs <- length(self$obs_names)
      }
      private$.n_obs
    },

    #' @description Number of variables in the AnnData object
    n_vars = function() {
      if (is.null(private$.n_vars)) {
        private$.n_vars <- length(self$var_names)
      }
      private$.n_vars
    }
  )
)

#' Convert an AnnData object to an HDF5AnnData object
#'
#' This function takes an AnnData object and converts it to an HDF5AnnData
#' object, loading all fields into memory.
#'
#' @param adata An AnnData object to be converted to HDF5AnnData.
#' @param file The filename (character) of the `.h5ad` file.
#'
#' @return An HDF5AnnData object with the same data as the input AnnData
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
#' to_HDF5AnnData(ad, "test.h5ad")
#' # remove file
#' file.remove("test.h5ad")
to_HDF5AnnData <- function(adata, file) { # nolint
  stopifnot(
    inherits(adata, "AbstractAnnData")
  )
  HDF5AnnData$new(
    file = file,
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    obs_names = adata$obs_names,
    var_names = adata$var_names,
    layers = adata$layers
  )
}
