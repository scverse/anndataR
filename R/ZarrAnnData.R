#' @title ZarrAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
ZarrAnnData <- R6::R6Class("ZarrAnnData", # nolint
  inherit = AbstractAnnData,
  private = list(
    zarr_store = NULL,
    zarr_root = NULL,
    .n_obs = NULL,
    .n_vars = NULL,
    .obs_names = NULL,
    .var_names = NULL,
    .compression = NULL
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_X, status=done
        read_zarr_element(private$zarr_store, "/X")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_X, status=done
        value <- private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        write_zarr_element(value, private$zarr_store, "/X", private$.compression)
      }
    },
    #' @field layers The layers slot. Must be NULL or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_layers, status=done
        read_zarr_element(private$zarr_store, "layers")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_layers, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        write_zarr_element(value, private$zarr_store, "/layers", private$.compression)
      }
    },
    #' @field obsm The obsm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `obs`.
    obsm = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsm, status=done
        read_zarr_element(private$zarr_store, "obsm")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obsm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = rownames(self)
        )
        write_zarr_element(value, private$zarr_store, "/obsm")
      }
    },
    #' @field varm The varm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `var`.
    varm = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varm, status=done
        read_zarr_element(private$zarr_store, "varm")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_varm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = colnames(self)
        )
        write_zarr_element(value, private$zarr_store, "/varm")
      }
    },
    #' @field obsp The obsp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `obs`.
    obsp = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsp, status=done
        read_zarr_element(private$zarr_store, "obsp")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obsp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = rownames(self),
          expected_colnames = rownames(self)
        )
        write_zarr_element(value, private$zarr_store, "/obsp")
      }
    },
    #' @field varp The varp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `var`.
    varp = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varp, status=done
        read_zarr_element(private$zarr_store, "varp")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_varp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = colnames(self),
          expected_colnames = colnames(self)
        )
        write_zarr_element(value, private$zarr_store, "/varp")
      }
    },

    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs, status=done
        read_zarr_element(private$zarr_store, "/obs", include_index = FALSE)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs, status=done
        value <- private$.validate_obsvar_dataframe(value, "obs")
        write_zarr_element(
          value,
          private$zarr_store,
          "/obs",
          private$.compression,
          index = self$obs_names
        )
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var, status=done
        read_zarr_element(private$zarr_store, "/var", include_index = FALSE)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var, status=done
        value <- private$.validate_obsvar_dataframe(value, "var")
        write_zarr_element(
          value,
          private$zarr_store,
          "/var",
          index = self$var_names
        )
      }
    },
    #' @field obs_names Names of observations
    obs_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs_names, status=done
        # obs names are cached to avoid reading all of obs whenever they are
        # accessed
        if (is.null(private$.obs_names)) {
          private$.obs_names <- read_zarr_data_frame_index(private$zarr_store, "obs")
        }
        private$.obs_names
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs_names, status=done
        value <- private$.validate_obsvar_names(value, "obs")
        write_zarr_data_frame_index(value, private$zarr_store, "obs", private$.compression, "_index")
        private$.obs_names <- value
      }
    },
    #' @field var_names Names of variables
    var_names = function(value) {
      # TODO: directly write to and read from /var/_index
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var_names, status=done
        # var names are cached to avoid reading all of var whenever they are
        # accessed
        if (is.null(private$.var_names)) {
          private$.var_names <- read_zarr_data_frame_index(private$zarr_store, "var")
        }
        private$.var_names
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var_names, status=done
        value <- private$.validate_obsvar_names(value, "var")
        write_zarr_data_frame_index(value, private$zarr_store, "var", private$.compression, "_index")
        private$.var_names <- value
      }
    },
    #' @field uns The uns slot. Must be `NULL` or a named list.
    uns = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_uns, status=done
        read_zarr_element(private$zarr_store, "uns")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_uns, status=done
        value <- private$.validate_named_list(value, "uns")
        write_zarr_element(value, private$zarr_store, "/uns")
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
    #' @param obsm The obsm slot is used to store multi-dimensional annotation
    #'   arrays. It must be either `NULL` or a named list, where each element is a
    #'   matrix with `n_obs` rows and an arbitrary number of columns.
    #' @param varm The varm slot is used to store multi-dimensional annotation
    #'   arrays. It must be either `NULL` or a named list, where each element is a
    #'   matrix with `n_vars` rows and an arbitrary number of columns.
    #' @param obsp The obsp slot is used to store sparse multi-dimensional
    #'   annotation arrays. It must be either `NULL` or a named list, where each
    #'   element is a sparse matrix where each dimension has length `n_obs`.
    #' @param varp The varp slot is used to store sparse multi-dimensional
    #'   annotation arrays. It must be either `NULL` or a named list, where each
    #'   element is a sparse matrix where each dimension has length `n_vars`.
    #' @param uns The uns slot is used to store unstructured annotation. It must
    #'   be either `NULL` or a named list.
    #' @param compression The compression algorithm to use when writing the
    #'  HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
    #' `"none"`.
    #'
    #' @details
    #' The constructor creates a new HDF5 AnnData interface object. This can
    #' either be used to either connect to an existing `.h5ad` file or to
    #' create a new one. To create a new file both `obs_names` and `var_names`
    #' must be specified. In both cases, any additional slots provided will be
    #' set on the created object. This will cause data to be overwritten if the
    #' file already exists.
    initialize = function(store,
                          obs_names = NULL,
                          var_names = NULL,
                          X = NULL,
                          obs = NULL,
                          var = NULL,
                          layers = NULL,
                          obsm = NULL,
                          varm = NULL,
                          obsp = NULL,
                          varp = NULL,
                          uns = NULL,
                          compression = c("none", "gzip", "lzf")) {
      if (!requireNamespace("pizzarr", quietly = TRUE)) {
        stop("The Zarr interface requires the 'pizzarr' package to be installed")
      }

      compression <- match.arg(compression)
      private$.compression <- compression
      
      root <- pizzarr::zarr_open_group(store, path = "/")


      if(length(root$get_attrs()$to_list()) == 0) {
        # Check obs_names and var_names have been provided
        if (is.null(obs_names)) {
          stop("When creating a new .h5ad file, `obs_names` must be defined.")
        }
        if (is.null(var_names)) {
          stop("When creating a new .h5ad file, `var_names` must be defined.")
        }

        # Create an empty H5AD using the provided obs/var names
        write_empty_zarr(store, obs_names, var_names, compression)

        # Set private object slots
        private$zarr_store <- store
        private$.n_obs <- length(obs_names)
        private$.n_vars <- length(var_names)
        private$.obs_names <- obs_names
        private$.var_names <- var_names
      } else {
        # Check the file is a valid H5AD
        
        attrs <- root$get_attrs()$to_list()

        if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
          stop(
            "H5AD encoding information is missing. ",
            "This file may have been created with Python anndata<0.8.0."
          )
        }

        # Set the file path
        private$zarr_store <- store
        private$zarr_root <- root

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

      if (!is.null(obsm)) {
        self$obsm <- obsm
      }

      if (!is.null(varm)) {
        self$varm <- varm
      }

      if (!is.null(obsp)) {
        self$obsp <- obsp
      }

      if (!is.null(varp)) {
        self$varp <- varp
      }

      if (!is.null(uns)) {
        self$uns <- uns
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
#' @param compression The compression algorithm to use when writing the
#'  HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
#' `"none"`.
#'
#' @return An HDF5AnnData object with the same data as the input AnnData
#'   object.
#'
#' @export
#'
#' @examples
#' ad <- AnnData(
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
to_ZarrAnnData <- function(adata, store, compression = c("none", "gzip", "lzf")) { # nolint
  stopifnot(
    inherits(adata, "AbstractAnnData")
  )
  ZarrAnnData$new(
    store = store,
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    obsm = adata$obsm,
    varm = adata$varm,
    obs_names = adata$obs_names,
    var_names = adata$var_names,
    layers = adata$layers,
    obsp = adata$obsp,
    varp = adata$varp,
    uns = adata$uns,
    compression = compression
  )
}
