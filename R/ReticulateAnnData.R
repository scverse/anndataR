#' ReticulateAnnData class
#'
#' An R6 class that provides an interface to interact with AnnData objects in Python
#' using the reticulate package.
ReticulateAnnData <- R6::R6Class("ReticulateAnnData", # nolint
  inherit = AbstractAnnData,
  private = list(
    .obj = NULL
  ),
  active = list(
    #' @field X NULL or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and
    #'   `var`.
    X = function(value) {
      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_X, status=done
        py_to_r_ifneedbe(private$.obj$X)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_X, status=done
        private$.obj$X <- private$.validate_matrix(value, "X")
        self
      }
    },
    #' @field layers NULL or a named list with all elements having the
    #'   dimensions consistent with `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_layers, status=done
        py_to_r_ifneedbe(private$.obj$layers)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_layers, status=done
        private$.obj$layers <- private$.validate_layers(value)
        self
      }
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_obs, status=done
        py_to_r_ifneedbe(private$.obj$obs)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obs, status=done
        obs_ <- private$.validate_obsvar_dataframe(value, "obs")
        reticulate::py_set_attr(private$.obj, "obs", obs_)
        self
      }
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_var, status=done
        py_to_r_ifneedbe(private$.obj$var)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_var, status=done
        var_ <- private$.validate_obsvar_dataframe(value, "var")
        reticulate::py_set_attr(private$.obj, "var", var_)
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
        # trackstatus: class=ReticulateAnnData, feature=get_obs_names, status=done
        py_to_r_ifneedbe(private$.obj$obs_names)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obs_names, status=done
        private$.obj$obs_names <- private$.validate_obsvar_names(value, "obs")
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
        # trackstatus: class=ReticulateAnnData, feature=get_var_names, status=done
        py_to_r_ifneedbe(private$.obj$var_names)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_var_names, status=done
        private$.obj$var_names <- private$.validate_obsvar_names(value, "var")
        self
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #'   Inherits from ReticulateAnnData.
    #'
    #' @param file The filename (character) of the `.h5ad` file. If this
    #'   file does not exist yet, `obs_names` and `var_names` must be provided.
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
    #' @param shape Shape tuple (#observations, #variables). Can be provided
    #'   if `X` or `obs` and `var` are not provided.
    initialize = function(file = NULL,
                          X = NULL,
                          obs = NULL,
                          var = NULL,
                          layers = NULL,
                          obsm = NULL,
                          varm = NULL,
                          obsp = NULL,
                          varp = NULL,
                          uns = NULL,
                          shape = NULL) {
      python_anndata <- .reticulate_load_anndata()

      # in case the user also loaded the 'anndata' package
      if (inherits(file, "AnnDataR6")) {
        file <- file$`.__enclos_env__`$private$.anndata
      }

      # file should be a Python AnnData object or a path
      if (
        !is.null(file) &&
          !inherits(file, "anndata._core.anndata.AnnData") &&
          !is.character(file)) {
        stop(
          "Argument 'file' should be NULL, a character path, or a Python AnnData created by ",
          "`anndata <- reticulate::anndata(\"anndata\"); file <- anndata$read_h5ad(...)`."
        )
      }

      if (is.null(file)) {
        # create a new h5ad from scratch

        # Determine initial obs and var
        shape <- get_shape(obs, var, X, shape)
        obs <- get_initial_obs(obs, X, shape)
        var <- get_initial_var(var, X, shape)

        # Create an empty H5AD
        private$.obj <- python_anndata$AnnData(
          obs = obs,
          var = var
        )

        # set other slots
        if (!is.null(X)) {
          self$X <- X
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
      } else {
        if (!is.null(obs)) {
          stop("obs must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(var)) {
          stop("var must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(X)) {
          stop("X must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(layers)) {
          stop("layers must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(obsm)) {
          stop("obsm must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(varm)) {
          stop("varm must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(obsp)) {
          stop("obsp must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(varp)) {
          stop("varp must be NULL when loading an existing .h5ad file")
        }
        if (!is.null(uns)) {
          stop("uns must be NULL when loading an existing .h5ad file")
        }

        if (is.character(file)) {
          if (!file.exists(file)) {
            stop("File '", file, "' does not exist.")
          }

          file <- python_anndata$read_h5ad(file)
        }

        # store file
        private$.obj <- file
      }
    },

    # copied from https://github.com/dynverse/anndata/blob/main/R/class_anndata.R#LL562C1-L601C7
    #' @description Write .h5ad-formatted hdf5 file.
    #'
    #' Generally, if you have sparse data that are stored as a dense matrix, you can
    #' dramatically improve performance and reduce disk space by converting to a csr_matrix:
    #'
    #' @param filename Filename of data file. Defaults to backing file.
    #' @param compression See the h5py
    #'   [filter pipeline](http://docs.h5py.org/en/latest/high/dataset.html#dataset-compression).
    #'   Options are `"gzip"`, `"lzf"` or `NULL`.
    #' @param compression_opts See the h5py
    #'   [filter pipeline](http://docs.h5py.org/en/latest/high/dataset.html#dataset-compression).
    #' @param as_dense Sparse in AnnData object to write as dense. Currently only supports `"X"` and `"raw/X"`.
    write_h5ad = function(filename, compression = NULL, compression_opts = NULL, as_dense = list()) {
      filename <- normalizePath(filename, mustWork = FALSE)
      invisible(py_to_r_ifneedbe(private$.obj$write_h5ad(
        filename = filename,
        compression = compression,
        compression_opts = compression_opts,
        as_dense = as_dense
      )))
    }
  )
)


#' Convert an AnnData object to a ReticulateAnnData object
#'
#' This function takes an AnnData object and converts it to a ReticulateAnnData
#' object, loading all fields into memory.
#'
#' @param adata An AnnData object to be converted to ReticulateAnnData.
#'
#' @return A ReticulateAnnData object with the same data as the input AnnData
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
#' to_ReticulateAnnData(ad)
to_ReticulateAnnData <- function(adata) { # nolint
  ReticulateAnnData$new(
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    obsm = adata$obsm,
    varm = adata$varm,
    layers = adata$layers,
    obsp = adata$obsp,
    varp = adata$varp,
    uns = adata$uns,
    shape = adata$shape()
  )
}

.reticulate_load_anndata <- function() {
  # check whether reticulate is installed
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The Reticulate interface requires the 'reticulate' package to be installed")
  }

  # check whether python anndata is installed
  tryCatch(
    {
      reticulate::import("anndata", convert = FALSE)
    },
    error = function(e) {
      stop(
        "Could not find the Python anndata package.\n  Please run ",
        "`reticulate::py_install(\"anndata\")` to solve this issue."
      )
    }
  )
}
