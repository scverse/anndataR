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
        obs_ <- py_to_r_ifneedbe(private$.obj$obs)
        rownames(obs_) <- NULL
        attr(obs_, "pandas.index") <- NULL # nolint
        obs_
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obs, status=done
        # orig_obs_names <- self$obs_names
        # private$.obj$obs <- private$.validate_obsvar_dataframe(value, "obs")
        # self$obs_names <- orig_obs_names
        obs_ <- private$.validate_obsvar_dataframe(value, "obs")
        reticulate::py_del_attr(private$.obj, "obs")
        for (name in colnames(obs_)) {
          reticulate::py_set_item(private$.obj$obs, name, obs_[[name]])
        }
        self
      }
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_var, status=done
        var_ <- py_to_r_ifneedbe(private$.obj$var)
        rownames(var_) <- NULL
        attr(var_, "pandas.index") <- NULL # nolint
        var_
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_var, status=done
        # orig_var_names <- self$var_names
        # private$.obj$var <- private$.validate_obsvar_dataframe(value, "var")
        # self$var_names <- orig_var_names
        var_ <- private$.validate_obsvar_dataframe(value, "var")
        reticulate::py_del_attr(private$.obj, "var")
        for (name in colnames(var_)) {
          reticulate::py_set_item(private$.obj$var, name, var_[[name]])
        }
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
    #' @param file The filename (character) of the `.h5ad` file. Alternatively,
    #' you can also provide an object created by
    #' `anndata <- reticulate::import("anndata"); anndata$read_h5ad(...)`. If this
    #' file does not exist yet, you must provide values for `X`, `obs`, `var` to
    #' create a new AnnData with.
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
    initialize = function(obs_names, var_names, file = NULL, X = NULL, obs = NULL, var = NULL, layers = NULL) {
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
        if (missing(X)) X <- NULL
        if (missing(obs)) obs <- NULL
        if (missing(var)) var <- NULL
        if (missing(obs_names)) stop("When creating a new .h5ad file, `obs_names` must be defined.")
        if (missing(var_names)) stop("When creating a new .h5ad file, `var_names` must be defined.")
        if (missing(layers)) layers <- NULL

        private$.obj <- python_anndata$AnnData(shape = c(length(obs_names), length(var_names)))

        # write obs and var first, because these are used by other validators
        self$obs_names <- obs_names
        self$var_names <- var_names

        # write other slots later
        self$obs <- obs
        self$var <- var
        self$X <- X
        self$layers <- layers
      } else {
        # check if other arguments are defined
        slots_missing <- missing(X) && missing(obs) && missing(var) &&
          missing(obs_names) && missing(var_names) && missing(layers)
        if (!slots_missing) {
          stop(
            "Failed to create ReticulateAnnData object. ",
            "Reason: If the file provided already exists, ",
            "then the other arguments (X, obs, var, ...) ",
            "should not be passed any value."
          )
        }

        # open h5 file if this isn't the already
        if (is.character(file)) {
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
    obs_names = adata$obs_names,
    var_names = adata$var_names,
    layers = adata$layers
  )
}

.reticulate_load_anndata <- function() {
  # check whether reticulate is installed
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The Reticulate interface requires the 'reticulate' package to be installed")
  }

  # check whether python anndata is installed
  tryCatch({
    reticulate::import("anndata", convert = FALSE)
  }, error = function(e) {
    stop("Could not find the Python anndata package.\n  Please run ",
    "`reticulate::py_install(\"anndata\")` to solve this issue.")
  })
}