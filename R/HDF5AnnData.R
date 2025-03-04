#' @title HDF5AnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#' @noRd
HDF5AnnData <- R6::R6Class(
  "HDF5AnnData", # nolint
  inherit = AbstractAnnData,
  cloneable = FALSE,
  private = list(
    .h5obj = NULL,
    .close_on_finalize = FALSE,
    .compression = NULL
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_X, status=done
        read_h5ad_element(private$.h5obj, "X")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_X, status=done
        value <- private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        write_h5ad_element(value, private$.h5obj, "X", private$.compression)
      }
    },
    #' @field layers The layers slot. Must be NULL or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_layers, status=done
        read_h5ad_element(private$.h5obj, "layers")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_layers, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "layers",
          private$.compression
        )
      }
    },
    #' @field obsm The obsm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `obs`.
    obsm = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsm, status=done
        read_h5ad_element(private$.h5obj, "obsm")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obsm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = rownames(self)
        )
        write_h5ad_element(value, private$.h5obj, "obsm")
      }
    },
    #' @field varm The varm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `var`.
    varm = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varm, status=done
        read_h5ad_element(private$.h5obj, "varm")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_varm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = colnames(self)
        )
        write_h5ad_element(value, private$.h5obj, "varm")
      }
    },
    #' @field obsp The obsp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `obs`.
    obsp = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsp, status=done
        read_h5ad_element(private$.h5obj, "obsp")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obsp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = rownames(self),
          expected_colnames = rownames(self)
        )
        write_h5ad_element(value, private$.h5obj, "obsp")
      }
    },
    #' @field varp The varp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `var`.
    varp = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varp, status=done
        read_h5ad_element(private$.h5obj, "varp")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_varp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = colnames(self),
          expected_colnames = colnames(self)
        )
        write_h5ad_element(value, private$.h5obj, "varp")
      }
    },

    #' @field obs The obs slot
    obs = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs, status=done
        read_h5ad_element(private$.h5obj, "obs")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs, status=done
        value <- private$.validate_obsvar_dataframe(value, "obs")
        write_h5ad_element(
          value,
          private$.h5obj,
          "obs",
          private$.compression
        )
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var, status=done
        read_h5ad_element(private$.h5obj, "var")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var, status=done
        value <- private$.validate_obsvar_dataframe(value, "var")
        write_h5ad_element(
          value,
          private$.h5obj,
          "var"
        )
      }
    },
    #' @field obs_names Names of observations
    obs_names = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs_names, status=done
        rownames(self$obs)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs_names, status=done
        rownames(self$obs) <- value
      }
    },
    #' @field var_names Names of variables
    var_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var_names, status=done
        rownames(self$var)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var_names, status=done
        rownames(self$var) <- value
      }
    },
    #' @field uns The uns slot. Must be `NULL` or a named list.
    uns = function(value) {
      if (!private$.h5obj$is_valid) stop("HDF5 file is closed")
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_uns, status=done
        read_h5ad_element(private$.h5obj, "uns")
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_uns, status=done
        value <- private$.validate_named_list(value, "uns")
        write_h5ad_element(value, private$.h5obj, "uns")
      }
    }
  ),
  public = list(
    #' @description HDF5AnnData constructor
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
    #' @param compression The compression algorithm to use when writing the
    #'  HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
    #' `"none"`.
    #' @param mode The mode to open the HDF5 file.
    #'
    #'   * `a` creates a new file or opens an existing one for read/write.
    #'   * `r` opens an existing file for reading.
    #'   * `r+` opens an existing file for read/write.
    #'   * `w` creates a file, truncating any existing ones.
    #'   * `w-`/`x` are synonyms, creating a file and failing if it already exists.
    #'
    #' @details
    #' The constructor creates a new HDF5 AnnData interface object. This can
    #' either be used to either connect to an existing `.h5ad` file or to
    #' create a new one. To create a new file both `obs_names` and `var_names`
    #' must be specified. In both cases, any additional slots provided will be
    #' set on the created object. This will cause data to be overwritten if the
    #' file already exists.
    initialize = function(
      file,
      X = NULL,
      obs = NULL,
      var = NULL,
      layers = NULL,
      obsm = NULL,
      varm = NULL,
      obsp = NULL,
      varp = NULL,
      uns = NULL,
      shape = NULL,
      mode = c("r", "r+", "a", "w", "w-", "x"),
      compression = c("none", "gzip", "lzf")
    ) {
      check_requires("HDF5AnnData", "hdf5r")

      # check arguments
      compression <- match.arg(compression)
      mode <- match.arg(mode)

      # store compression for later use
      private$.compression <- compression

      if (is.character(file) && !file.exists(file)) {
        # store private values
        private$.h5obj <- hdf5r::H5File$new(file, mode = "w-")
        private$.close_on_finalize <- TRUE

        # Determine initial obs and var
        shape <- get_shape(obs, var, X, shape)
        obs <- get_initial_obs(obs, X, shape)
        var <- get_initial_var(var, X, shape)

        # Create an empty H5AD
        write_empty_h5ad(private$.h5obj, obs, var, compression)

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
        open_hdf5_file <- is.character(file)
        if (open_hdf5_file) {
          file <- hdf5r::H5File$new(file, mode = mode)
        }

        if (!inherits(file, "H5File")) {
          stop("file must be a character string or an H5File object")
        }

        # Check the file is a valid H5AD
        attrs <- hdf5r::h5attributes(file)

        if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
          stop(
            "H5AD encoding information is missing. ",
            "This file may have been created with Python anndata<0.8.0."
          )
        }

        # Set the file path
        private$.h5obj <- file
        private$.close_on_finalize <- open_hdf5_file

        # assert other arguments are NULL
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
      }
    },

    #' @description Close the HDF5 file when the object is garbage collected
    finalize = function() {
      if (private$.close_on_finalize) {
        self$close()
      }
      return(invisible(self))
    },

    #' @description Close the HDF5 file
    close = function() {
      if (private$.h5obj$is_valid) {
        private$.h5obj$close()
      }
    },

    #' @description Number of observations in the AnnData object
    n_obs = function() {
      nrow(self$obs)
    },

    #' @description Number of variables in the AnnData object
    n_vars = function() {
      nrow(self$var)
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
#' @param mode The mode to open the HDF5 file.
#'
#'   * `a` creates a new file or opens an existing one for read/write.
#'   * `r` opens an existing file for reading.
#'   * `r+` opens an existing file for read/write.
#'   * `w` creates a file, truncating any existing ones.
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists.
#'
#' @return An HDF5AnnData object with the same data as the input AnnData
#'   object.
#'
#' @noRd
#'
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5),
#' )
#' to_HDF5AnnData(ad, "test.h5ad")
#' # remove file
#' file.remove("test.h5ad")
# nolint start: object_name_linter
to_HDF5AnnData <- function(
  # nolint end: object_name_linter
  adata,
  file,
  compression = c("none", "gzip", "lzf"),
  mode = c("w-", "r", "r+", "a", "w", "x")
) {
  stopifnot(
    inherits(adata, "AbstractAnnData")
  )
  mode <- match.arg(mode)
  HDF5AnnData$new(
    file = file,
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    obsm = adata$obsm,
    varm = adata$varm,
    layers = adata$layers,
    obsp = adata$obsp,
    varp = adata$varp,
    uns = adata$uns,
    compression = compression,
    shape = adata$shape(),
    mode = mode
  )
}

# nolint start: object_name_linter
cleanup_HDF5AnnData <- function(...) {
  # nolint end: object_name_linter
  args <- list(...)

  if (
    !is.null(args$file) && is.character(args$file) && file.exists(args$file)
  ) {
    cli::cli_alert("Removing file: ", args$file)
    unlink(args$file)
  }
}
