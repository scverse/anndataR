#' @title ReticulateAnnData
#'
#' @description
#' Implementation of an AnnData object that wraps a Python anndata.AnnData object
#' using reticulate. This allows direct interaction with Python AnnData objects
#' while maintaining the R interface.
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @return A `ReticulateAnnData` object
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
ReticulateAnnData <- R6::R6Class(
  "ReticulateAnnData", # nolint
  inherit = AbstractAnnData,
  cloneable = FALSE,
  private = list(
    .py_anndata = NULL,

    .check_py_object_valid = function() {
      if (is.null(private$.py_anndata)) {
        cli_abort("Python AnnData object is not initialized")
      }
    }
  ),
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_X, status=done
        py_to_r(private$.py_anndata$X)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_X, status=done
        value <- private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        private$.py_anndata$X <- r_to_py(value)
        self
      }
    },

    #' @field layers See [AnnData-usage]
    layers = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_layers, status=done
        bi <- reticulate::import_builtins()
        out <- list()
        keys <- bi$list(private$.py_anndata$layers$keys())
        for (name in keys) {
          out[[name]] <- py_to_r(private$.py_anndata$layers[[name]])
        }
        out
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_layers, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        private$.py_anndata$layers <- r_to_py(value)
        self
      }
    },

    #' @field obs See [AnnData-usage]
    obs = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_obs, status=done
        py_to_r(private$.py_anndata$obs)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obs, status=done
        value <- private$.validate_obsvar_dataframe(value, "obs")
        private$.py_anndata$obs <- r_to_py(value)
        self
      }
    },

    #' @field var See [AnnData-usage]
    var = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_var, status=done
        py_to_r(private$.py_anndata$var)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_var, status=done
        value <- private$.validate_obsvar_dataframe(value, "var")
        private$.py_anndata$var <- r_to_py(value)
        self
      }
    },

    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_obs_names, status=done
        bi <- reticulate::import_builtins()
        reticulate::py_to_r(bi$list(private$.py_anndata$obs_names))
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obs_names, status=done
        private$.py_anndata$obs_names <- reticulate::r_to_py(value)
        self
      }
    },

    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_var_names, status=done
        bi <- reticulate::import_builtins()
        reticulate::py_to_r(bi$list(private$.py_anndata$var_names))
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_var_names, status=done
        private$.py_anndata$var_names <- reticulate::r_to_py(value)
        self
      }
    },

    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_obsm, status=done
        py_to_r(private$.py_anndata$obsm)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obsm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = rownames(self)
        )
        private$.py_anndata$obsm <- r_to_py(value)
        self
      }
    },

    #' @field varm See [AnnData-usage]
    varm = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_varm, status=done
        py_to_r(private$.py_anndata$varm)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_varm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = colnames(self)
        )
        private$.py_anndata$varm <- r_to_py(value)
        self
      }
    },

    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_obsp, status=done
        py_to_r(private$.py_anndata$obsp)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_obsp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = rownames(self),
          expected_colnames = rownames(self)
        )
        private$.py_anndata$obsp <- r_to_py(value)
        self
      }
    },

    #' @field varp See [AnnData-usage]
    varp = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_varp, status=done
        py_to_r(private$.py_anndata$varp)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_varp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = colnames(self),
          expected_colnames = colnames(self)
        )
        private$.py_anndata$varp <- r_to_py(value)
        self
      }
    },

    #' @field uns See [AnnData-usage]
    uns = function(value) {
      private$.check_py_object_valid()

      if (missing(value)) {
        # trackstatus: class=ReticulateAnnData, feature=get_uns, status=done
        py_to_r(private$.py_anndata$uns)
      } else {
        # trackstatus: class=ReticulateAnnData, feature=set_uns, status=done
        value <- private$.validate_named_list(value, "uns")
        private$.py_anndata$uns <- r_to_py(value)
        self
      }
    }
  ),
  public = list(
    #' @description
    #' `ReticulateAnnData` constructor
    #'
    #' @param py_anndata A Python AnnData object created using reticulate, or
    #'   NULL to create a new empty Python AnnData object
    #' @param X See the `X` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param layers See the `layers` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param obs See the `obs` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param var See the `var` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param obsm See the `obsm` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param varm See the `varm` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param obsp See the `obsp` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param varp See the `varp` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param uns See the `uns` slot in [AnnData-usage] (only used if py_anndata is NULL)
    #' @param shape Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if
    #'   both `X` or `obs` and `var` are not provided. (only used if py_anndata is NULL)
    #'
    #' @details
    #' The constructor creates a new ReticulateAnnData interface object that
    #' wraps a Python AnnData object. If `py_anndata` is provided, it must be
    #' a valid Python AnnData object. If NULL, a new Python AnnData object
    #' will be created using the other provided arguments.
    initialize = function(
      py_anndata = NULL,
      X = NULL,
      obs = NULL,
      var = NULL,
      layers = NULL,
      obsm = NULL,
      varm = NULL,
      obsp = NULL,
      varp = NULL,
      uns = NULL,
      shape = NULL
    ) {
      check_requires("ReticulateAnnData", "reticulate", where = "CRAN")
      check_requires("ReticulateAnnData", "anndata", where = "Python")

      if (!is.null(py_anndata)) {
        # Use existing Python AnnData object
        if (!inherits(py_anndata, "anndata._core.anndata.AnnData")) {
          cli_abort(
            "{.arg py_anndata} must be a Python AnnData object from {.pkg reticulate}"
          )
        }

        private$.py_anndata <- py_anndata
      } else {
        # Create new Python AnnData object
        anndata <- reticulate::import("anndata", convert = FALSE)

        if (!is.null(X)) {
          # Convert X to Python and create AnnData with proper obs/var names
          py_x <- r_to_py(X)

          # Create obs and var DataFrames with proper row names
          n_obs <- nrow(X)
          n_vars <- ncol(X)

          # Get row and column names, or create defaults
          obs_names <- if (!is.null(rownames(X))) {
            rownames(X)
          } else {
            paste0("cell", seq_len(n_obs))
          }
          var_names <- if (!is.null(colnames(X))) {
            colnames(X)
          } else {
            paste0("gene", seq_len(n_vars))
          }

          # Create DataFrames
          pandas <- reticulate::import("pandas", convert = FALSE)
          py_obs <- pandas$DataFrame(index = reticulate::r_to_py(obs_names))
          py_var <- pandas$DataFrame(index = reticulate::r_to_py(var_names))

          private$.py_anndata <- anndata$AnnData(
            X = py_x,
            obs = py_obs,
            var = py_var
          )
        } else {
          shape <- get_shape(obs, var, X, shape)
          obs <- get_initial_obs(obs, X, shape)
          var <- get_initial_var(var, X, shape)

          # Create empty AnnData with specified obs and var
          pandas <- reticulate::import("pandas", convert = FALSE)

          # For Python AnnData, we need proper index names
          py_obs <- r_to_py(obs)
          py_var <- r_to_py(var)

          # Set proper index names if they don't exist or are auto-generated
          if (
            is.null(rownames(obs)) ||
              all(rownames(obs) == "") ||
              all(rownames(obs) == as.character(seq_len(nrow(obs))))
          ) {
            py_obs$index <- reticulate::r_to_py(paste0(
              "obs",
              seq_len(nrow(obs))
            ))
          }
          if (
            is.null(rownames(var)) ||
              all(rownames(var) == "") ||
              all(rownames(var) == as.character(seq_len(nrow(var))))
          ) {
            py_var$index <- reticulate::r_to_py(paste0(
              "var",
              seq_len(nrow(var))
            ))
          }

          private$.py_anndata <- anndata$AnnData(
            obs = py_obs,
            var = py_var
          )
        }

        # Set other slots if provided (but not obs/var since they're already set)
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
        if (!is.null(uns)) self$uns <- uns
      }

      self
    },

    #' @description See the `n_obs` field in [AnnData-usage]
    n_obs = function() {
      private$.check_py_object_valid()
      as.integer(reticulate::py_to_r(private$.py_anndata$n_obs))
    },

    #' @description See the `n_vars` field in [AnnData-usage]
    n_vars = function() {
      private$.check_py_object_valid()
      as.integer(reticulate::py_to_r(private$.py_anndata$n_vars))
    },

    #' @description Get the underlying Python AnnData object
    #' @return The Python AnnData object wrapped by this ReticulateAnnData
    py_anndata = function() {
      private$.check_py_object_valid()
      private$.py_anndata
    }
  )
)

#' Convert an `AnnData` to a `ReticulateAnnData`
#'
#' Convert another `AnnData` object to a [`ReticulateAnnData`] object
#'
#' @param adata An `AnnData` object to be converted to [`ReticulateAnnData`]
#'
#' @return A [`ReticulateAnnData`] object with the same data as the input
#'   `AnnData` object
#' @keywords internal
#'
#' @family object converters
#'
#' @examples
#' \dontrun{
#' # Requires Python anndata to be installed
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' ad$as_ReticulateAnnData()
#' }
# nolint start: object_name_linter
as_ReticulateAnnData <- function(adata) {
  # nolint end: object_name_linter
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  ReticulateAnnData$new(
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    layers = adata$layers,
    obsm = adata$obsm,
    varm = adata$varm,
    obsp = adata$obsp,
    varp = adata$varp,
    uns = adata$uns,
    shape = adata$shape()
  )
}
