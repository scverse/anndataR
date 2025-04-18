#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
#'
#' @importFrom Matrix as.matrix
#'
#' @noRd
#'
#' @examples
#' ## complete example
#' ad <- AnnData(
#'   X = matrix(1:15, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' ad
#'
#' ## minimum example
#' AnnData(
#'   obs = data.frame(row.names = letters[1:10]),
#'   var = data.frame(row.names = LETTERS[1:5])
#' )
InMemoryAnnData <- R6::R6Class(
  "InMemoryAnnData", # nolint
  inherit = AbstractAnnData,
  private = list(
    .X = NULL,
    .layers = NULL,
    .obs = NULL,
    .var = NULL,
    .obsm = NULL,
    .varm = NULL,
    .obsp = NULL,
    .varp = NULL,
    .uns = NULL
  ),
  active = list(
    #' @field X NULL or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and
    #'   `var`.
    X = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_X, status=done
        private$.X
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_X, status=done
        private$.X <- private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        self
      }
    },
    #' @field layers NULL or a named list with all elements having the
    #'   dimensions consistent with `obs` and `var`.
    layers = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_layers, status=done
        private$.layers
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_layers, status=done
        private$.layers <- private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        self
      }
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_obs, status=done
        private$.obs
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_obs, status=done
        private$.obs <- private$.validate_obsvar_dataframe(value, "obs")
        self
      }
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_var, status=done
        private$.var
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_var, status=done
        private$.var <- private$.validate_obsvar_dataframe(value, "var")
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
        # trackstatus: class=InMemoryAnnData, feature=get_obs_names, status=done
        rownames(private$.obs)
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_obs_names, status=done
        rownames(private$.obs) <- private$.validate_obsvar_names(value, "obs")
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
        # trackstatus: class=InMemoryAnnData, feature=get_var_names, status=done
        rownames(private$.var)
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_var_names, status=done
        rownames(private$.var) <- private$.validate_obsvar_names(value, "var")
        self
      }
    },
    #' @field obsm The obsm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `obs`.
    obsm = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_obsm, status=done
        private$.obsm
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_obsm, status=done
        private$.obsm <- private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = rownames(self)
        )
        self
      }
    },
    #' @field varm The varm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `var`.
    varm = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_varm, status=done
        private$.varm
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_varm, status=done
        private$.varm <- private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = colnames(self)
        )
        self
      }
    },
    #' @field obsp The obsp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `obs`.
    obsp = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_obsp, status=done
        private$.obsp
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_obsp, status=done
        private$.obsp <- private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = rownames(self),
          expected_colnames = rownames(self)
        )
        self
      }
    },
    #' @field varp The varp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `var`.
    varp = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_varp, status=done
        private$.varp
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_varp, status=done
        private$.varp <- private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = colnames(self),
          expected_colnames = colnames(self)
        )
        self
      }
    },
    #' @field uns The uns slot. Must be `NULL` or a named list.
    uns = function(value) {
      if (missing(value)) {
        # trackstatus: class=InMemoryAnnData, feature=get_uns, status=done
        private$.uns
      } else {
        # trackstatus: class=InMemoryAnnData, feature=set_uns, status=done
        private$.uns <- private$.validate_named_list(value, "uns")
        self
      }
    }
  ),
  public = list(
    #' @description Creates a new instance of an in memory AnnData object.
    #'   Inherits from AbstractAnnData.
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
    #' @param uns The uns slot is used to store unstructured annotation.
    #'   It must be either `NULL` or a named list.
    #' @param shape Shape tuple (#observations, #variables). Can be provided
    #'   if `X` or `obs` and `var` are not provided.
    initialize = function(
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
      # Determine initial obs and var
      shape <- get_shape(obs, var, X, shape)
      obs <- get_initial_obs(obs, X, shape)
      var <- get_initial_var(var, X, shape)

      # set obs and var first
      if (!is.data.frame(obs)) {
        cli_abort("{.arg obs} must be a {.cls data.frame}")
      }
      if (!is.data.frame(var)) {
        cli_abort("{.arg var} must be a {.cls data.frame}")
      }
      private$.obs <- obs
      private$.var <- var

      # write other slots later
      self$X <- X
      self$layers <- layers
      self$obsm <- obsm
      self$varm <- varm
      self$obsp <- obsp
      self$varp <- varp
      self$uns <- uns
    }
  )
)

#' Convert an AnnData object to an InMemoryAnnData object
#'
#' This function takes an AnnData object and converts it to an InMemoryAnnData
#' object, loading all fields into memory.
#'
#' @param adata An AnnData object to be converted to InMemoryAnnData.
#'
#' @return An InMemoryAnnData object with the same data as the input AnnData
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
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' to_InMemoryAnnData(ad)
# nolint start object_name_linter
to_InMemoryAnnData <- function(adata) {
  # nolint end object_name_linter
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  InMemoryAnnData$new(
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
