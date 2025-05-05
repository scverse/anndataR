#' @title InMemoryAnnData
#'
#' @description
#' Implementation of an in-memory AnnData object.
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
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
#'
#' @importFrom Matrix as.matrix
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
    #' @field X See [AnnData-usage]
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
    #' @field layers See [AnnData-usage]
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
    #' @field obs See [AnnData-usage]
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
    #' @field var See [AnnData-usage]
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
    #' @field obs_names See [AnnData-usage]
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
    #' @field var_names See [AnnData-usage]
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
    #' @field obsm See [AnnData-usage]
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
    #' @field varm See [AnnData-usage]
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
    #' @field obsp See [AnnData-usage]
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
    #' @field varp See [AnnData-usage]
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
    #' @field uns See [AnnData-usage]
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
    #' @description Creates a new instance of an in-memory `AnnData` object.
    #'   Inherits from [AbstractAnnData].
    #'
    #' @param X See the `X` slot in [AnnData-usage]
    #' @param layers See the `layers` slot in [AnnData-usage]
    #' @param obs See the `obs` slot in [AnnData-usage]
    #' @param var See the `var` slot in [AnnData-usage]
    #' @param obsm See the `obsm` slot in [AnnData-usage]
    #' @param varm See the `varm` slot in [AnnData-usage]
    #' @param obsp See the `obsp` slot in [AnnData-usage]
    #' @param varp See the `varp` slot in [AnnData-usage]
    #' @param uns See the `uns` slot in [AnnData-usage]
    #' @param shape Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if
    #'   both `X` or `obs` and `var` are not provided.
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

# See as_InMemoryAnnData() for function documentation
to_InMemoryAnnData <- function(adata) {
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
