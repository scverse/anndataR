#' @title A View of an AnnData Object
#'
#' @description A lazy view of an AnnData object that allows applying subsetting operations
#' without immediately executing them. Subsetting is applied when converting to
#' a concrete AnnData implementation (InMemoryAnnData, HDF5AnnData) or other
#' formats (SingleCellExperiment, Seurat).
#'
#' @return A `AnnDataView` object
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
#'
#' @examples
#' # Create a base AnnData object
#' ad <- AnnData(
#'   X = matrix(1:15, 3L, 5L),
#'   obs = data.frame(row.names = LETTERS[1:3], cell_type = c("A", "B", "A")),
#'   var = data.frame(row.names = letters[1:5], gene_type = c("X", "Y", "X", "Y", "Z"))
#' )
#'
#' # Create a view with lazy subsetting using S3 [ method
#' view <- ad[ad$obs$cell_type == "A", ad$var$gene_type %in% c("X", "Y")]
#'
#' # Apply subsetting by converting to a concrete implementation
#' result <- view$as_InMemoryAnnData()
#'
#' @export
AnnDataView <- R6::R6Class(
  "AnnDataView",
  inherit = AbstractAnnData,
  private = list(
    .base_adata = NULL,
    .obs_subset = NULL,
    .var_subset = NULL,

    # Apply subsetting to matrices/data frames with explicit dimension control
    .apply_subset = function(
      x,
      obs_on_rows = FALSE,
      var_on_rows = FALSE,
      obs_on_cols = FALSE,
      var_on_cols = FALSE
    ) {
      if (is.null(x)) {
        return(x)
      }

      # Check if x is a subsettable matrix-like object (data.frame, matrix, or Matrix)
      if (is.data.frame(x) || is.matrix(x) || inherits(x, "Matrix")) {
        result <- x

        # Apply obs subsetting to rows if requested
        if (obs_on_rows && !is.null(private$.obs_subset)) {
          result <- result[private$.obs_subset, , drop = FALSE]
        }

        # Apply var subsetting to rows if requested
        if (var_on_rows && !is.null(private$.var_subset)) {
          result <- result[private$.var_subset, , drop = FALSE]
        }

        # Apply obs subsetting to columns if requested
        if (obs_on_cols && !is.null(private$.obs_subset)) {
          result <- result[, private$.obs_subset, drop = FALSE]
        }

        # Apply var subsetting to columns if requested
        if (var_on_cols && !is.null(private$.var_subset)) {
          result <- result[, private$.var_subset, drop = FALSE]
        }

        return(result)
      }

      x
    },

    # Apply subsetting to vectors
    .apply_vector_subset = function(x, use_obs = FALSE, use_var = FALSE) {
      if (is.null(x) || !is.vector(x)) {
        return(x)
      }

      if (use_obs && !is.null(private$.obs_subset)) {
        return(x[private$.obs_subset])
      }

      if (use_var && !is.null(private$.var_subset)) {
        return(x[private$.var_subset])
      }

      x
    },

    # Helper to show error for setter operations
    .abort_setter = function(field_name) {
      cli_abort(
        paste0(
          "Cannot set {.field {field_name}} on an {.cls AnnDataView} object. Convert to a concrete implementation ",
          "first using any of the {.code as_*()} methods ",
          "(e.g., {.code adata$as_InMemoryAnnData()}, {.code adata$as_HDF5AnnData()})."
        )
      )
    },

    # Override class name to show "View of BaseClass"
    .class_name = function() {
      base_class <- class(private$.base_adata)[1]
      paste("View of", base_class)
    }
  ),
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      if (!missing(value)) {
        private$.abort_setter("X")
      }

      private$.apply_subset(
        private$.base_adata$X,
        obs_on_rows = TRUE,
        var_on_cols = TRUE
      )
    },

    #' @field layers See [AnnData-usage]
    layers = function(value) {
      if (!missing(value)) {
        private$.abort_setter("layers")
      }

      layers <- private$.base_adata$layers
      if (is.null(layers)) {
        return(NULL)
      }

      # Apply subsetting to each layer (obs on rows, var on columns)
      lapply(layers, function(x) {
        private$.apply_subset(x, obs_on_rows = TRUE, var_on_cols = TRUE)
      })
    },

    #' @field obs See [AnnData-usage]
    obs = function(value) {
      if (!missing(value)) {
        private$.abort_setter("obs")
      }

      private$.apply_subset(private$.base_adata$obs, obs_on_rows = TRUE)
    },

    #' @field var See [AnnData-usage]
    var = function(value) {
      if (!missing(value)) {
        private$.abort_setter("var")
      }

      private$.apply_subset(private$.base_adata$var, var_on_rows = TRUE)
    },

    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      if (!missing(value)) {
        private$.abort_setter("obs_names")
      }

      private$.apply_vector_subset(
        private$.base_adata$obs_names,
        use_obs = TRUE
      )
    },

    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      if (!missing(value)) {
        private$.abort_setter("var_names")
      }

      private$.apply_vector_subset(
        private$.base_adata$var_names,
        use_var = TRUE
      )
    },

    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      if (!missing(value)) {
        private$.abort_setter("obsm")
      }

      obsm <- private$.base_adata$obsm
      if (is.null(obsm)) {
        return(NULL)
      }

      # Apply obs subsetting to each matrix (filter rows only)
      lapply(obsm, function(x) {
        private$.apply_subset(x, obs_on_rows = TRUE)
      })
    },

    #' @field varm See [AnnData-usage]
    varm = function(value) {
      if (!missing(value)) {
        private$.abort_setter("varm")
      }

      varm <- private$.base_adata$varm
      if (is.null(varm)) {
        return(NULL)
      }

      # Apply var subsetting to each matrix (filter rows only)
      lapply(varm, function(x) {
        private$.apply_subset(x, var_on_rows = TRUE)
      })
    },

    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      if (!missing(value)) {
        private$.abort_setter("obsp")
      }

      obsp <- private$.base_adata$obsp
      if (is.null(obsp)) {
        return(NULL)
      }

      # Apply obs subsetting to each matrix (both rows and columns)
      lapply(obsp, function(x) {
        private$.apply_subset(x, obs_on_rows = TRUE, obs_on_cols = TRUE)
      })
    },

    #' @field varp See [AnnData-usage]
    varp = function(value) {
      if (!missing(value)) {
        private$.abort_setter("varp")
      }

      varp <- private$.base_adata$varp
      if (is.null(varp)) {
        return(NULL)
      }

      # Apply var subsetting to each matrix (both rows and columns)
      lapply(varp, function(x) {
        private$.apply_subset(x, var_on_rows = TRUE, var_on_cols = TRUE)
      })
    },

    #' @field uns See [AnnData-usage]
    uns = function(value) {
      if (!missing(value)) {
        private$.abort_setter("uns")
      }

      # uns is not affected by subsetting
      private$.base_adata$uns
    }
  ),
  public = list(
    #' @description
    #' Create a new AnnDataView object
    #'
    #' @param base_adata An existing AnnData object to create a view of
    #' @param i Optional initial obs subset (logical, integer, or character vector)
    #' @param j Optional initial var subset (logical, integer, or character vector)
    #'
    #' @return A new `AnnDataView` object
    initialize = function(base_adata, i, j) {
      if (!inherits(base_adata, "AbstractAnnData")) {
        cli_abort("{.arg base_adata} must be an {.cls AbstractAnnData} object")
      }

      private$.base_adata <- base_adata

      private$.obs_subset <- convert_to_indices(
        i,
        base_adata$n_obs(),
        base_adata$obs_names,
        "observations"
      )
      private$.var_subset <- convert_to_indices(
        j,
        base_adata$n_vars(),
        base_adata$var_names,
        "variables"
      )

      invisible(self)
    },
    #' Subset the AnnDataView
    #'
    #' @param i Row indices (observations). Can be numeric, logical, or character.
    #' @param j Column indices (variables). Can be numeric, logical, or character.
    #'
    #' @return A new `AnnDataView` object
    subset = function(i, j) {
      if (missing(i)) {
        i <- NULL
      }
      if (missing(j)) {
        j <- NULL
      }

      # processing observation indices
      processed_i <- convert_to_indices(
        i,
        self$n_obs(),
        self$obs_names,
        "observations"
      )
      new_i <- if (is.null(processed_i)) {
        private$.obs_subset
      } else if (is.null(private$.obs_subset)) {
        processed_i
      } else {
        private$.obs_subset[processed_i]
      }

      # processing variable indices
      processed_j <- convert_to_indices(
        j,
        self$n_vars(),
        self$var_names,
        "variables"
      )
      new_j <- if (is.null(processed_j)) {
        private$.var_subset
      } else if (is.null(private$.var_subset)) {
        processed_j
      } else {
        private$.var_subset[processed_j]
      }

      AnnDataView$new(private$.base_adata, new_i, new_j)
    }
  )
)

#' Convert various subset types to integer indices
#'
#' @param subset The subset condition (logical, integer, or character vector)
#' @param max_length The maximum allowed length/index
#' @param names_vector The names to match against for character subsetting
#' @param context_name Name for error messages ("observations" or "variables")
#'
#' @return Integer vector of indices, or NULL if subset is NULL
#' @keywords internal
#' @noRd
convert_to_indices <- function(
  subset,
  max_length,
  names_vector,
  context_name = "items"
) {
  if (missing(subset) || is.null(subset)) {
    return(NULL)
  }

  if (is.logical(subset)) {
    if (length(subset) != max_length) {
      cli_abort(
        "Logical subset of {context_name} must have length {max_length}"
      )
    }

    which(subset)
  } else if (is.integer(subset) || is.numeric(subset)) {
    subset <- as.integer(subset)
    if (any(subset < 1 | subset > max_length)) {
      cli_abort(
        "Integer indices for {context_name} must be between 1 and {max_length}"
      )
    }

    subset
  } else if (is.character(subset)) {
    if (is.null(names_vector)) {
      cli_abort(c(
        "Cannot use character names to subset {context_name}",
        "i" = "The {context_name} do not have names ({.field obs_names}/{.field var_names} are NULL)",
        "i" = "Use integer indices instead, or set names first"
      ))
    }
    indices <- match(subset, names_vector)
    if (any(is.na(indices))) {
      # nolint start: object_usage_linter
      missing_names <- subset[is.na(indices)]
      # nolint end: object_usage_linter
      cli_abort(
        "Names of {context_name} not found: {paste(missing_names, collapse = ', ')}"
      )
    }

    indices
  } else {
    cli_abort(
      "Subset of {context_name} must be logical, integer, or character vector"
    )
  }
}
