#' @title AnnDataView
#'
#' @description
#' A lazy view of an AnnData object that allows applying transformations
#' (subsetting, renaming keys, or transformation functions) without immediately
#' executing them. The transformations are applied when converting to a concrete
#' AnnData implementation (InMemoryAnnData, HDF5AnnData) or other formats
#' (SingleCellExperiment, Seurat).
#'
#' This class is useful for chaining multiple operations on large datasets
#' without creating intermediate copies of the data.
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
#' # Create a view with lazy transformations using AbstractAnnData methods
#' view <- ad$
#'   subset_obs(cell_type == "A")$
#'   subset_var(gene_type %in% c("X", "Y"))$
#'   rename_obs(cell_category = "cell_type")
#'
#' # Apply transformations by converting to concrete implementation
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
    .obs_renames = NULL,
    .var_renames = NULL,
    .obsm_renames = NULL,
    .varm_renames = NULL,
    .obsp_renames = NULL,
    .varp_renames = NULL,
    .layers_renames = NULL,
    .uns_renames = NULL,
    .transformation_functions = NULL,

    # Apply obs subsetting to a data frame or vector
    .apply_obs_subset = function(x) {
      if (is.null(private$.obs_subset) || is.null(x)) {
        return(x)
      }
      if (is.data.frame(x)) {
        return(x[private$.obs_subset, , drop = FALSE])
      } else if (is.matrix(x)) {
        return(x[private$.obs_subset, , drop = FALSE])
      } else if (is.vector(x)) {
        return(x[private$.obs_subset])
      }
      return(x)
    },

    # Apply var subsetting to a data frame or matrix
    .apply_var_subset = function(x) {
      if (is.null(private$.var_subset) || is.null(x)) {
        return(x)
      }
      if (is.data.frame(x)) {
        return(x[private$.var_subset, , drop = FALSE])
      } else if (is.matrix(x)) {
        return(x[, private$.var_subset, drop = FALSE])
      } else if (is.vector(x)) {
        return(x[private$.var_subset])
      }
      return(x)
    },

    # Apply both obs and var subsetting to matrices
    .apply_obs_var_subset = function(x) {
      if (is.null(x)) {
        return(x)
      }

      result <- x
      if (!is.null(private$.obs_subset)) {
        result <- result[private$.obs_subset, , drop = FALSE]
      }
      if (!is.null(private$.var_subset)) {
        result <- result[, private$.var_subset, drop = FALSE]
      }
      return(result)
    },

    # Apply renaming to a named list
    .apply_renames = function(x, rename_map) {
      if (is.null(x) || is.null(rename_map)) {
        return(x)
      }

      if (is.list(x)) {
        old_names <- names(x)
        new_names <- old_names
        for (old_name in names(rename_map)) {
          if (old_name %in% old_names) {
            new_names[old_names == old_name] <- rename_map[[old_name]]
          }
        }
        names(x) <- new_names
      } else if (is.data.frame(x)) {
        old_names <- colnames(x)
        new_names <- old_names
        for (old_name in names(rename_map)) {
          if (old_name %in% old_names) {
            new_names[old_names == old_name] <- rename_map[[old_name]]
          }
        }
        colnames(x) <- new_names
      }

      return(x)
    },

    # Apply transformation functions
    .apply_transformations = function(slot_name, x) {
      if (is.null(private$.transformation_functions) || is.null(x)) {
        return(x)
      }

      transformations <- private$.transformation_functions[[slot_name]]
      if (is.null(transformations)) {
        return(x)
      }

      for (transform_func in transformations) {
        x <- transform_func(x)
      }

      return(x)
    }
  ),
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set X on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      x <- private$.base_adata$X
      x <- private$.apply_obs_var_subset(x)
      x <- private$.apply_transformations("X", x)
      return(x)
    },

    #' @field layers See [AnnData-usage]
    layers = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set layers on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      layers <- private$.base_adata$layers
      if (!is.null(layers)) {
        # Apply subsetting to each layer
        layers <- lapply(layers, private$.apply_obs_var_subset)
        # Apply renaming
        layers <- private$.apply_renames(layers, private$.layers_renames)
        # Apply transformations
        layers <- private$.apply_transformations("layers", layers)
      }
      return(layers)
    },

    #' @field obs See [AnnData-usage]
    obs = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set obs on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      obs <- private$.base_adata$obs
      obs <- private$.apply_obs_subset(obs)
      obs <- private$.apply_renames(obs, private$.obs_renames)
      obs <- private$.apply_transformations("obs", obs)
      return(obs)
    },

    #' @field var See [AnnData-usage]
    var = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set var on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      var <- private$.base_adata$var
      var <- private$.apply_var_subset(var)
      var <- private$.apply_renames(var, private$.var_renames)
      var <- private$.apply_transformations("var", var)
      return(var)
    },

    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set obs_names on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      names <- private$.base_adata$obs_names
      names <- private$.apply_obs_subset(names)
      return(names)
    },

    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set var_names on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      names <- private$.base_adata$var_names
      names <- private$.apply_var_subset(names)
      return(names)
    },

    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set obsm on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      obsm <- private$.base_adata$obsm
      if (!is.null(obsm)) {
        # Apply obs subsetting to each matrix
        obsm <- lapply(obsm, private$.apply_obs_subset)
        # Apply renaming
        obsm <- private$.apply_renames(obsm, private$.obsm_renames)
        # Apply transformations
        obsm <- private$.apply_transformations("obsm", obsm)
      }
      return(obsm)
    },

    #' @field varm See [AnnData-usage]
    varm = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set varm on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      varm <- private$.base_adata$varm
      if (!is.null(varm)) {
        # Apply var subsetting to each matrix
        varm <- lapply(varm, private$.apply_var_subset)
        # Apply renaming
        varm <- private$.apply_renames(varm, private$.varm_renames)
        # Apply transformations
        varm <- private$.apply_transformations("varm", varm)
      }
      return(varm)
    },

    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set obsp on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      obsp <- private$.base_adata$obsp
      if (!is.null(obsp)) {
        # Apply obs subsetting to each matrix (both dimensions)
        obsp <- lapply(obsp, function(x) {
          if (!is.null(private$.obs_subset) && !is.null(x)) {
            x[private$.obs_subset, private$.obs_subset, drop = FALSE]
          } else {
            x
          }
        })
        # Apply renaming
        obsp <- private$.apply_renames(obsp, private$.obsp_renames)
        # Apply transformations
        obsp <- private$.apply_transformations("obsp", obsp)
      }
      return(obsp)
    },

    #' @field varp See [AnnData-usage]
    varp = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set varp on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      varp <- private$.base_adata$varp
      if (!is.null(varp)) {
        # Apply var subsetting to each matrix (both dimensions)
        varp <- lapply(varp, function(x) {
          if (!is.null(private$.var_subset) && !is.null(x)) {
            x[private$.var_subset, private$.var_subset, drop = FALSE]
          } else {
            x
          }
        })
        # Apply renaming
        varp <- private$.apply_renames(varp, private$.varp_renames)
        # Apply transformations
        varp <- private$.apply_transformations("varp", varp)
      }
      return(varp)
    },

    #' @field uns See [AnnData-usage]
    uns = function(value) {
      if (!missing(value)) {
        cli_abort(
          "Cannot set uns on a AnnDataView object. Convert to a concrete implementation first."
        )
      }

      uns <- private$.base_adata$uns
      uns <- private$.apply_renames(uns, private$.uns_renames)
      uns <- private$.apply_transformations("uns", uns)
      return(uns)
    }
  ),
  public = list(
    #' @description
    #' Create a new AnnDataView object
    #'
    #' @param base_adata An existing AnnData object to create a view of
    #'
    #' @return A new `AnnDataView` object
    initialize = function(base_adata) {
      if (!inherits(base_adata, "AbstractAnnData")) {
        cli_abort("base_adata must be an AbstractAnnData object")
      }

      private$.base_adata <- base_adata
      private$.obs_subset <- NULL
      private$.var_subset <- NULL
      private$.obs_renames <- list()
      private$.var_renames <- list()
      private$.obsm_renames <- list()
      private$.varm_renames <- list()
      private$.obsp_renames <- list()
      private$.varp_renames <- list()
      private$.layers_renames <- list()
      private$.uns_renames <- list()
      private$.transformation_functions <- list()

      invisible(self)
    },

    #' @description
    #' Subset observations based on a logical condition
    #'
    #' @param condition A logical condition to apply to observations.
    #'   Can be a logical vector of length n_obs or an expression that
    #'   evaluates to such a vector in the context of the obs data frame.
    #'
    #' @return The AnnDataView object (for method chaining)
    subset_obs = function(condition) {
      if (is.logical(condition)) {
        if (length(condition) != self$n_obs()) {
          cli_abort("Logical condition must have length equal to n_obs")
        }
        subset_indices <- which(condition)
      } else {
        # Evaluate condition in the context of obs
        obs_data <- private$.base_adata$obs
        if (!is.null(private$.obs_subset)) {
          obs_data <- obs_data[private$.obs_subset, , drop = FALSE]
        }

        condition_result <- eval(
          substitute(condition),
          obs_data,
          parent.frame()
        )
        if (!is.logical(condition_result)) {
          cli_abort("Condition must evaluate to a logical vector")
        }

        current_indices <- private$.obs_subset %||%
          seq_len(private$.base_adata$n_obs())
        subset_indices <- current_indices[which(condition_result)]
      }

      private$.obs_subset <- subset_indices
      invisible(self)
    },

    #' @description
    #' Subset variables based on a logical condition
    #'
    #' @param condition A logical condition to apply to variables.
    #'   Can be a logical vector of length n_vars or an expression that
    #'   evaluates to such a vector in the context of the var data frame.
    #'
    #' @return The AnnDataView object (for method chaining)
    subset_var = function(condition) {
      if (is.logical(condition)) {
        if (length(condition) != self$n_vars()) {
          cli_abort("Logical condition must have length equal to n_vars")
        }
        subset_indices <- which(condition)
      } else {
        # Evaluate condition in the context of var
        var_data <- private$.base_adata$var
        if (!is.null(private$.var_subset)) {
          var_data <- var_data[private$.var_subset, , drop = FALSE]
        }

        condition_result <- eval(
          substitute(condition),
          var_data,
          parent.frame()
        )
        if (!is.logical(condition_result)) {
          cli_abort("Condition must evaluate to a logical vector")
        }

        current_indices <- private$.var_subset %||%
          seq_len(private$.base_adata$n_vars())
        subset_indices <- current_indices[which(condition_result)]
      }

      private$.var_subset <- subset_indices
      invisible(self)
    },

    #' @description
    #' Rename columns in obs
    #'
    #' @param ... Named arguments where names are new column names and values are old column names
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_obs = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(invisible(self))
      }

      # Reverse the mapping: old_name -> new_name
      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.obs_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename columns in var
    #'
    #' @param ... Named arguments where names are new column names and values are old column names
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_var = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(self)
      }

      # Reverse the mapping: old_name -> new_name
      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.var_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename keys in obsm
    #'
    #' @param ... Named arguments where names are new keys and values are old keys
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_obsm = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(invisible(self))
      }

      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.obsm_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename keys in varm
    #'
    #' @param ... Named arguments where names are new keys and values are old keys
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_varm = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(invisible(self))
      }

      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.varm_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename keys in obsp
    #'
    #' @param ... Named arguments where names are new keys and values are old keys
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_obsp = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(self)
      }

      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.obsp_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename keys in varp
    #'
    #' @param ... Named arguments where names are new keys and values are old keys
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_varp = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(self)
      }

      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.varp_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename keys in layers
    #'
    #' @param ... Named arguments where names are new keys and values are old keys
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_layers = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(self)
      }

      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.layers_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Rename keys in uns
    #'
    #' @param ... Named arguments where names are new keys and values are old keys
    #'
    #' @return The AnnDataView object (for method chaining)
    rename_uns = function(...) {
      renames <- list(...)
      if (length(renames) == 0) {
        return(self)
      }

      for (new_name in names(renames)) {
        old_name <- renames[[new_name]]
        private$.uns_renames[[old_name]] <- new_name
      }

      invisible(self)
    },

    #' @description
    #' Add a transformation function to be applied to a specific slot
    #'
    #' @param slot_name The name of the slot to apply the transformation to
    #'   (e.g., "X", "obs", "var", "layers", "obsm", "varm", "obsp", "varp", "uns")
    #' @param transform_func A function that takes the slot data and returns transformed data
    #'
    #' @return The AnnDataView object (for method chaining)
    add_transformation = function(slot_name, transform_func) {
      if (!is.function(transform_func)) {
        cli_abort("transform_func must be a function")
      }

      valid_slots <- c(
        "X",
        "obs",
        "var",
        "layers",
        "obsm",
        "varm",
        "obsp",
        "varp",
        "uns"
      )
      if (!slot_name %in% valid_slots) {
        cli_abort(
          "slot_name must be one of: {paste(valid_slots, collapse = ', ')}"
        )
      }

      if (is.null(private$.transformation_functions[[slot_name]])) {
        private$.transformation_functions[[slot_name]] <- list()
      }

      private$.transformation_functions[[slot_name]] <- append(
        private$.transformation_functions[[slot_name]],
        transform_func
      )

      invisible(self)
    }
  )
)
