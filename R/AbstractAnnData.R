.abstract_function <- function(fun_name = NA_character_) {
  fun_name <- if (is.na(fun_name)) "This function" else fun_name
  stop(fun_name, " is an abstract function.")
}

#' @title AbstractAnnData
#'
#' @description
#'   Abstract [R6][R6::R6Class] class representing an AnnData
#'   object. Defines the interface.
#' @importFrom R6 R6Class
AbstractAnnData <- R6::R6Class("AbstractAnnData",
  active = list(
    #' @field X NULL or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and `var`.
    X = function(value) {
      .abstract_function("ad$X")
    },
    #' @field layers The layers slot. Must be NULL or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      .abstract_function("ad$layers")
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the AnnData object.
    obs = function(value) {
      .abstract_function("ad$obs")
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the AnnData object.
    var = function(value) {
      .abstract_function("ad$var")
    },
    #' @field obs_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into
    #'   the observation dimension of the AnnData object. For
    #'   compatibility with *R* representations, `obs_names` should be a
    #'   character vector.
    obs_names = function(value) {
      .abstract_function("ad$obs_names")
    },
    #' @field var_names Either NULL or a vector of unique identifiers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the AnnData object.. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    var_names = function(value) {
      .abstract_function("ad$var_names")
    }
  ),
  public = list(
    #' @description Print a summary of the AnnData object. `print()`
    #'   methods should be implemented so that they are not
    #'   computationally expensive.
    #' @param ... Optional arguments to print method.
    print = function(...) {
      x_info <- if (!is.null(self$X)) {
        class(self$X)[[1]]
      } else {
        NULL
      }
      cat(
        "class: ", class(self)[[1]], "\n",
        "dim: ", self$n_obs(), " obs x ", self$n_vars(), " var\n",
        "X: ", x_info, "\n",
        pretty_print("layers", self$layers_keys()), "\n",
        pretty_print("obs", self$obs_keys()), "\n",
        pretty_print("var", self$var_keys()), "\n",
        sep = ""
      )
    },

    #' @description Dimensions (observations x variables) of the AnnData object.
    shape = function() {
      c(
        self$n_obs(),
        self$n_vars()
      )
    },
    #' @description Number of observations in the AnnData object.
    n_obs = function() {
      nrow(self$obs)
    },
    #' @description Number of variables in the AnnData object.
    n_vars = function() {
      nrow(self$var)
    },
    #' @description Keys ('column names') of `obs`.
    obs_keys = function() {
      names(self$obs)
    },
    #' @description Keys ('column names') of `var`.
    var_keys = function() {
      names(self$var)
    },
    #' @description Keys (element names) of `layers`.
    layers_keys = function() {
      names(self$layers)
    },
    #' @description Convert to SingleCellExperiment
    to_SingleCellExperiment = function() {
      to_SingleCellExperiment(self)
    },
    #' @description Convert to Seurat
    to_Seurat = function() {
      to_Seurat(self)
    },
    #' @description Convert to an InMemory AnnData
    to_InMemoryAnnData = function() {
      to_InMemoryAnnData(self)
    }
  ),
  private = list(
    # @description `.validate_matrix()` checks that dimensions are
    #   consistent with `obs` and `var`, and removes dimnames if
    #   present.
    # @param mat A matrix to validate
    # @param label Must be `"X"` or `"layer[[...]]"` where `...` is
    #   the name of a layer.
    .validate_matrix = function(mat, label) {
      if (!is.null(mat)) {
        if (nrow(mat) != nrow(self$obs)) {
          stop("nrow(", label, ") should be the same as nrow(obs)")
        }
        if (ncol(mat) != nrow(self$var)) {
          stop("ncol(", label, ") should be the same as nrow(var)")
        }

        if (!is.null(rownames(mat))) {
          warning(wrap_message(
            "rownames(", label, ") should be NULL, removing them from the ",
            "matrix"
          ))
          rownames(mat) <- NULL
        }

        if (!is.null(colnames(mat))) {
          warning(wrap_message(
            "colnames(", label, ") should be NULL, removing them from the ",
            "matrix"
          ))
          colnames(mat) <- NULL
        }
      }

      mat
    },

    # @description `.validate_layers()` checks for named lists and
    #   correct dimensions on elements.
    # @param layers A named list of 0 or more matrix elements with
    #   dimensions consistent with `obs` and `var`.
    .validate_layers = function(layers) {
      if (is.null(layers)) {
        return(layers)
      }

      ## layers and names
      layer_names <- names(layers)
      if (!is.list(layers) || is.null(layer_names)) {
        stop("'layers' must must be a named list")
      }
      if (any(!nzchar(layer_names))) {
        stop("all 'layers' elements must have non-trivial names")
      }

      ## layer elements
      for (layer in layer_names) {
        layer_name <- paste0("layers[[", layer, "]]")
        private$.validate_matrix(layers[[layer]], layer_name)
      }

      layers
    },

    # @description `.validate_obsvar_dataframe()` checks that the
    #   object is a data.frame and removes explicit dimnames.
    # @param df A data frame to validate. Should be an obs or a var.
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_dataframe = function(df, label) {
      if (is.null(df)) stop(label, " should be a data frame")
      if (.row_names_info(df) > 0) {
        warning(wrap_message(
          "'", label, "' should not have any dimnames, removing them from ",
          "the matrix"
        ))
        rownames(df) <- NULL
      }
      df
    },

    # @description `.validate_obsvar_names()` checks that `*_names()`
    #   are NULL or consistent with the dimensions of `obs` or `var`.
    # @param names A vector to validate
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_names = function(names, label) {
      if (!is.null(names) && length(names) != nrow(self[[label]])) {
        stop(wrap_message(
          "length(", label, "_names) should be the same as ",
          "nrow(", label, ")"
        ))
      }
      names
    }
  )
)
