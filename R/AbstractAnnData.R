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
AbstractAnnData <- R6::R6Class("AbstractAnnData", # nolint
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
    #'   the variable dimension of the AnnData object. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    var_names = function(value) {
      .abstract_function("ad$var_names")
    },
    #' @field obsm The obsm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `obs`.
    obsm = function(value) {
      .abstract_function("ad$obsm")
    },
    #' @field varm The varm slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `var`.
    varm = function(value) {
      .abstract_function("ad$varm")
    },
    #' @field obsp The obsp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `obs`.
    obsp = function(value) {
      .abstract_function("ad$obsp")
    },
    #' @field varp The varp slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `var`.
    varp = function(value) {
      .abstract_function("ad$varp")
    },
    #' @field uns The uns slot. Must be `NULL` or a named list.
    uns = function(value) {
      .abstract_function("ad$uns")
    }
  ),
  public = list(
    #' @description Print a summary of the AnnData object. `print()`
    #'   methods should be implemented so that they are not
    #'   computationally expensive.
    #' @param ... Optional arguments to print method.
    print = function(...) {
      cat("AnnData object with n_obs \u00D7 n_vars = ", self$n_obs(), " \u00D7 ", self$n_vars(), "\n", sep = "")

      for (attribute in c(
        "obs",
        "var",
        # "uns", # TODO: remove this when uns is implemented
        "obsm",
        "varm",
        "layers",
        "obsp",
        "varp"
      )) {
        key_fun <- self[[paste0(attribute, "_keys")]]
        keys <-
          if (!is.null(key_fun)) {
            key_fun()
          } else {
            NULL
          }
        if (length(keys) > 0) {
          cat(
            "    ", attribute, ":",
            paste("'", keys, "'", collapse = ", "),
            "\n",
            sep = ""
          )
        }
      }
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
      length(self$obs_names)
    },
    #' @description Number of variables in the AnnData object.
    n_vars = function() {
      length(self$var_names)
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
    #' @description Keys (element names) of `obsm`.
    obsm_keys = function() {
      names(self$obsm)
    },
    #' @description Keys (element names) of `varm`.
    varm_keys = function() {
      names(self$varm)
    },
    #' @description Keys (element names) of `obsp`.
    obsp_keys = function() {
      names(self$obsp)
    },
    #' @description Keys (element names) of `varp`.
    varp_keys = function() {
      names(self$varp)
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
    },
    #' @description Convert to an HDF5 Backed AnnData
    #' @param path The path to the HDF5 file
    to_HDF5AnnData = function(path) {
      to_HDF5AnnData(self, path)
    }
  ),
  private = list(
    # @description `.validate_aligned_array()` checks that dimensions are
    #   consistent with the anndata object.
    # @param mat A matrix to validate
    # @param label Must be `"X"` or `"layer[[...]]"` where `...` is
    #   the name of a layer.
    # @param shape Expected dimensions of matrix
    # @param expected_rownames
    # @param excepted_colnames
    .validate_aligned_array = function(mat, label, shape, expected_rownames = NULL, expected_colnames = NULL) {
      if (is.null(mat)) {
        return(mat)
      }
      mat_dims <- dim(mat)
      for (i in seq_along(shape)) {
        expected_dim <- shape[i]
        found_dim <- mat_dims[i]
        if (found_dim != expected_dim) {
          stop("dim(", label, ")[", i, "] should have shape: ", expected_dim, ", found: ", found_dim, ".")
        }
      }
      if (!is.null(expected_rownames) & !has_row_names(mat)) {
        if (!identical(rownames(mat), expected_rownames)) {
          stop("rownames(", label, ") should be the same as expected_rownames")
        }
        rownames(mat) <- NULL
      }
      if (!is.null(expected_colnames) & !is.null(colnames(mat))) {
        if (!identical(colnames(mat), expected_colnames)) {
          stop("colnames(", label, ") should be the same as expected_colnames")
        }
        colnames(mat) <- NULL
      }

      mat
    },
    # @description `.validate_aligned_mapping()` checks for named lists and
    #   correct dimensions on elements.
    # @param collection A named list of 0 or more matrix elements with
    #   whose entries will be validated
    # @param label The label of the collection, used for error messages
    # @param shape Expected dimensions of arrays. Arrays may have more dimensions than specified here
    # @param expected_rownames Expected row names
    # @param expected_colnames Expected column names
    .validate_aligned_mapping = function(collection, label, shape, expected_rownames = NULL, expected_colnames = NULL) {
      if (is.null(collection)) {
        return(collection)
      }

      collection_names <- names(collection)
      if (!is.list(collection) || ((length(collection) != 0) && is.null(collection_names))) {
        stop(paste0(label, " must be a named list, was ", class(collection)))
      }

      for (mtx_name in collection_names) {
        collection_name <- paste0(label, "[['", mtx_name, "']]")
        private$.validate_aligned_array(
          collection[[mtx_name]],
          collection_name,
          shape = shape,
          expected_rownames = expected_rownames,
          expected_colnames = expected_colnames
        )
      }

      collection
    },

    # @description `.validate_named_list()` checks for whether a value
    #   is NULL or a named list and throws an error if it is not.
    .validate_named_list = function(collection, label) {
      if (is.null(collection)) {
        return(collection)
      }

      collection_names <- names(collection)
      if (!is.list(collection) || ((length(collection) != 0) && is.null(collection_names))) {
        stop(paste0(label, " must be a named list, was ", class(collection)))
      }

      collection
    },

    # @description `.validate_obsvar_dataframe()` checks that the
    #   object is a data.frame and removes explicit dimnames.
    # @param df A data frame to validate. Should be an obs or a var.
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_dataframe = function(df, label = c("obs", "var")) {
      label <- match.arg(label)

      expected_nrow <- switch(label,
        obs = self$n_obs(),
        var = self$n_vars()
      )

      if (is.null(df)) {
        # create empty data frame
        df <- data.frame(i = seq_len(expected_nrow))[, -1, drop = FALSE]
      }

      if (!is.data.frame(df)) {
        stop(label, " should be a data frame")
      }

      if (nrow(df) != expected_nrow) {
        stop(wrap_message(
          "nrow(df) should match the number of ", label, ". ",
          "Expected nrow: ", expected_nrow, ". ",
          "Observed nrow: ", nrow(df), "."
        ))
      }

      if (has_row_names(df)) {
        warning(wrap_message(
          "'", label, "' should not have any rownames, removing them from the data frame."
        ))
        rownames(df) <- NULL
      }

      df
    },

    # @description `.validate_obsvar_names()` checks that `*_names()`
    #   are NULL or consistent with the dimensions of `obs` or `var`.
    # @param names A vector to validate
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_names = function(names, label = c("obs", "var"), check_length = TRUE) {
      label <- match.arg(label)

      if (is.null(names)) {
        stop(wrap_message(label, "_names should be defined."))
      }

      # only check whether sizes match if the obsvar names has already been defined
      prev_names <- attr(self, paste0(label, "_names"))
      if (!is.null(prev_names)) {
        expected_len <- length(prev_names)

        if (length(names) != expected_len) {
          size_check_label <- if (label == "obs") "n_obs" else "n_vars"
          stop(wrap_message(
            "length(", label, "_names) should be the same as ad$", size_check_label, "()"
          ))
        }
      }

      names
    }
  )
)
