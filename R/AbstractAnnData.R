.abstract_function <- function(fun_name = NA_character_) {
  fun_str <- if (is.na(fun_name)) {
    "This function"
  } else {
    "{.fun {fun_name}}"
  }
  cli_abort(
    paste(fun_str, "is an abstract function."),
    call = rlang::caller_env()
  )
}

.anndata_slots <- c(
  "X",
  "obs",
  "var",
  "uns",
  "obsm",
  "varm",
  "layers",
  "obsp",
  "varp"
)

#' @title Abstract AnnData class
#'
#' @description
#' This class is an abstract representation of an `AnnData` object. It is
#' intended to be used as a base class for concrete implementations of
#' `AnnData` objects, such as [InMemoryAnnData] or [HDF5AnnData].
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
AbstractAnnData <- R6::R6Class(
  "AbstractAnnData",
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      .abstract_function("ad$X")
    },
    #' @field layers See [AnnData-usage]
    layers = function(value) {
      .abstract_function("ad$layers")
    },
    #' @field obs See [AnnData-usage]
    obs = function(value) {
      .abstract_function("ad$obs")
    },
    #' @field var See [AnnData-usage]
    var = function(value) {
      .abstract_function("ad$var")
    },
    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      .abstract_function("ad$obs_names")
    },
    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      .abstract_function("ad$var_names")
    },
    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      .abstract_function("ad$obsm")
    },
    #' @field varm See [AnnData-usage]
    varm = function(value) {
      .abstract_function("ad$varm")
    },
    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      .abstract_function("ad$obsp")
    },
    #' @field varp See [AnnData-usage]
    varp = function(value) {
      .abstract_function("ad$varp")
    },
    #' @field uns See [AnnData-usage]
    uns = function(value) {
      .abstract_function("ad$uns")
    }
  ),
  public = list(
    #' @description See [AnnData-usage]
    #'
    #' @param ... Optional arguments to print method
    print = function(...) {
      cat(
        "AnnData object with n_obs \u00D7 n_vars = ",
        self$n_obs(),
        " \u00D7 ",
        self$n_vars(),
        "\n",
        sep = ""
      )

      for (attribute in .anndata_slots[-1]) {
        key_fun <- self[[paste0(attribute, "_keys")]]
        keys <-
          if (!is.null(key_fun)) {
            key_fun()
          } else {
            NULL
          }
        if (length(keys) > 0) {
          cat(
            "    ",
            attribute,
            ": ",
            paste(paste0("'", keys, "'"), collapse = ", "),
            "\n",
            sep = ""
          )
        }
      }
    },

    #' @description See [AnnData-usage]
    shape = function() {
      c(
        self$n_obs(),
        self$n_vars()
      )
    },
    #' @description See [AnnData-usage]
    n_obs = function() {
      nrow(self$obs)
    },
    #' @description See [AnnData-usage]
    n_vars = function() {
      nrow(self$var)
    },
    #' @description See [AnnData-usage]
    obs_keys = function() {
      names(self$obs)
    },
    #' @description See [AnnData-usage]
    var_keys = function() {
      names(self$var)
    },
    #' @description See [AnnData-usage]
    layers_keys = function() {
      names(self$layers)
    },
    #' @description See [AnnData-usage]
    obsm_keys = function() {
      names(self$obsm)
    },
    #' @description See [AnnData-usage]
    varm_keys = function() {
      names(self$varm)
    },
    #' @description See [AnnData-usage]
    obsp_keys = function() {
      names(self$obsp)
    },
    #' @description See [AnnData-usage]
    varp_keys = function() {
      names(self$varp)
    },
    #' @description See [AnnData-usage]
    uns_keys = function() {
      names(self$uns)
    },
    #' @description
    #' Convert to `SingleCellExperiment`
    #'
    #' See [as_SingleCellExperiment()] for more details on the conversion
    #'
    #' @param x_mapping See [as_SingleCellExperiment()]
    #' @param assays_mapping See [as_SingleCellExperiment()]
    #' @param colData_mapping See [as_SingleCellExperiment()]
    #' @param rowData_mapping See [as_SingleCellExperiment()]
    #' @param reducedDims_mapping See [as_SingleCellExperiment()]
    #' @param colPairs_mapping See [as_SingleCellExperiment()]
    #' @param rowPairs_mapping See [as_SingleCellExperiment()]
    #' @param metadata_mapping See [as_SingleCellExperiment()]
    #'
    #' @return A `SingleCellExperiment` object
    # nolint start: object_name_linter
    as_SingleCellExperiment = function(
      x_mapping = NULL,
      assays_mapping = NULL,
      colData_mapping = NULL,
      rowData_mapping = NULL,
      reducedDims_mapping = NULL,
      colPairs_mapping = NULL,
      rowPairs_mapping = NULL,
      metadata_mapping = NULL
    ) {
      # nolint end: object_name_linter
      to_SingleCellExperiment(
        self,
        x_mapping = x_mapping,
        assays_mapping = assays_mapping,
        colData_mapping = colData_mapping,
        rowData_mapping = rowData_mapping,
        reducedDims_mapping = reducedDims_mapping,
        colPairs_mapping = colPairs_mapping, # nolint
        rowPairs_mapping = rowPairs_mapping, # nolint
        metadata_mapping = metadata_mapping
      )
    },
    #' @description
    #' `r lifecycle::badge('deprecated')`
    #'
    #' Deprecated, use `as_SingleCellExperiment()` instead
    #'
    #' @param ... Arguments passed to `adata$as_SingleCellExperiment()`
    #'
    #' @return A `SingleCellExperiment` object
    to_SingleCellExperiment = function(...) {
      lifecycle::deprecate_warn(
        "0.99.0",
        "to_SingleCellExperiment()",
        "as_SingleCellExperiment()"
      )

      self$as_SingleCellExperiment(...)
    },
    #' @description
    #' Convert to `Seurat`
    #'
    #' See [as_Seurat()] for more details on the conversion
    #'
    #' @param assay_name See [as_Seurat()]
    #' @param x_mapping See [as_Seurat()]
    #' @param layers_mapping See [as_Seurat()]
    #' @param object_metadata_mapping See [as_Seurat()]
    #' @param assay_metadata_mapping See [as_Seurat()]
    #' @param reduction_mapping See [as_Seurat()]
    #' @param graph_mapping See [as_Seurat()]
    #' @param misc_mapping See [as_Seurat()]
    #'
    #' @return A `Seurat` object
    as_Seurat = function(
      assay_name = "RNA",
      x_mapping = NULL,
      layers_mapping = NULL,
      object_metadata_mapping = NULL,
      assay_metadata_mapping = NULL,
      reduction_mapping = NULL,
      graph_mapping = NULL,
      misc_mapping = NULL
    ) {
      to_Seurat(
        self,
        assay_name = assay_name,
        x_mapping = x_mapping,
        layers_mapping = layers_mapping,
        object_metadata_mapping = object_metadata_mapping,
        assay_metadata_mapping = assay_metadata_mapping,
        reduction_mapping = reduction_mapping,
        graph_mapping = graph_mapping,
        misc_mapping = misc_mapping
      )
    },
    #' @description
    #' `r lifecycle::badge('deprecated')`
    #'
    #' Deprecated, use `as_Seurat()` instead
    #'
    #' @param ... Arguments passed to `adata$as_Seurat()`
    #'
    #' @return A `Seurat` object
    to_Seurat = function(...) {
      lifecycle::deprecate_warn(
        "0.99.0",
        "to_Seurat()",
        "as_Seurat()"
      )

      self$as_Seurat(...)
    },
    #' @description
    #' Convert to an [`InMemoryAnnData`]
    #'
    #' See [as_InMemoryAnnData()] for more details on the conversion
    #'
    #' @return An [InMemoryAnnData] object
    as_InMemoryAnnData = function() {
      to_InMemoryAnnData(self)
    },
    #' @description
    #' `r lifecycle::badge('deprecated')`
    #'
    #' Deprecated, use `as_InMemoryAnnData()` instead
    #'
    #' @return An [`InMemoryAnnData`] object
    to_InMemoryAnnData = function() {
      lifecycle::deprecate_warn(
        "0.99.0",
        "to_InMemoryAnnData()",
        "as_InMemoryAnnData()"
      )

      self$as_InMemoryAnnData()
    },
    #' @description
    #' Convert to an [`HDF5AnnData`]
    #'
    #' See [as_HDF5AnnData()] for more details on the conversion
    #'
    #' @param file See [as_HDF5AnnData()]
    #' @param compression See [as_HDF5AnnData()]
    #' @param mode See [as_HDF5AnnData()]
    #'
    #' @return An [`HDF5AnnData`] object
    as_HDF5AnnData = function(
      file,
      compression = c("none", "gzip", "lzf"),
      mode = c("w-", "r", "r+", "a", "w", "x")
    ) {
      to_HDF5AnnData(
        adata = self,
        file = file,
        compression = compression,
        mode = mode
      )
    },
    #' @description
    #' `r lifecycle::badge('deprecated')`
    #'
    #' Deprecated, use `as_HDF5AnnData()` instead
    #'
    #' @param ... Arguments passed to `adata$as_HDF5AnnData()`
    #'
    #' @return An [`HDF5AnnData`] object
    to_HDF5AnnData = function(...) {
      lifecycle::deprecate_warn(
        "0.99.0",
        "to_HDF5AnnDAta()",
        "as_HDF5AnnData()",
      )

      self$as_HDF5AnnData(...)
    },
    #' @description
    #' Write the `AnnData` object to an H5AD file
    #'
    #' See [write_h5ad()] for details
    #'
    #' @param path See [write_h5ad()]
    #' @param compression See [write_h5ad()]
    #' @param mode See [write_h5ad()]
    #'
    #' @return `path` invisibly
    write_h5ad = function(
      path,
      compression = c("none", "gzip", "lzf"),
      mode = c("w-", "r", "r+", "a", "w", "x")
    ) {
      write_h5ad(object = self, path, compression = compression, mode = mode)
    }
  ),
  private = list(
    # @description `.validate_aligned_array()` checks that dimensions are
    #   consistent with the anndata object.
    #
    # @param mat A matrix to validate
    # @param label Must be `"X"` or `"layer[[...]]"` where `...` is
    #   the name of a layer.
    # @param shape Expected dimensions of matrix
    # @param expected_rownames Expected row names
    # @param expected_colnames Expected column names
    .validate_aligned_array = function(
      mat,
      label,
      shape,
      expected_rownames = NULL,
      expected_colnames = NULL
    ) {
      if (is.null(mat)) {
        return(mat)
      }
      mat_dims <- dim(mat)

      if (
        length(shape) > length(mat_dims) ||
          any(shape != mat_dims[seq_along(shape)])
      ) {
        cli_abort(
          c(
            "Unexpected shape for {.field {label}}",
            "i" = paste0(
              "Expected [",
              paste(shape, collapse = ", "),
              "], got [",
              paste(mat_dims, collapse = ", "),
              "]"
            )
          ),
          call = rlang::caller_env()
        )
      }

      if (!is.null(expected_rownames) & has_row_names(mat)) {
        if (!identical(rownames(mat), expected_rownames)) {
          cli_abort(
            c(
              "{.code rownames({label})} is not as expected",
              "i" = "Expected row names: {style_vec(expected_colnames)}",
              "i" = "Provided row names: {style_vec(rownames(mat))}"
            ),
            call = rlang::caller_env()
          )
        }
        rownames(mat) <- NULL
      }

      if (!is.null(expected_colnames) & !is.null(colnames(mat))) {
        if (!identical(colnames(mat), expected_colnames)) {
          cli_abort(
            c(
              "{.code colnames({label})} is not as expected",
              "i" = "Expected column names: {style_vec(expected_colnames)}",
              "i" = "Provided column names: {style_vec(colnames(mat))}"
            ),
            call = rlang::caller_env()
          )
        }
        colnames(mat) <- NULL
      }

      mat
    },
    # @description `.validate_aligned_mapping()` checks for named lists and
    #   correct dimensions on elements.
    #
    # @param collection A named list of 0 or more matrix elements with
    #   whose entries will be validated
    # @param label The label of the collection, used for error messages
    # @param shape Expected dimensions of arrays. Arrays may have more dimensions than specified here
    # @param expected_rownames Expected row names
    # @param expected_colnames Expected column names
    .validate_aligned_mapping = function(
      collection,
      label,
      shape,
      expected_rownames = NULL,
      expected_colnames = NULL
    ) {
      if (is.null(collection)) {
        return(collection)
      }

      collection <- private$.validate_named_list(collection, label)
      collection_names <- names(collection)

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
      if (
        !is.list(collection) ||
          ((length(collection) != 0) && is.null(collection_names))
      ) {
        cli_abort(
          "{.field {label}} must be a named {.cls list}, got {.cls {class(collection)}}",
          call = rlang::caller_env()
        )
      }

      collection
    },

    # @description `.validate_obsvar_dataframe()` checks that the
    #   object is a data.frame and removes explicit dimnames.
    #
    # @param df A data frame to validate. Should be an obs or a var.
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_dataframe = function(df, label = c("obs", "var")) {
      label <- match.arg(label)

      expected_nrow <- switch(label, obs = self$n_obs(), var = self$n_vars())

      if (is.null(df)) {
        # create empty data frame
        df <- data.frame(i = seq_len(expected_nrow))[, -1, drop = FALSE]
      }

      if (!is.data.frame(df)) {
        cli_abort(
          "{.field label} must be a {.cls data.frame}, got {.cls {class(df)}}",
          call = rlang::caller_env()
        )
      }

      if (nrow(df) != expected_nrow) {
        cli_abort(
          paste(
            "{.code nrow({label})} should equal {.val {expected_nrow}},",
            "got {.val {nrow(df)}}"
          ),
          call = rlang::caller_env()
        )
      }

      df
    },

    # @description `.validate_obsvar_names()` checks that `*_names()`
    #   are NULL or consistent with the dimensions of `obs` or `var`.
    #
    # @param names A vector to validate
    # @param label Must be `"obs"` or `"var"`
    .validate_obsvar_names = function(
      names,
      label = c("obs", "var"),
      check_length = TRUE
    ) {
      label <- match.arg(label)

      if (is.null(names)) {
        cli_abort(
          "{.field {label}_names} should be defined",
          call = rlang::caller_env()
        )
      }

      # only check whether sizes match if the obsvar names has already been defined
      prev_names <- attr(self, paste0(label, "_names"))
      if (!is.null(prev_names)) {
        expected_len <- length(prev_names)

        if (length(names) != expected_len) {
          size_check_label <- if (label == "obs") "n_obs" else "n_vars"
          cli_abort(
            paste(
              "{.code length({label}_names)} should match",
              "{.code ad${size_check_label}}"
            ),
            call = rlang::caller_env()
          )
        }
      }

      names
    }
  )
)
