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
#' This class is an abstract representation of an AnnData object. It is
#' intended to be used as a base class for concrete implementations of
#' AnnData objects, such as [InMemoryAnnData] or [HDF5AnnData].
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @seealso
#'
#' - [AnnData-usage]: For details on creating and using  `AnnData` objects
#' - [InMemoryAnnData]: For the in-memory `AnnData` class
#' - [HDF5AnnData]: For the HDF5-backed `AnnData` class
AbstractAnnData <- R6::R6Class(
  "AbstractAnnData",
  active = list(
    #' @field X `NULL` or an observation x variable matrix (without
    #'   dimnames) consistent with the number of rows of `obs` and `var`.
    X = function(value) {
      .abstract_function("ad$X")
    },
    #' @field layers The `layers` slot. Must be `NULL` or a named list
    #'   with with all elements having the dimensions consistent with
    #'   `obs` and `var`.
    layers = function(value) {
      .abstract_function("ad$layers")
    },
    #' @field obs A `data.frame` with columns containing information
    #'   about observations. The number of rows of `obs` defines the
    #'   observation dimension of the `AnnData` object.
    obs = function(value) {
      .abstract_function("ad$obs")
    },
    #' @field var A `data.frame` with columns containing information
    #'   about variables. The number of rows of `var` defines the variable
    #'   dimension of the `AnnData` object.
    var = function(value) {
      .abstract_function("ad$var")
    },
    #' @field obs_names Either `NULL` or a vector of unique identifiers
    #'   used to identify each row of `obs` and to act as an index into
    #'   the observation dimension of the `AnnData` object. For
    #'   compatibility with *R* representations, `obs_names` should be a
    #'   character vector.
    obs_names = function(value) {
      .abstract_function("ad$obs_names")
    },
    #' @field var_names Either `NULL` or a vector of unique identifiers
    #'   used to identify each row of `var` and to act as an index into
    #'   the variable dimension of the `AnnData` object. For compatibility
    #'   with *R* representations, `var_names` should be a character
    #'   vector.
    var_names = function(value) {
      .abstract_function("ad$var_names")
    },
    #' @field obsm The `obsm` slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `obs`.
    obsm = function(value) {
      .abstract_function("ad$obsm")
    },
    #' @field varm The `varm` slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows as `var`.
    varm = function(value) {
      .abstract_function("ad$varm")
    },
    #' @field obsp The `obsp` slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `obs`.
    obsp = function(value) {
      .abstract_function("ad$obsp")
    },
    #' @field varp The `varp` slot. Must be `NULL` or a named list with
    #'   with all elements having the same number of rows and columns as `var`.
    varp = function(value) {
      .abstract_function("ad$varp")
    },
    #' @field uns The `uns` slot. Must be `NULL` or a named list.
    uns = function(value) {
      .abstract_function("ad$uns")
    }
  ),
  public = list(
    #' @description
    #' Print a summary of the AnnData object. `print()` methods should be
    #' implemented so that they are not computationally expensive.
    #'
    #' @param ... Optional arguments to print method.
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

    #' @description Dimensions (observations x variables) of the `AnnData` object.
    shape = function() {
      c(
        self$n_obs(),
        self$n_vars()
      )
    },
    #' @description Number of observations in the `AnnData` object.
    n_obs = function() {
      nrow(self$obs)
    },
    #' @description Number of variables in the `AnnData` object.
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
    #' @description Keys (element names) of `uns`.
    uns_keys = function() {
      names(self$uns)
    },
    #' @description
    #' Convert to `SingleCellExperiment`
    #'
    #' See [to_SingleCellExperiment()] for more details on the conversion.
    #'
    #' @param assays_mapping A named list mapping `AnnData` layers to
    #'   `SingleCellExperiment` assays
    #' @param colData_mapping A named list mapping `AnnData` `obs` to
    #'   `SingleCellExperiment` `colData`
    #' @param rowData_mapping A named list mapping `AnnData` `var` to
    #'   `SingleCellExperiment` `rowData`
    #' @param reduction_mapping A named list mapping `AnnData` `obsm`/`varm` to
    #'   `SingleCellExperiment` reductions
    #' @param colPairs_mapping A named list mapping `AnnData` `obsp`/`varp` to
    #'   `SingleCellExperiment` `colPairs`
    #' @param rowPairs_mapping A named list mapping `AnnData` `obsp`/`varp` to
    #'   `SingleCellExperiment` `rowPairs`
    #' @param metadata_mapping A named list mapping `AnnData` `uns` to
    #'   `SingleCellExperiment` `metadata`
    #'
    #' @return A `SingleCellExperiment` object
    to_SingleCellExperiment = function(
      assays_mapping = NULL,
      colData_mapping = NULL, # nolint
      rowData_mapping = NULL, # nolint
      reduction_mapping = NULL,
      colPairs_mapping = NULL, # nolint
      rowPairs_mapping = NULL, # nolint
      metadata_mapping = NULL
    ) {
      to_SingleCellExperiment(
        self,
        assays_mapping = assays_mapping,
        colData_mapping = colData_mapping,
        rowData_mapping = rowData_mapping,
        reduction_mapping = reduction_mapping,
        colPairs_mapping = colPairs_mapping, # nolint
        rowPairs_mapping = rowPairs_mapping, # nolint
        metadata_mapping = metadata_mapping
      )
    },
    #' @description
    #' Convert to `Seurat`
    #'
    #' See [to_Seurat()] for more details on the conversion and each of the
    #' parameters.
    #'
    #' @param assay_name The name of the `assay` to store data in
    #' @param layers_mapping A named list mapping `Seurat` layers to `AnnData`
    #'   layers
    #' @param object_metadata_mapping A named list mapping `AnnData` `obs` to
    #'   `Seurat` `metadata`
    #' @param assay_metadata_mapping A named list mapping `AnnData` `var` to
    #'   `Seurat` assay metadata
    #' @param reduction_mapping A named list mapping `AnnData` `obsm`/`varm`
    #'   to `Seurat` reductions
    #' @param graph_mapping A named list mapping `AnnData` `obsp`/`varp` to
    #'   `Seurat` graphs
    #' @param misc_mapping A named list mapping `AnnData` `uns` to `Seurat`
    #'   `misc`
    #'
    #' @return A `Seurat` object
    to_Seurat = function(
      assay_name = "RNA",
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
        layers_mapping = layers_mapping,
        object_metadata_mapping = object_metadata_mapping,
        assay_metadata_mapping = assay_metadata_mapping,
        reduction_mapping = reduction_mapping,
        graph_mapping = graph_mapping,
        misc_mapping = misc_mapping
      )
    },
    #' @description Convert to an in-memory `AnnData`
    #'
    #' @return An [InMemoryAnnData] object
    to_InMemoryAnnData = function() {
      to_InMemoryAnnData(self)
    },
    #' @description
    #' Convert to an HDF5-backed `AnnData`
    #'
    #' @param file The path to the HDF5 file
    #' @param compression The compression algorithm to use when writing the
    #'   HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
    #'   `"none"`.
    #' @param mode The mode to open the HDF5 file.
    #'   * `a` creates a new file or opens an existing one for read/write.
    #'   * `r+` opens an existing file for read/write.
    #'   * `w` creates a file, truncating any existing ones
    #'   * `w-`/`x` are synonyms creating a file and failing if it already
    #'     exists.
    #'
    #' @return An [HDF5AnnData] object
    to_HDF5AnnData = function(
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
    #' Write the `AnnData` object to an H5AD file.
    #'
    #' @param path The path to the H5AD file
    #' @param compression The compression algorithm to use when writing the
    #'   HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
    #'   `"none"`.
    #' @param mode The mode to open the HDF5 file.
    #'   * `a` creates a new file or opens an existing one for read/write.
    #'   * `r+` opens an existing file for read/write.
    #'   * `w` creates a file, truncating any existing ones
    #'   * `w-`/`x` are synonyms creating a file and failing if it already
    #'     exists.
    #'
    #' @return `path` invisibly
    #'
    #' @examples
    #' adata <- AnnData(
    #'   X = matrix(1:5, 3L, 5L),
    #'   layers = list(
    #'     A = matrix(5:1, 3L, 5L),
    #'     B = matrix(letters[1:5], 3L, 5L)
    #'   ),
    #'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
    #'   var = data.frame(row.names = letters[1:5], gene = 1:5)
    #' )
    #' h5ad_file <- tempfile(fileext = ".h5ad")
    #' adata$write_h5ad(h5ad_file)
    write_h5ad = function(path, compression = "gzip", mode = "w") {
      self$to_HDF5AnnData(path, compression = compression, mode = mode)
      path
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
    # @param expected_rownames
    # @param expected_colnames
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
