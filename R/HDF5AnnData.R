#' @title HDF5AnnData
#'
#' @description
#' Implementation of an HDF5-backed `AnnData` object. This class provides an
#' interface to a H5AD file and minimal data is stored in memory until it is
#' requested by the user. It is primarily designed as an intermediate object
#' when reading/writing H5AD files but can be useful for accessing parts of
#' large files.
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @return An `HDF5AnnData` object
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
HDF5AnnData <- R6::R6Class(
  "HDF5AnnData", # nolint
  inherit = AbstractAnnData,
  cloneable = FALSE,
  private = list(
    .hdf5_file = NULL,
    .mode = NULL,
    .backed = NULL,
    .compression = NULL,
    .chunk_size = "auto",

    .check_mode_writeable = function() {
      if (!is.null(private$.mode) && private$.mode == "r") {
        cli_abort(
          paste(
            "The object is opened in read-only mode, cannot modify data"
          )
        )
      }
    },

    finalize = function() {
      self$close()
    }
  ),
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_X, status=done
        read_h5ad_element(private$.hdf5_file, "X", backed = private$.backed) |>
          private$.add_matrix_dimnames("X")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_X, status=done
        private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = self$obs_names,
          expected_colnames = self$var_names
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "X",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    },
    #' @field layers See [AnnData-usage]
    layers = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_layers, status=done
        read_h5ad_element(
          private$.hdf5_file,
          "layers",
          backed = private$.backed
        ) |>
          private$.add_mapping_dimnames("layers")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_layers, status=done
        private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = self$obs_names,
          expected_colnames = self$var_names
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "layers",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    },
    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsm, status=done
        read_h5ad_element(
          private$.hdf5_file,
          "obsm",
          backed = private$.backed
        ) |>
          private$.add_mapping_dimnames("obsm")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_obsm, status=done
        private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = self$obs_names,
          strip_rownames = TRUE,
          strip_colnames = FALSE,
          warn_colnames = TRUE
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "obsm",
            private$.compression,
            chunk_size = private$.chunk_size,
            index = self$obs_names
          )
      }
    },
    #' @field varm See [AnnData-usage]
    varm = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varm, status=done
        read_h5ad_element(
          private$.hdf5_file,
          "varm",
          backed = private$.backed
        ) |>
          private$.add_mapping_dimnames("varm")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_varm, status=done
        private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = self$var_names,
          strip_rownames = TRUE,
          strip_colnames = FALSE,
          warn_colnames = TRUE
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "varm",
            private$.compression,
            chunk_size = private$.chunk_size,
            index = self$var_names
          )
      }
    },
    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsp, status=done
        read_h5ad_element(
          private$.hdf5_file,
          "obsp",
          backed = private$.backed
        ) |>
          private$.add_mapping_dimnames("obsp")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_obsp, status=done
        private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = self$obs_names,
          expected_colnames = self$obs_names
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "obsp",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    },
    #' @field varp See [AnnData-usage]
    varp = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varp, status=done
        read_h5ad_element(
          private$.hdf5_file,
          "varp",
          backed = private$.backed
        ) |>
          private$.add_mapping_dimnames("varp")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_varp, status=done
        private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = self$var_names,
          expected_colnames = self$var_names
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "varp",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    },
    #' @field obs See [AnnData-usage]
    obs = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs, status=done
        read_h5ad_element(private$.hdf5_file, "obs", backed = private$.backed)
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_obs, status=done
        private$.validate_obsvar_dataframe(value, "obs") |>
          write_h5ad_element(
            private$.hdf5_file,
            "obs",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    },
    #' @field var See [AnnData-usage]
    var = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var, status=done
        read_h5ad_element(private$.hdf5_file, "var", backed = private$.backed)
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_var, status=done
        private$.validate_obsvar_dataframe(value, "var") |>
          write_h5ad_element(
            private$.hdf5_file,
            "var",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    },
    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs_names, status=done
        read_h5ad_element_keys(private$.hdf5_file, "obs", dim = "rows")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_obs_names, status=done
        write_h5ad_data_frame_index(
          value,
          private$.hdf5_file,
          "obs",
          private$.compression,
          chunk_size = private$.chunk_size
        )
      }
    },
    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var_names, status=done
        read_h5ad_element_keys(private$.hdf5_file, "var", dim = "rows")
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_var_names, status=done
        write_h5ad_data_frame_index(
          value,
          private$.hdf5_file,
          "var",
          private$.compression,
          chunk_size = private$.chunk_size
        )
      }
    },
    #' @field uns See [AnnData-usage]
    uns = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_uns, status=done
        read_h5ad_element(private$.hdf5_file, "uns", backed = private$.backed)
      } else {
        private$.check_mode_writeable()

        # trackstatus: class=HDF5AnnData, feature=set_uns, status=done
        private$.validate_named_list(
          value,
          "uns",
          warn_matrix_dimnames = TRUE
        ) |>
          write_h5ad_element(
            private$.hdf5_file,
            "uns",
            private$.compression,
            chunk_size = private$.chunk_size
          )
      }
    }
  ),
  public = list(
    #' @description
    #' `HDF5AnnData` constructor
    #'
    #' @param file The file name (character) of the `.h5ad` file. If this file
    #'   already exits, other arguments must be `NULL`.
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
    #' @param mode The mode to open the HDF5 file. See [as_HDF5AnnData()] for
    #'   details
    #' @param backed Whether the object is disk backed. See [as_HDF5AnnData()]
    #'   for details
    #' @param compression The compression algorithm to use. See
    #'   [as_HDF5AnnData()] for details
    #' @param chunk_size The target chunk size in bytes. See
    #'   [as_HDF5AnnData()] for details
    #'
    #' @details
    #' The constructor creates a new HDF5 `AnnData` interface object. This can
    #' either be used to either connect to an existing `.h5ad` file or to
    #' create a new one. If any additional slot arguments are set an existing
    #' file will be overwritten.
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
      mode = c("a", "r", "r+", "w", "w-", "x"),
      backed = FALSE,
      compression = c("none", "gzip", "lzf"),
      chunk_size = "auto"
    ) {
      check_requires("HDF5AnnData", "rhdf5", where = "Bioc")
      check_requires("HDF5AnnData", "withr", where = "CRAN")

      compression <- match.arg(compression)
      mode <- match.arg(mode)

      private$.mode <- mode
      private$.compression <- compression
      private$.chunk_size <- chunk_size

      is_readonly <- if (mode == "r") {
        TRUE
      } else {
        FALSE
      }

      # Set the auto mode
      if (mode == "a") {
        if (file.exists(file)) {
          mode <- "r+"
        } else {
          mode <- "w-"
        }
      }

      if (backed && mode != "r") {
        cli_warn(
          paste(
            "{.arg backed} can only be {.val TRUE} when {.arg mode} is {.val 'r'}.",
            "Setting {.arg backed} to {.val FALSE}."
          )
        )
        backed <- FALSE
      }
      private$.backed <- backed

      if (isTRUE(private$.backed)) {
        check_requires(
          "Reading a backed HDF5AnnData",
          c("HDF5Array", "DelayedArray"),
          where = "Bioc"
        )
      }

      # Fail is the file does not exist and in read mode
      if (!file.exists(file) && mode %in% c("r", "r+")) {
        cli_abort(
          paste(
            "File {.file {file}} does not exist but mode is set to {.val {mode}}.",
            "If you want to create a new file, use a different mode (e.g. 'w-').",
            "See {.help read_h5ad} or {.help write_h5ad} for more information."
          ),
          call = rlang::caller_env()
        )
      }

      # Fail if the file exists not allowed to overwrite
      if (file.exists(file) && mode %in% c("w-", "x")) {
        cli_abort(
          paste(
            "File {.file {file}} already exists but mode is set to {.val {mode}}.",
            "If you want to overwrite the file, use a different mode (e.g. 'w').",
            "See {.help read_h5ad} or {.help write_h5ad} for more information."
          ),
          call = rlang::caller_env()
        )
      }

      # Create/truncate the file
      if (mode %in% c("w", "w-", "x")) {
        h5file <- rhdf5::H5Fcreate(
          file,
          flags = "H5F_ACC_TRUNC",
          native = FALSE
        )
        rhdf5::H5Fclose(h5file)
      }

      if (!rhdf5::H5Fis_hdf5(file)) {
        cli_abort("File {.file {file}} is not a valid HDF5 file.")
      }

      # Set the HDF5File
      private$.hdf5_file <- HDF5File$new(file)

      is_empty <- nrow(rhdf5::h5ls(private$.hdf5_file$path)) == 0L

      if (!is_readonly) {
        if (!is_empty) {
          cli_warn(
            paste(
              "An non-empty file is opened in read/write mode.",
              "Use with caution, as this can lead to data corruption."
            )
          )
        } else {
          shape <- get_shape(obs, var, X, shape)
          obs <- get_initial_obs(obs, X, shape)
          var <- get_initial_var(var, X, shape)
          write_empty_h5ad(
            private$.hdf5_file,
            obs,
            var,
            compression,
            chunk_size
          )
        }
      }

      # File is supposed to exist by now. Check if it is a valid H5AD file
      attrs <- rhdf5::h5readAttributes(private$.hdf5_file$path, "/")
      if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
        cli_abort(c(
          "File {.file {file}} is not a valid H5AD file.",
          i = "Either the file is not an H5AD file or it was created with {.pkg anndata<0.8.0}."
        ))
      }

      if (is_readonly) {
        # if any of these variables are not NULL, throw an error
        are_null <- vapply(
          .anndata_slots,
          function(x) is.null(get(x)),
          logical(1)
        )
        if (!all(are_null)) {
          cli_abort(
            paste0(
              "Error trying to write data (",
              paste(.anndata_slots[!are_null], collapse = ", "),
              ") to an H5AD file opened in read-only mode."
            )
          )
        }
      } else {
        # Otherwise, write data to slots
        for (slot in .anndata_slots) {
          value <- get(slot)
          if (!is.null(value)) {
            self[[slot]] <- value
          }
        }
      }

      self
    },

    #' @description See [AnnData-usage]
    obs_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "obs", dim = "cols")
    },
    #' @description See [AnnData-usage]
    var_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "var", dim = "cols")
    },
    #' @description See [AnnData-usage]
    layers_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "layers")
    },
    #' @description See [AnnData-usage]
    obsm_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "obsm")
    },
    #' @description See [AnnData-usage]
    varm_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "varm")
    },
    #' @description See [AnnData-usage]
    obsp_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "obsp")
    },
    #' @description See [AnnData-usage]
    varp_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "varp")
    },
    #' @description See [AnnData-usage]
    uns_keys = function() {
      read_h5ad_element_keys(private$.hdf5_file, "uns")
    },

    #' @description See the `n_obs` field in [AnnData-usage]
    n_obs = function() {
      length(self$obs_names)
    },

    #' @description See the `n_vars` field in [AnnData-usage]
    n_vars = function() {
      length(self$var_names)
    },

    #' @description Open the HDF5 file handle
    #'
    #' @param readonly Whether to open in read-only mode
    open = function(readonly = FALSE) {
      if (private$.mode == "r" && !readonly) {
        cli_abort(
          paste(
            "This {.cls HDF5AnnData} is opened in read-only mode,",
            "cannot open the associated HDF5 file in read/write mode"
          )
        )
      }

      private$.hdf5_file$open(readonly = readonly)
    },

    #' @description Close the HDF5 file handle
    close = function() {
      private$.hdf5_file$close()
    },

    #' @description Convert to an [`InMemoryAnnData`]. Any backed
    #'   (`DelayedArray`) slots are materialized into ordinary in-memory
    #'   matrices: an in-memory object should not stay tied to an on-disk file.
    as_InMemoryAnnData = function() {
      if (isTRUE(private$.backed)) {
        prev <- private$.backed
        private$.backed <- FALSE
        on.exit(private$.backed <- prev, add = TRUE)
        # Hold the file open for the whole multi-slot read (single handle).
        private$.hdf5_file$open_and_defer_close(readonly = TRUE)
      }
      super$as_InMemoryAnnData()
    }
  )
)

#' Convert an `AnnData` to an `HDF5AnnData`
#'
#' Convert another `AnnData` object to an [`HDF5AnnData`] object
#'
#' @param adata An `AnnData` object to be converted to [`HDF5AnnData`]
#' @param file The file name (character) of the `.h5ad` file
#' @param compression The compression algorithm to use when writing the
#'   HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
#'   `"none"`.
#' @param chunk_size The target chunk size in bytes to use when writing
#'   HDF5 datasets. When `"auto"` (default), the chunk size is determined
#'   automatically using an algorithm that mimics h5py's auto-chunking
#'   behaviour. Set to `NULL` to disable chunking (contiguous storage,
#'   the rhdf5 default), or a number to use a specific target size in bytes.
#' @param mode The mode to open the HDF5 file:
#'
#'   * `a` creates a new file or opens an existing one for read/write
#'   * `r` opens an existing file for reading
#'   * `r+` opens an existing file for read/write
#'   * `w` creates a file, truncating any existing ones
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists
#' @param backed Whether the object is disk backed and returns
#'   [DelayedArray::DelayedArray] object for matrix data. Can only be `TRUE`
#'   when `mode == "r"`.
#'
#' @return An [`HDF5AnnData`] object with the same data as the input `AnnData`
#'   object.
#' @keywords internal
#'
#' @family object converters
#'
# nolint start: object_name_linter
as_HDF5AnnData <- function(
  # nolint end: object_name_linter
  adata,
  file,
  compression = c("none", "gzip", "lzf"),
  chunk_size = "auto",
  mode = c("w-", "r", "r+", "a", "w", "x"),
  backed = FALSE
) {
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

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
    shape = adata$shape(),
    mode = mode,
    compression = compression,
    chunk_size = chunk_size,
    backed = backed
  )
}
