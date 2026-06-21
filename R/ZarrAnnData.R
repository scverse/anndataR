#' @title ZarrAnnData
#'
#' @description
#' Implementation of a Zarr-backed `AnnData` object. This class provides an
#' interface to a Zarr file and minimal data is stored in memory until it is
#' requested by the user. It is primarily designed as an intermediate object
#' when reading/writing Zarr files but can be useful for accessing parts of
#' large files.
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @return A `ZarrAnnData` object
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
ZarrAnnData <- R6::R6Class(
  "ZarrAnnData", # nolint
  inherit = AbstractAnnData,
  cloneable = FALSE,
  private = list(
    .zarrobj = NULL,
    .compression = NULL,
    .readonly = NULL,

    .check_file_valid = function() {
      if (!zarr_path_exists(private$.zarrobj, "/")) {
        cli_abort(
          "The Zarr path does not exist or is not a valid Zarr store"
        )
      }
    },

    .check_writeable = function() {
      if (isTRUE(private$.readonly)) {
        cli_abort(
          "Cannot write to a Zarr store opened in read-only mode.",
          call = rlang::caller_env()
        )
      }
    }
  ),
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_X, status=done
        read_zarr_element(private$.zarrobj, "X") |>
          private$.add_matrix_dimnames("X")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_X, status=done
        private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = self$obs_names,
          expected_colnames = self$var_names
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "X",
            private$.compression
          )
      }
    },
    #' @field layers See [AnnData-usage]
    layers = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_layers, status=done
        read_zarr_element(private$.zarrobj, "layers") |>
          private$.add_mapping_dimnames("layers")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_layers, status=done
        private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = self$obs_names,
          expected_colnames = self$var_names
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "layers",
            private$.compression
          )
      }
    },
    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_obsm, status=done
        read_zarr_element(private$.zarrobj, "obsm") |>
          private$.add_mapping_dimnames("obsm")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_obsm, status=done
        private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = self$obs_names,
          strip_rownames = TRUE,
          strip_colnames = FALSE,
          warn_colnames = TRUE
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "obsm",
            private$.compression
          )
      }
    },
    #' @field varm See [AnnData-usage]
    varm = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_varm, status=done
        read_zarr_element(private$.zarrobj, "varm") |>
          private$.add_mapping_dimnames("varm")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_varm, status=done
        private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = self$var_names,
          strip_rownames = TRUE,
          strip_colnames = FALSE,
          warn_colnames = TRUE
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "varm",
            private$.compression
          )
      }
    },
    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_obsp, status=done
        read_zarr_element(private$.zarrobj, "obsp") |>
          private$.add_mapping_dimnames("obsp")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_obsp, status=done
        private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = self$obs_names,
          expected_colnames = self$obs_names
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "obsp",
            private$.compression
          )
      }
    },
    #' @field varp See [AnnData-usage]
    varp = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_varp, status=done
        read_zarr_element(private$.zarrobj, "varp") |>
          private$.add_mapping_dimnames("varp")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_varp, status=done
        private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = self$var_names,
          expected_colnames = self$var_names
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "varp",
            private$.compression
          )
      }
    },
    #' @field obs See [AnnData-usage]
    obs = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_obs, status=done
        read_zarr_element(private$.zarrobj, "obs")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_obs, status=done
        private$.validate_obsvar_dataframe(value, "obs") |>
          write_zarr_element(
            private$.zarrobj,
            "obs",
            private$.compression
          )
      }
    },
    #' @field var See [AnnData-usage]
    var = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_var, status=done
        read_zarr_element(private$.zarrobj, "var")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_var, status=done
        private$.validate_obsvar_dataframe(value, "var") |>
          write_zarr_element(
            private$.zarrobj,
            "var",
            private$.compression
          )
      }
    },
    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_obs_names, status=done
        read_zarr_element_keys(private$.zarrobj, "obs", dim = "rows")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_obs_names, status=done
        rownames(self$obs) <- value
      }
    },
    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_var_names, status=done
        read_zarr_element_keys(private$.zarrobj, "var", dim = "rows")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_var_names, status=done
        rownames(self$var) <- value
      }
    },
    #' @field uns See [AnnData-usage]
    uns = function(value) {
      private$.check_file_valid()

      if (missing(value)) {
        # trackstatus: class=ZarrAnnData, feature=get_uns, status=done
        read_zarr_element(private$.zarrobj, "uns")
      } else {
        private$.check_writeable()
        # trackstatus: class=ZarrAnnData, feature=set_uns, status=done
        private$.validate_named_list(
          value,
          "uns",
          warn_matrix_dimnames = TRUE
        ) |>
          write_zarr_element(
            private$.zarrobj,
            "uns",
            private$.compression
          )
      }
    }
  ),
  public = list(
    #' @description
    #' `ZarrAnnData` constructor
    #'
    #' @param file The file name (character) of the `.zarr` file. If this file
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
    #' @param mode The mode to open the Zarr file. See [as_ZarrAnnData()] for
    #'   details
    #' @param compression The compression algorithm to use. See
    #'   [as_ZarrAnnData()] for details
    #'
    #' @details
    #' The constructor creates a new Zarr `AnnData` interface object. This can
    #' either be used to either connect to an existing `.zarr` file or to
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
      compression = c(
        "none",
        "gzip",
        "blosc",
        "zstd",
        "lzma",
        "bz2",
        "zlib",
        "lz4"
      )
    ) {
      check_requires("ZarrAnnData", "Rarr", where = "Bioc")

      compression <- match.arg(compression)
      mode <- match.arg(mode)

      private$.compression <- compression

      is_readonly <- FALSE

      if (is.character(file)) {
        if (mode == "a") {
          if (dir.exists(file)) {
            mode <- "r+"
          } else {
            mode <- "w-"
          }
        }

        if (!dir.exists(file) && mode %in% c("r", "r+")) {
          cli_abort(
            paste(
              "File {.file {file}} does not exist but mode is set to {.val {mode}}.",
              "If you want to create a new file, use a different mode (e.g. 'w-').",
              "See {.help read_zarr} or {.help write_zarr} for more information."
            ),
            call = rlang::caller_env()
          )
        }

        if (dir.exists(file) && mode %in% c("w-", "x")) {
          cli_abort(
            paste(
              "File {.file {file}} already exists but mode is set to {.val {mode}}.",
              "If you want to overwrite the file, use a different mode (e.g. 'w').",
              "See {.help read_zarr} or {.help write_zarr} for more information."
            ),
            call = rlang::caller_env()
          )
        }

        if (mode %in% c("w", "w-", "x")) {
          create_zarr(file)
        } else if (mode == "r") {
          is_readonly <- TRUE
        }
      } else {
        cli_abort(
          paste(
            "{.arg file} must be a {.cls character}"
          )
        )
      }

      if (!zarr_path_exists(file, "/")) {
        cli_abort(
          paste(
            "{.arg file} must be a valid zarr store/file"
          )
        )
      }

      is_empty <- is_zarr_empty(file)

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
          write_empty_zarr(file, obs, var, compression)
        }
      }

      # File is supposed to exist by now. Check if it is a valid Zarr file
      attrs <- Rarr::read_zarr_attributes(file)
      if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
        cli_abort(c(
          "File {.file {file}} is not a valid AnnData-Zarr file."
        ))
      }

      # Set the file path
      private$.zarrobj <- file
      private$.readonly <- is_readonly

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
              ") to a Zarr file opened in read-only mode."
            )
          )
        }
      } else {
        for (slot in .anndata_slots) {
          value <- get(slot)
          if (!is.null(value)) {
            self[[slot]] <- value
          }
        }
      }

      self
    },

    #' @description See the `n_obs` field in [AnnData-usage]
    n_obs = function() {
      length(self$obs_names)
    },

    #' @description See the `n_vars` field in [AnnData-usage]
    n_vars = function() {
      length(self$var_names)
    },

    #' @description See [AnnData-usage]
    obs_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "obs", dim = "cols")
    },
    #' @description See [AnnData-usage]
    var_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "var", dim = "cols")
    },
    #' @description See [AnnData-usage]
    layers_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "layers")
    },
    #' @description See [AnnData-usage]
    obsm_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "obsm")
    },
    #' @description See [AnnData-usage]
    varm_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "varm")
    },
    #' @description See [AnnData-usage]
    obsp_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "obsp")
    },
    #' @description See [AnnData-usage]
    varp_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "varp")
    },
    #' @description See [AnnData-usage]
    uns_keys = function() {
      read_zarr_element_keys(private$.zarrobj, "uns")
    }
  )
)

#' Convert an `AnnData` to an `ZarrAnnData`
#'
#' Convert another `AnnData` object to an [`ZarrAnnData`] object
#'
#' @param adata An `AnnData` object to be converted to [`ZarrAnnData`]
#' @param file The file name (character) of the `.zarr` file
#' @param compression The compression algorithm to use when writing the
#'   Zarr file. Can be one of `"none"`, `"gzip"`, `"blosc"`, `"zstd"`,
#'   `"lzma"`, `"bz2"`, `"zlib"` or `"lz4"`. Defaults to `"none"`.
#' @param mode The mode to open the Zarr file:
#'
#'   * `a` creates a new file or opens an existing one for read/write
#'   * `r` opens an existing file for reading
#'   * `r+` opens an existing file for read/write
#'   * `w` creates a file, truncating any existing ones
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists
#'
#' @return A [`ZarrAnnData`] object with the same data as the input `AnnData`
#'   object.
#' @keywords internal
#'
#' @family object converters
#'
# nolint start: object_name_linter
as_ZarrAnnData <- function(
  # nolint end: object_name_linter
  adata,
  file,
  compression = c(
    "none",
    "gzip",
    "blosc",
    "zstd",
    "lzma",
    "bz2",
    "zlib",
    "lz4"
  ),
  mode = c("w-", "r", "r+", "a", "w", "x")
) {
  if (!(inherits(adata, "AbstractAnnData"))) {
    cli_abort(
      "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
    )
  }

  mode <- match.arg(mode)
  ZarrAnnData$new(
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
    compression = compression
  )
}
