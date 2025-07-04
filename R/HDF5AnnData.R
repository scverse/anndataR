#' @title HDF5AnnData
#'
#' @description
#' Implementation of an HDF5-backed `AnnData` object.
#'
#' See [AnnData-usage] for details on creating and using `AnnData` objects.
#'
#' @seealso [AnnData-usage] for details on creating and using `AnnData` objects
#'
#' @family AnnData classes
HDF5AnnData <- R6::R6Class(
  "HDF5AnnData", # nolint
  inherit = AbstractAnnData,
  cloneable = FALSE,
  private = list(
    .h5obj = NULL,
    .close_on_finalize = FALSE,
    .compression = NULL,
    .rhdf5 = FALSE,

    #' @description Close the HDF5 file when the object is garbage collected
    finalize = function() {
      if (private$.close_on_finalize) {
        self$close()
      }
      invisible(self)
    }
  ),
  active = list(
    #' @field X See [AnnData-usage]
    X = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_X, status=done
        read_h5ad_element(private$.h5obj, "X", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_X, status=done
        value <- private$.validate_aligned_array(
          value,
          "X",
          shape = c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "X",
          private$.compression,
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field layers See [AnnData-usage]
    layers = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_layers, status=done
        read_h5ad_element(private$.h5obj, "layers", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_layers, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "layers",
          c(self$n_obs(), self$n_vars()),
          expected_rownames = rownames(self),
          expected_colnames = colnames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "layers",
          private$.compression,
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field obsm See [AnnData-usage]
    obsm = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsm, status=done
        read_h5ad_element(private$.h5obj, "obsm", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obsm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsm",
          c(self$n_obs()),
          expected_rownames = rownames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "obsm",
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field varm See [AnnData-usage]
    varm = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varm, status=done
        read_h5ad_element(private$.h5obj, "varm", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_varm, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varm",
          c(self$n_vars()),
          expected_rownames = colnames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "varm",
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field obsp See [AnnData-usage]
    obsp = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obsp, status=done
        read_h5ad_element(private$.h5obj, "obsp", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obsp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "obsp",
          c(self$n_obs(), self$n_obs()),
          expected_rownames = rownames(self),
          expected_colnames = rownames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "obsp",
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field varp See [AnnData-usage]
    varp = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_varp, status=done
        read_h5ad_element(private$.h5obj, "varp", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_varp, status=done
        value <- private$.validate_aligned_mapping(
          value,
          "varp",
          c(self$n_vars(), self$n_vars()),
          expected_rownames = colnames(self),
          expected_colnames = colnames(self)
        )
        write_h5ad_element(
          value,
          private$.h5obj,
          "varp",
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field obs See [AnnData-usage]
    obs = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs, status=done
        read_h5ad_element(private$.h5obj, "obs", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs, status=done
        value <- private$.validate_obsvar_dataframe(value, "obs")
        write_h5ad_element(
          value,
          private$.h5obj,
          "obs",
          private$.compression,
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field var See [AnnData-usage]
    var = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var, status=done
        read_h5ad_element(private$.h5obj, "var", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var, status=done
        value <- private$.validate_obsvar_dataframe(value, "var")
        write_h5ad_element(
          value,
          private$.h5obj,
          "var",
          rhdf5 = private$.rhdf5
        )
      }
    },
    #' @field obs_names See [AnnData-usage]
    obs_names = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_obs_names, status=done
        rownames(self$obs)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_obs_names, status=done
        rownames(self$obs) <- value
      }
    },
    #' @field var_names See [AnnData-usage]
    var_names = function(value) {
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_var_names, status=done
        rownames(self$var)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_var_names, status=done
        rownames(self$var) <- value
      }
    },
    #' @field uns See [AnnData-usage]
    uns = function(value) {
      if (!private$.rhdf5 && !private$.h5obj$is_valid) {
        cli_abort("HDF5 file is closed")
      }
      if (missing(value)) {
        # trackstatus: class=HDF5AnnData, feature=get_uns, status=done
        read_h5ad_element(private$.h5obj, "uns", rhdf5 = private$.rhdf5)
      } else {
        # trackstatus: class=HDF5AnnData, feature=set_uns, status=done
        value <- private$.validate_named_list(value, "uns")
        write_h5ad_element(value, private$.h5obj, "uns", rhdf5 = private$.rhdf5)
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
    #' @param compression The compression algorithm to use. See
    #'   [as_HDF5AnnData()] for details
    #' @param mode The mode to open the HDF5 file. See [as_HDF5AnnData()] for
    #'   details
    #' @param rhdf5 Whether to use the `rhdf5` package for HDF5 operations
    #'
    #' @details
    #' The constructor creates a new HDF5 `AnnData` interface object. This can
    #' either be used to either connect to an existing `.h5ad` file or to
    #' create a new one. To create a new file both `obs_names` and `var_names`
    #' must be specified. In both cases, any additional slots provided will be
    #' set on the created object. This will cause data to be overwritten if the
    #' file already exists.
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
      mode = c("r", "r+", "a", "w", "w-", "x"),
      compression = c("none", "gzip", "lzf"),
      rhdf5 = FALSE
    ) {
      if (rhdf5) {
        self$rhdf5_initialize(
          file = file,
          X = X,
          obs = obs,
          var = var,
          layers = layers,
          obsm = obsm,
          varm = varm,
          obsp = obsp,
          varp = varp,
          uns = uns,
          shape = shape,
          compression = compression
        )
        return(self)
      }

      check_requires("HDF5AnnData", "hdf5r")

      # check arguments
      compression <- match.arg(compression)
      mode <- match.arg(mode)

      # store compression for later use
      private$.compression <- compression

      # check whether file needs to be closed on finalize
      private$.close_on_finalize <- is.character(file)

      # check whether file already exists
      if (is.character(file)) {
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

        if (file.exists(file) && mode %in% c("w-", "x")) {
          cli_abort(
            paste(
              "File {.file {file}} does not exist but mode is set to {.val {mode}}.",
              "If you want to overwrite the file, use a different mode (e.g. 'w').",
              "See {.help read_h5ad} or {.help write_h5ad} for more information."
            ),
            call = rlang::caller_env()
          )
        }

        file <- hdf5r::H5File$new(file, mode = mode)
      }

      # type check
      if (!inherits(file, "H5File")) {
        cli_abort(
          paste(
            "{.arg file} must be a {.cls character} or a {.cls hdf5r::H5File} object,",
            "but is a {.cls {class(file)}}"
          )
        )
      }

      # detect file mode
      is_readonly <- file$get_intent() == "H5F_ACC_RDONLY"
      is_empty <- length(hdf5r::list.groups(file)) == 0L &&
        length(hdf5r::list.datasets(file)) == 0L &&
        length(hdf5r::list.attributes(file)) == 0L &&
        length(hdf5r::list.objects(file)) == 0L

      if (!is_readonly) {
        if (!is_empty) {
          cli_warn(
            paste(
              "File is opened in read/write mode.",
              "Use with caution, as this can lead to data corruption."
            )
          )
        } else {
          # initialise AnnData
          shape <- get_shape(obs, var, X, shape)
          obs <- get_initial_obs(obs, X, shape)
          var <- get_initial_var(var, X, shape)
          write_empty_h5ad(file, obs, var, compression)
        }
      }

      # File is supposed to exist by now. Check if it is a valid H5AD file
      attrs <- hdf5r::h5attributes(file)
      if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
        cli_abort(c(
          "File {.file {file}} is not a valid H5AD file.",
          i = "Either the file is not an H5AD file or it was created with {.pkg anndata<0.8.0}."
        ))
      }

      # Set the file path
      private$.h5obj <- file

      if (is_readonly) {
        # if any of these variables are not NULL, throw an error
        are_null <- sapply(.anndata_slots, function(x) is.null(get(x)))
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
        for (slot in .anndata_slots) {
          value <- get(slot)
          if (!is.null(value)) {
            self[[slot]] <- value
          }
        }
      }

      self
    },

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
    #' @param compression The compression algorithm to use. See
    #'   [as_HDF5AnnData()] for details
    #'
    #' @details
    #' The constructor creates a new HDF5 `AnnData` interface object. This can
    #' either be used to either connect to an existing `.h5ad` file or to
    #' create a new one. If any additional slot arguments are set an existing
    #' file will be overwritten.
    rhdf5_initialize = function(
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
      compression = c("none", "gzip", "lzf")
    ) {
      check_requires("HDF5AnnData", "rhdf5")

      compression <- match.arg(compression)

      private$.compression <- compression

      private$.close_on_finalize <- is.character(file)

      write_mode <- any(
        !is.null(X),
        !is.null(obs),
        !is.null(var),
        !is.null(layers),
        !is.null(obsm),
        !is.null(varm),
        !is.null(obsp),
        !is.null(varp),
        !is.null(uns),
        !is.null(shape)
      )

      write_empty <- write_mode

      file_exists <- (is.character(file) && file.exists(file)) ||
        inherits(file, "H5IdComponent")

      if (write_mode && file_exists) {
        if (is.character(file)) {
          unlink(file)
          file <- rhdf5::H5Fcreate(file, native = FALSE)
        } else {
          rhdf5_hdf5_clear_file(file)
        }
      } else if (write_mode && !file_exists) {
        file <- rhdf5::H5Fcreate(file, native = FALSE)
      } else if (!write_mode && file_exists) {
        if (is.character(file)) {
          file <- rhdf5::H5Fopen(file, native = FALSE)
        }
      } else if (!write_mode && !file_exists) {
        file <- rhdf5::H5Fcreate(file, native = FALSE)
        write_empty <- TRUE
      }

      if (!inherits(file, "H5IdComponent")) {
        cli_abort(
          paste(
            "{.arg file} must be a {.cls character} or a {.cls rhdf5::H5IdComponent} object,",
            "but is a {.cls {class(file)}}"
          )
        )
      }

      if (write_empty) {
        shape <- get_shape(obs, var, X, shape)
        obs <- get_initial_obs(obs, X, shape)
        var <- get_initial_var(var, X, shape)
        rhdf5_write_empty_h5ad(file, obs, var, compression)
      }

      # File is supposed to exist by now. Check if it is a valid H5AD file
      attrs <- rhdf5::h5readAttributes(file, "/")
      if (!all(c("encoding-type", "encoding-version") %in% names(attrs))) {
        path <- rhdf5::H5Fget_name(file)
        cli_abort(c(
          "File {.file {path}} is not a valid H5AD file.",
          i = "Either the file is not an H5AD file or it was created with {.pkg anndata<0.8.0}."
        ))
      }

      # Set the file path
      private$.h5obj <- file
      private$.rhdf5 <- TRUE

      if (write_mode) {
        for (slot in .anndata_slots) {
          value <- get(slot)
          if (!is.null(value)) {
            self[[slot]] <- value
          }
        }
      }

      self
    },

    #' @description Close the HDF5 file
    close = function() {
      if (!private$.rhdf5 && private$.h5obj$is_valid) {
        private$.h5obj$close_all()
        cleanup_HDF5AnnData(file = private$.h5obj$file)
        gc()
      }

      if (private$.rhdf5) {
        h5obj <- private$.h5obj
        rhdf5::H5Fclose(private$.h5obj)
        rhdf5::H5garbage_collect()
        gc()
      }
    },

    #' @description See the `n_obs` field in [AnnData-usage]
    n_obs = function() {
      nrow(self$obs)
    },

    #' @description See the `n_vars` field in [AnnData-usage]
    n_vars = function() {
      nrow(self$var)
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
#' @param mode The mode to open the HDF5 file:
#'
#'   * `a` creates a new file or opens an existing one for read/write
#'   * `r` opens an existing file for reading
#'   * `r+` opens an existing file for read/write
#'   * `w` creates a file, truncating any existing ones
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists
#' @param rhdf5 Whether to use the `rhdf5` package for HDF5 operations
#'
#' @return An [`HDF5AnnData`] object with the same data as the input `AnnData`
#'   object.
#' @keywords internal
#'
#' @family object converters
#'
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5),
#' )
#' ad$as_HDF5AnnData("test.h5ad")
#' # remove file
#' unlink("test.h5ad")
# nolint start: object_name_linter
as_HDF5AnnData <- function(
  # nolint end: object_name_linter
  adata,
  file,
  compression = c("none", "gzip", "lzf"),
  mode = c("w-", "r", "r+", "a", "w", "x"),
  rhdf5 = FALSE
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
    compression = compression,
    shape = adata$shape(),
    mode = mode,
    rhdf5 = rhdf5
  )
}

# nolint start: object_name_linter
cleanup_HDF5AnnData <- function(...) {
  # nolint end: object_name_linter
  args <- list(...)

  if (
    !is.null(args$file) && is.character(args$file) && file.exists(args$file)
  ) {
    cli::cli_alert("Removing file: ", args$file)
    unlink(args$file)
  }
}
