#' @title HDF5File
#'
#' @description
#' A container class for HDF5 file handles. It stores the path to the HDF5 file
#' as well as the handle from [rhdf5::H5Fopen()] if the file is currently
#' opened, simplifying managing the file handle.
#'
#' @return An `HDF5File` object
#' @noRd
HDF5File <- R6::R6Class(
  "HDF5File", # nolint
  cloneable = FALSE,
  public = list(
    path = NULL,
    handle = NULL,
    is_open = NULL,
    is_readonly = FALSE,
    #' @description
    #' `HDF5File` constructor
    #'
    #' @param file The file name (character) of the `.hdf5` file. The file must
    #' already exist and be a valid HDF5 file.
    initialize = function(file) {
      if (!file.exists(file)) {
        cli_abort("File {.file {file}} does not exist")
      }

      if (!rhdf5::H5Fis_hdf5(file)) {
        cli_abort("File {.file {file}} is not a valid HDF5 file")
      }

      self$path <- file
      self$is_open <- FALSE
    },
    #' @description Format the `HDF5File` object
    #'
    #' @param ... Optional arguments to format method
    format = function(...) {
      status <- if (self$is_open) cli::col_green("OPEN") else cli::col_red("CLOSED")
      readonly <- if (!is.null(self$is_readonly) && self$is_readonly) {
        paste0("(", cli::col_yellow("read-only"), ")")
      } else {
        ""
      }
      cli::cli_format_method({
        cli::cli_text("{.cls HDF5File} {.file {self$path}} ({status}) {readonly}")
      })
    },
    #' @description Open the HDF5 file handle if needed
    #'
    #' @param readonly Whether to open in read-only mode.
    #'
    #' @return `TRUE` if the handle was successfully opened (invisibly)
    open = function(readonly = FALSE) {
      if (
        self$is_open &&
          private$.handle_is_valid() &&
          (self$is_readonly == readonly)
      ) {
        return(invisible(TRUE))
      }

      if (self$is_readonly != readonly) {
        self$close()
      }

      flags <- if (readonly) {
        "H5F_ACC_RDONLY"
      } else {
        "H5F_ACC_RDWR"
      }

      self$handle <- rhdf5::H5Fopen(self$path, flags = flags)
      self$is_open <- TRUE
      self$is_readonly <- readonly

      invisible(TRUE)
    },
    #' @description Close the HDF5 file handle if needed
    #'
    #' @return `TRUE` if the handle was successfully closed (invisibly)
    close = function(readonly = FALSE) {
      if (!self$is_open && !private$.handle_is_valid()) {
        return(invisible(TRUE))
      }

      if (private$.handle_is_valid()) {
        rhdf5::H5Fclose(self$handle)
      }
      self$is_open <- FALSE
      self$handle <- NULL

      invisible(TRUE)
    },
    #' @description Open the HDF5 file handle and automatically close it at the
    #' end of the calling function. If the file is already open, it is assumed
    #' the caller will close it and the deferred close call is not added.
    #'
    #' @param readonly Whether to open in read-only mode.
    #'
    #' @return `TRUE` if the handle was successfully opened (invisibly)
    open_and_defer_close = function(readonly = FALSE) {
      if (!self$is_open) {
        withr::defer_parent(self$close())
      }

      self$open(readonly = readonly)
    },
    #' @description Close the HDF5 file handle and automatically open it at the
    #' end of the calling function. If the file is already closed, it is assumed
    #' the caller will open it and the deferred open call is not added.
    #'
    #' @return `TRUE` if the handle was successfully closed (invisibly)
    close_and_defer_open = function() {
      if (self$is_open) {
        withr::defer_parent(self$open(readonly = self$is_readonly))
      }

      self$close()
    }
  ),
  private = list(
    finalize = function() {
      self$close()
    },
    .handle_is_valid = function() {
      if (is.null(self$handle)) {
        return(FALSE)
      }

      rhdf5::H5Iis_valid(self$handle)
    }
  )
)
