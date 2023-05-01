#' ReticulateAnnData class
#'
#' An R6 class that provides an interface to interact with AnnData objects in Python using the reticulate package.
ReticulateAnnData <- R6::R6Class("ReticulateAnnData",
  inherit = AbstractAnnData,
  private = list(
    .obj = NULL,

    py_to_r_index = function(x) {
      requireNamespace("reticulate")
      python_builtins <- reticulate::import_builtins()
      out <- python_builtins$list(x)
      attr(out, "name") <- py_to_r_ifneedbe(x$name)
      out
    }
  ),
  active = list(
    #' @field X A matrix representing the X slot.
    X = function(value) {
      if (missing(value)) {
        py_to_r_ifneedbe(private$.obj$X)
      } else {
        private$.obj$X <- value
        self
      }
    },
    #' @field obs A data.frame representing the obs slot.
    obs = function(value) {
      if (missing(value)) {
        py_to_r_ifneedbe(private$.obj$obs)
      } else {
        private$.obj$obs <- value
        self
      }
    },
    #' @field var A data.frame representing the var slot.
    var = function(value) {
      if (missing(value)) {
        py_to_r_ifneedbe(private$.obj$var)
      } else {
        private$.obj$var <- value
        self
      }
    },
    #' @field obs_names A vector representing the obs_names slot.
    obs_names = function(value) {
      if (missing(value)) {
        private$py_to_r_index(private$.obj$obs_names)
      } else {
        private$.obj$obs_names <- value
        self
      }
    },
    #' @field var_names A vector representing the var_names slot.
    var_names = function(value) {
      if (missing(value)) {
        private$py_to_r_index(private$.obj$var_names)
      } else {
        private$.obj$var_names <- value
        self
      }
    },
    #' @field layers A list representing the layers slot.
    layers = function(value) {
      if (missing(value)) {
        py_to_r_ifneedbe(private$.obj$layers)
      } else {
        private$.obj$layers <- value
        self
      }
    }
  ),
  public = list(
    #' @description ReticulateAnnData constructor
    #'
    #' @param obj The reticulate AnnData object
    initialize = function(obj) {
      private$.obj <- obj
    },

    # copied from https://github.com/dynverse/anndata/blob/main/R/class_anndata.R#LL562C1-L601C7
    #' @description Write .h5ad-formatted hdf5 file.
    #'
    #' Generally, if you have sparse data that are stored as a dense matrix, you can
    #' dramatically improve performance and reduce disk space by converting to a csr_matrix:
    #'
    #' @param filename Filename of data file. Defaults to backing file.
    #' @param compression See the h5py
    #'   [filter pipeline](http://docs.h5py.org/en/latest/high/dataset.html#dataset-compression).
    #'   Options are `"gzip"`, `"lzf"` or `NULL`.
    #' @param compression_opts See the h5py
    #'   [filter pipeline](http://docs.h5py.org/en/latest/high/dataset.html#dataset-compression).
    #' @param as_dense Sparse in AnnData object to write as dense. Currently only supports `"X"` and `"raw/X"`.
    write_h5ad = function(filename, compression = NULL, compression_opts = NULL, as_dense = list()) {
      filename <- normalizePath(filename, mustWork = FALSE)
      invisible(py_to_r_ifneedbe(private$.obj$write_h5ad(
        filename = filename,
        compression = compression,
        compression_opts = compression_opts,
        as_dense = as_dense
      )))
    }
  )
)

py_to_r_ifneedbe <- function(x) {
  requireNamespace("reticulate")
  if (inherits(x, "python.builtin.object")) {
    reticulate::py_to_r(x)
  } else {
    x
  }
}

# copied from https://github.com/dynverse/anndata/blob/main/R/class_anndata.R#L127
to_reticulate <- function(adata) {
  requireNamespace("reticulate")
  python_anndata <- reticulate::import("anndata", convert = FALSE)

  ad <- python_anndata$AnnData(
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    layers = adata$layers
  )
  ad$obs_names <- adata$obs_names
  ad$var_names <- adata$var_names

  ReticulateAnnData$new(ad)
}

# copied from https://github.com/dynverse/anndata/blob/main/R/read_h5ad.R
ReticulateAnnData$read_h5ad <- function(
  filename,
  backed = NULL
) {
  requireNamespace("reticulate")
  python_anndata <- reticulate::import("anndata", convert = FALSE)
  filename <- normalizePath(filename, mustWork = FALSE)
  ReticulateAnnData$new(python_anndata$read_h5ad(
    filename = filename,
    backed = backed
  ))
}