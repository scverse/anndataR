#' @title HDF5AnnData
#'
#' @description
#' Implementation of an in memory AnnData object.
HDF5AnnData <- R6::R6Class("HDF5AnnData",
  inherit = AbstractAnnData,
  private = list(
    .h5obj = NULL
  ),
  active = list(
    #' @field X The X slot
    X = function(value) {
      if (missing(value)) {
        # return X
      } else {
        # set X
      }
    },
    #' @field obs The obs slot
    obs = function(value) {
      if (missing(value)) {
        # return obs
      } else {
        # set obs
      }
    },
    #' @field var The var slot
    var = function(value) {
      if (missing(value)) {
        # return var
      } else {
        # set var
      }
    }
  ),
  public = list(
    #' @description HDF5AnnData constructor
    #'
    #' @param h5obj The rhdf5 object
    initialize = function(h5obj) {
      private$.h5obj <- h5obj
    }
  )
)

.create_empty_HDF5AnnData <- function(path, n_obs, n_vars) { # nolint
  # try to replicate this:
  # nolint start
  # adata <- anndata::AnnData(shape = c(10, 20))
  # adata$write_h5ad("test.h5ad")
  # h5obj <- rhdf5::H5Fopen("test.h5ad")
  # rhdf5::h5ls(h5obj)
  requireNamespace("rhdf5")
  empty <- rhdf5::H5Fcreate(path)

  # todo: uncomment once hdf5 functions have been merged
  # empty <- write_h5ad_mapping(empty, "layers", list())
  # empty <- write_h5ad_data_frame(empty, "obs", data.frame(row.names = as.character(seq_len(n_obs))))
  # empty <- write_h5ad_data_frame(empty, "var", data.frame(row.names = as.character(seq_len(n_vars))))
  # nolint end

  # todo: add obsm, obsp, uns, varm, varp

  HDF5AnnData$new(empty)
}

HDF5AnnData$create_empty <- .create_empty_HDF5AnnData