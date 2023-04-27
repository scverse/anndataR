.abstract_function <- function() {
  stop("This is an abstract method which should not be called directly.")
}

#' @title AbstractAnnData
#'
#' @description
#' Abstract class representing an anndata object. Defines the interface.
#' @importFrom R6 R6Class
AbstractAnnData <- R6::R6Class("AbstractAnnData",
  public = list(
    #' @description Write h5ad
    #' @param filename Filename of data file. Defaults to backing file.
    #' @param compression compression
    #' @param compression_opts compression options
    #' @param as_dense Sparse arrays in AnnData obect to write as dense.
    write_h5ad = function(filename = NULL, compression = NULL, compression_opts = NULL, as_dense = NULL){
      .abstract_function()
    },

    X = .abstract_function,
    obs = .abstract_function,
    var = .abstract_function
  )
)