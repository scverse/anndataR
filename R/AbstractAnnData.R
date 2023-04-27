#' @title AbstractAnnData
#'
#' @description
#' Abstract class representing an anndata object. Defines the interface.
#' @importFrom R6 R6Class
AbstractAnnData <- R6::R6Class("AbstractAnnData",
  public = list(


    #' @description Abstract class, representing an AnnData object.
    initialize = function() {
    },


    #' Write h5ad
    #' @param filename Filename of data file. Defaults to backing file.
    #' @param compression compression
    #' @param compression_opts compression options
    #' @param as_dense Sparse arrays in AnnData obect to write as dense.
    write_h5ad = function(filename = NULL, compression = NULL, compression_opts = NULL, as_dense = NULL){
    },

    X = function(){
      stop("This is an abstract class. What you're trying to do makes no sense.")
    },
    
    obs = function(){
      stop("This is an abstract class. What you're trying to do makes no sense.")
    },
    
    var = function(){
      stop("This is an abstract class. What you're trying to do makes no sense.")
    }


  ))
