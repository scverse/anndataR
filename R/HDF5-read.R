#' Read HDF5 encoding
#'
#' Read the encoding and version of an element in a HDF5 file
#'
#' @param file Path to a HDF5 file or an open HDF5 handle
#' @param name Name of the element within the HDF5 file
#'
#' @return A named list with the encoding and version
read_hdf5_encoding <- function(file, path) {
  tryCatch(
    error = function(cnd){
      print(paste0("An error occurred when reading ", name, " in file ", file, "."))
      conditionMessage(cnd)
    },
    h5readAttributes(test_file, "X")
  )
}

#' Read HDF5 dense array
#'
#' Read a dense array from an HDF5 file
#'
#' @param file Path to a HDF5 file or an open HDF5 handle
#' @param name Name of the element within the HDF5 file
#' @param version Encoding version of the dense array
#'
#' @return a Matrix/sparse matrix/DelayedArray???
read_hdf5_dense_array <- function(file, name, version = c("0.2.0")) {
  tryCatch(
    error = function(cnd){
      print(paste0("An error occurred when reading ", name, " in file."))
      conditionMessage(cnd)
    },
    {
      if(version == c("0.2.0")){
        h5read(file, name)
      }
    }
  )
}


