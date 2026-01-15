S3AnnDataSlotList <- function(
    get_keys_fn,
    set_keys_fn,
    get_value_fn,
    set_value_fn,
    set_values_fn,
    get_rownames_fn = NULL,
    get_colnames_fn = NULL
  ) {
  structure(
    list(
      get_keys_fn = get_keys_fn,
      set_keys_fn = set_keys_fn,
      get_value_fn = get_value_fn,
      set_value_fn = set_value_fn,
      set_values_fn = set_values_fn,
      get_rownames_fn = get_rownames_fn,
      get_colnames_fn = get_colnames_fn
    ),
    class = c("S3AnnDataSlotList", "list")
  )
}

interact <- function(x, val, ...) {
  if(missing(val)){
    
  }
}



#' @method [ S3AnnDataSlotList
#' @export
`[.S3AnnDataSlotList` <- function(x, i, ...) {
  if (missing(i)) {
    # return all values when no index is provided
    keys <- x$get_keys_fn()
    result <- lapply(keys, x$get_value_fn)
    names(result) <- keys
    result
  } else {
    res <- list(x$get_value_fn(i))
    names(res) <- i
    res
  }
}

get_keys_fn.S3AnnDataSlotList <- function(x) {
  x$get_keys_fn()
}

# #' @method [[ S3AnnDataSlotList
# #' @export
# `[[.S3AnnDataSlotList` <- function(x, i, ...) {
#   x$get_value_fn(i)
# }

# #' @method [[<- S3AnnDataSlotList
# #' @export
# `[[<-.S3AnnDataSlotList` <- function(x, i, value) {
#   x$set_value_fn(i, value)
#   x
# }

# #' @method $ S3AnnDataSlotList
# #' @export
# `$.S3AnnDataSlotList` <- function(x, name) {
#   x$get_value_fn(name)
# }
