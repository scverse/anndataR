#' S3 Methods for AbstractAnnData Objects

#' @rdname AnnDataSlotList-s3methods
#' @method names AnnDataSlotList
#' @export
names.AnnDataSlotList <- function(x) {
  x$get_keys_fn()
}

#' @rdname AnnDataSlotList-s3methods
#' @method names<- AnnDataSlotList
#' @export
`names<-.AnnDataSlotList` <- function(x, value) {
  x$set_keys_fn(value)
  x
}

#' @rdname AnnDataSlotList-s3methods
#' @method [ AnnDataSlotList
#' @export
`[.AnnDataSlotList` <- function(x, i, ...) {
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

#' @rdname AnnDataSlotList-s3methods
#' @method [[ AnnDataSlotList
#' @export
`[[.AnnDataSlotList` <- function(x, i, ...) {
  x$get_value_fn(i)
}

#' @rdname AnnDataSlotList-s3methods
#' @method [[<- AnnDataSlotList
#' @export
`[[<-.AnnDataSlotList` <- function(x, i, value) {
  x$set_value_fn(i, value)
  x
}
