#' @title Lazy Named List
#'
#' @description A lazy named list that loads elements on-demand to avoid
#' materializing large objects unnecessarily. Used internally for efficient
#' access to layers, obsm, varm, obsp, varp, uns slots.
#'
#' @keywords internal
#' @export
AnnDataSlotList <- function(
  get_keys_fn,
  set_keys_fn,
  get_value_fn,
  set_value_fn,
  set_values_fn,
  get_rownames_fn = NULL,
  get_colnames_fn = NULL
) {
  self <- new.env(parent = emptyenv())
  self$get_keys_fn     <- get_keys_fn
  self$set_keys_fn     <- set_keys_fn
  self$get_value_fn    <- get_value_fn
  self$set_value_fn    <- set_value_fn
  self$set_values_fn   <- set_values_fn
  self$get_rownames_fn <- get_rownames_fn
  self$get_colnames_fn <- get_colnames_fn
  class(self) <- "AnnDataSlotList"
  self
}

#' S3 Methods for AnnDataSlotList Objects
#' @rdname AnnDataSlotList-s3methods
#' @method names AnnDataSlotList
#' @export
names.AnnDataSlotList <- function(x) {
  .subset2(x, "get_keys_fn")()
}

#' @rdname AnnDataSlotList-s3methods
#' @method names<- AnnDataSlotList
#' @export
`names<-.AnnDataSlotList` <- function(x, value) {
  .subset2(x, "set_keys_fn")(value)
  x
}

#' @rdname AnnDataSlotList-s3methods
#' @method [ AnnDataSlotList
#' @export
`[.AnnDataSlotList` <- function(x, i, ...) {
  get_value_fn <- .subset2(x, "get_value_fn")
  if (missing(i)) {
    keys <- .subset2(x, "get_keys_fn")()
    result <- lapply(keys, get_value_fn)
    names(result) <- keys
    result
  } else {
    res <- list(get_value_fn(i))
    names(res) <- i
    res
  }
}

#' @rdname AnnDataSlotList-s3methods
#' @method [[ AnnDataSlotList
#' @export
`[[.AnnDataSlotList` <- function(x, i, ...) {
  .subset2(x, "get_value_fn")(i)
}

#' @rdname AnnDataSlotList-s3methods
#' @method [[<- AnnDataSlotList
#' @export
`[[<-.AnnDataSlotList` <- function(x, i, value) {
  .subset2(x, "set_value_fn")(i, value)
  x
}

#' @rdname AnnDataSlotList-s3methods
#' @method $ AnnDataSlotList
#' @export
`$.AnnDataSlotList` <- function(x, name) {
  internal_fields <- ls(x)
  if (name %in% internal_fields) {
    return(.subset2(x, name))
  }
  keys <- .subset2(x, "get_keys_fn")()
  if (!(name %in% keys)) {
    stop(sprintf(
      "Key '%s' not found in AnnDataSlotList (available keys: %s)",
      name, paste(keys, collapse = ", ")
    ), call. = FALSE)
  }
  .subset2(x, "get_value_fn")(name)
}

#' @rdname AnnDataSlotList-s3methods
#' @method $<- AnnDataSlotList
#' @export
`$<-.AnnDataSlotList` <- function(x, name, value) {
  internal_fields <- ls(x)
  if (name %in% internal_fields) {
    assign(name, value, envir = x)
    return(x)
  }
  .subset2(x, "set_value_fn")(name, value)
  x
}

#' @rdname AnnDataSlotList-s3methods
#' @method print AnnDataSlotList
#' @export
print.AnnDataSlotList <- function(x, ...) {
  keys <- .subset2(x, "get_keys_fn")()
  cat("AnnDataSlotList with", length(keys), "element(s)")
  if (length(keys) > 0) {
    cat(":", paste(keys, collapse = ", "))
  }
  cat("\n")
  invisible(x)
}
