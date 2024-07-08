
# copied from https://github.com/dynverse/anndata/blob/main/R/reticulate_conversions.R#L1/
py_to_r_ifneedbe <- function(x) {
  requireNamespace("reticulate")
  if (inherits(x, "python.builtin.object")) {
    reticulate::py_to_r(x)
  } else {
    x
  }
}

#' Convert between Python and R objects
#'
#' @param x A Python object.
#' @param name A name
#' @param value A value
#'
#' @return An \R object, as converted from the Python object.
#'
#' @name r-py-conversion
#' @export
# copied from https://github.com/dynverse/anndata/blob/main/R/reticulate_conversions.R#L61
py_to_r.pandas.core.indexes.base.Index <- function(x) { # nolint
  requireNamespace("reticulate")
  python_builtins <- reticulate::import_builtins()
  out <- python_builtins$list(x)
  attr(out, "name") <- py_to_r_ifneedbe(x$name)
  out
}

#' @name r-py-conversion
#' @export
# copied from https://github.com/dynverse/anndata/blob/main/R/reticulate_conversions.R#LL77C38-L90C2
py_to_r.collections.abc.Mapping <- function(x) { # nolint
  requireNamespace("reticulate")
  python_builtins <- reticulate::import_builtins()

  x_list <- python_builtins$dict(x)

  # convert members of x_list if need be
  for (i in seq_along(x_list)) {
    if (inherits(x_list[[i]], "python.builtin.object")) {
      x_list[[i]] <- py_to_r_ifneedbe(x_list[[i]])
    }
  }

  x_list
}

# `py_to_r.collections.abc.Set` ?
# `py_to_r.collections.abc.KeysView` ?
# `py_to_r.scipy.sparse.csc.csc_matrix` ?