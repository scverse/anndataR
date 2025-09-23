py_to_r.collections.abc.Mapping <- function(x, ...) {
  check_requires("Converting Python to R", "reticulate")

  out <- list()
  bi <- reticulate::import_builtins()
  keys <- bi$list(x$keys())
  for (name in keys) {
    out[[name]] <- py_to_r(x[[name]], ...)
  }
  out
}
