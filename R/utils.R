wrap_message <- function(...) {
  txt <- paste0(..., collapse = "")
  paste(strwrap(txt, exdent = 2L), collapse = "\n")
}

has_row_names <- function(x) {
  if (is.data.frame(x)) {
    .row_names_info(x) > 0
  } else {
    !is.null(dimnames(x)[[1]])
  }
}
