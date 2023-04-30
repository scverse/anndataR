pretty_print <- function(label, value) {
  txt <- paste0(label, ": ", paste(value, collapse = " "))
  paste0(strwrap(txt, indent = 0, exdent = 2), collapse = "\n")
}

wrap_message <- function(...) {
  txt <- paste0(..., collapse = "")
  paste(strwrap(txt, exdent = 2L), collapse = "\n")
}
