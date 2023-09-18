pretty_print <- function(label, value) {
  paste0(label, ": ", paste0("'", value, "'", collapse = ", "))
}

wrap_message <- function(...) {
  txt <- paste0(..., collapse = "")
  paste(strwrap(txt, exdent = 2L), collapse = "\n")
}
