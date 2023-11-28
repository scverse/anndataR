wrap_message <- function(...) {
  txt <- paste0(..., collapse = "")
  paste(strwrap(txt, exdent = 2L), collapse = "\n")
}
