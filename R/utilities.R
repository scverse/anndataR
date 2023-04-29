pretty_print <- function(label, value) {
  txt <- paste0(label, ": ", paste(value, collapse = " "))
  paste0(strwrap(txt, indent = 0, exdent = 2), collapse = "\n")
}
