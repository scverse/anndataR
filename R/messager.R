#' Print messages 
#' 
#' Conditionally print messages.
#'  Allows developers to easily control verbosity of functions, 
#'  and meet Bioconductor requirements that dictate the message 
#'  must first be stored to a variable before passing to \link[base]{message}. 
#'  
#' @param v Whether to print messages or not.
#' @param parallel Whether to enable message print when wrapped 
#' in parallelised functions.
#' @inheritParams base::paste
#' @return Null 
#' @keywords internal 
messager <- function(..., 
                     v = Sys.getenv("ANNDATAR_VERBOSE")!=FALSE, 
                     sep=" ",
                     collapse=NULL) { 
  msg <- paste(..., sep=sep, collapse=collapse)
  if (v) try({message(msg)})
}
