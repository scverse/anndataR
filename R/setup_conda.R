#' Setup conda
#' 
#' Install miniconda and setup a conda env to use python dependencies.
#' @inheritParams reticulate::install_miniconda
#' @inheritParams reticulate::py_install
#' @inheritDotParams reticulate::py_install
#' @returns Null
#' 
#' @export
#' @examples
#' \dontrun{
#' setup_conda()
#' }
setup_conda <- function(path = reticulate::miniconda_path(),
                        update = TRUE, 
                        force = FALSE,
                        packages = c("anndata", "scanpy"),
                        ...){
  
  requireNamespace("reticulate")
  if(!reticulate::py_module_available(module = "anndata")){
    reticulate::install_miniconda(path = path, 
                                  update = update, 
                                  force = force)
    reticulate::py_install(packages = packages, 
                           pip = TRUE,
                           ...)
  } 
}
