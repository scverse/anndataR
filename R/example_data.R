#' Example AnnData
#' 
#' Create a small example AnnData object.
#' @param n_obs Number of observations.
#' @param n_var Number of variables.
#' @param output_class Name of the AnnData class. 
#' Must be one of:
#' \itemize{
#' \item{"InMemoryAnnData": }{Produces \link[anndataR]{InMemoryAnnData}
#'  (default).}
#' \item{"HDF5AnnData": }{Produces \link[anndataR]{HDF5AnnData}.}
#' \item{"Seurat": }{Produces \link[SeuratObject]{Seurat}.}
#' }
#' @param ... Additional arguments passed to the generator function.
#' See the "Details" section for more information on which parameters.
#' @returns \link[anndataR]{InMemoryAnnData} or 
#' \link[anndataR]{HDF5AnnData} or link[SeuratObject]{Seurat}.
#' @export
#' @examples
#' example_anndata()
example_data <- function(n_obs=3L,
                         n_var=5L,
                         output_class=c("InMemoryAnnData",
                                        "HDF5AnnData",
                                        "Seurat"),
                         file=tempfile(),
                         ...){
  
  output_class <- output_class[1]
  generator <- get_generator(
    if(output_class=="Seurat") {
      eval(formals(example_data)$output_class)[1]
    } else {
      output_class
    }
  )
  
  if(output_class=="HDF5AnnData"){
    ad <- generator$new(
      X = matrix(seq(n_var), n_obs, n_var),
      obs = data.frame(cell = seq(n_obs)),
      obs_names = letters[seq(n_obs)],
      var = data.frame(gene = seq(n_var)),
      var_names = letters[seq(n_var)],
      file=file,
      ...
    )
  } else {
    ad <- generator$new(
      X = matrix(seq(n_var), n_obs, n_var),
      obs = data.frame(cell = seq(n_obs)),
      obs_names = letters[seq(n_obs)],
      var = data.frame(gene = seq(n_var)),
      var_names = letters[seq(n_var)], 
      ...
    )
  }
  #### Return ####
  if(output_class=="Seurat"){ 
    return(to_Seurat(ad))
  } else {
    return(ad)
  }
}