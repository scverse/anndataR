#' Generator Types
#'
#' A list of available generator types for different components of a dataset
#' that can be passed to arguments of [generate_dataset()]
#'
#' @format A named list with the following elements:
#'
#' - `X`: Types of matrices that can be used for the main data matrix
#' - `layers`: Types of matrices that can be used for additional layers
#' - `obs`: Types of vectors that can be used for observation metadata
#' - `var`: Types of vectors that can be used for variable metadata
#' - `obsm`: Types of matrices or vectors that can be used for observation
#'   multi-dimensional annotations
#' - `varm`: Types of matrices or vectors that can be used for variable
#'   multi-dimensional annotations
#' - `obsp`: Types of matrices that can be used for pairwise observation
#'   annotations
#' - `varp`: Types of matrices that can be used for pairwise variable
#'   annotations
#' - `uns`: Types of miscellaneous data that can be used for unstructured
#'   annotations
#'
#' @details
#' The `generator_types` list contains all available types and
#' `example_generator_types` contains only those that are used when
#' `example = TRUE` is passed to [generate_dataset()].
#'
#' @seealso [generate_dataset()] for the generator function
#'
#' @examples
#' generator_types
#'
#' @export
generator_types <- list(
  X = names(matrix_generators),
  layers = names(matrix_generators),
  obs = names(vector_generators),
  var = names(vector_generators),
  obsm = c(names(matrix_generators), names(vector_generators)),
  varm = c(names(matrix_generators), names(vector_generators)),
  obsp = names(matrix_generators),
  varp = names(matrix_generators),
  uns = c(
    paste0("scalar_", names(vector_generators)),
    paste0("vec_", names(vector_generators)),
    paste0("df_", names(vector_generators)),
    paste0("mat_", names(matrix_generators)),
    "list"
  )
)

#' @rdname generator_types
#' @export
#'
#' @format NULL
#'
#' @examples
#' example_generator_types
example_generator_types <- list(
  X = "numeric_matrix",
  layers = c("numeric_matrix", "numeric_dense", "numeric_csparse"),
  obs = c("character", "integer", "factor"),
  var = c("character", "integer", "factor"),
  obsm = c("numeric_matrix", "numeric_dense", "numeric_csparse"),
  varm = c("numeric_matrix", "numeric_dense", "numeric_csparse"),
  obsp = c("numeric_matrix", "numeric_dense", "numeric_csparse"),
  varp = c("numeric_matrix", "numeric_dense", "numeric_csparse"),
  uns = c(
    "scalar_character",
    "vec_integer",
    "df_logical",
    "mat_numeric_matrix"
  )
)
