#' @description
#' Use `get_generator_types()` to get the available types for each slot.
#'
#' @param slot Which slot to return types for, if `NULL` a named list of all
#'   slots is returned
#'
#' @details
#' Use `get_generator_types()` to get a list of the available types for each
#' slot, or for a specific slot by setting `slot`. If `example = TRUE`, only the
#' example types are returned.
#'
#' @returns For `get_generator_types()`, a named list of character vectors or a
#'   single character vector if `slot` is not `NULL`
#' @export
#' @rdname generate_dataset
#'
#' @examples
#' # Get all available generator types
#' get_generator_types()
#'
#' # Get generator types for a specific slot
#' get_generator_types(slot = "obs")
#'
#' # Get generator types used when example = TRUE
#' get_generator_types(example = TRUE)
get_generator_types <- function(example = FALSE, slot = NULL) {
  if (example) {
    types <- .example_generator_types
  } else {
    types <- .generator_types
  }

  if (!is.null(slot)) {
    if (!(slot %in% names(types))) {
      cli_abort(c(
        "Invalid slot: {.val {slot}}",
        "i" = "Must be one of {.or {.val {names(types)}}}"
      ))
    }
    types <- types[[slot]]
  }

  types
}

.generator_types <- list(
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

.example_generator_types <- list(
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
