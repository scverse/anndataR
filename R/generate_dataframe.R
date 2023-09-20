#' Generate a dataframe
#'
#' Generate a dataframe with different types of columns
#'
#' @param num_rows Number of rows to generate
#' @param types Types of columns to generate
#'
#' @return A dataframe with the generated columns
#'
#' @examples
#' generate_dataframe(10L)
generate_dataframe <- function(num_rows, types = names(vector_generators)) {
  data <- lapply(types, generate_vector, n = num_rows)
  names(data) <- types
  as.data.frame(data)
}
