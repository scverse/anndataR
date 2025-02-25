vector_generators <- list(
  character = function(n) paste0("value_", seq(from = 0, to = n - 1)),
  integer = function(n) seq(from = 0, to = n - 1),
  factor = function(n) factor(rep(c("Value1", "Value2"), length.out = n)),
  factor_ordered = function(n)
    factor(rep(c("Value1", "Value2"), length.out = n), ordered = TRUE),
  logical = function(n) sample(c(TRUE, FALSE), n, replace = TRUE),
  numeric = function(n) seq(from = 0.5, to = n),
  character_with_nas = function(n) {
    x <- paste0("value", seq_len(n))
    x[seq(1, n, by = 2)] <- NA_character_
    x
  },
  integer_with_nas = function(n) {
    x <- seq(from = 0, to = n - 1)
    x[1] <- NA_integer_
    x
  },
  factor_with_nas = function(n) {
    x <- factor(rep(c("Value1", "Value2"), length.out = n))
    x[1] <- NA_character_
    x
  },
  factor_ordered_with_nas = function(n) {
    x <- factor(rep(c("Value1", "Value2"), length.out = n), ordered = TRUE)
    x[1] <- NA_character_
    x
  },
  logical_with_nas = function(n) {
    x <- sample(c(TRUE, FALSE), n, replace = TRUE)
    x[seq(1, n, by = 2)] <- NA
    x
  },
  numeric_with_nas = function(n) {
    x <- runif(n)
    x[seq(1, n, by = 2)] <- NA_real_
    x
  }
)

#' Generate a vector
#'
#' Generate a vector of a given type
#'
#' @param n Number of elements to generate
#' @param type Type of vector to generate
#'
#' @return A vector of the given type
#'
#' @noRd
#'
#' @examples
#' generate_vector(10L)
generate_vector <- function(n, type = names(vector_generators)) {
  type <- match.arg(type)
  vector_generators[[type]](n)
}
