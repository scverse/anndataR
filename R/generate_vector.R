vector_generators <- list(
  character = function(n) paste0("value", seq_len(n)),
  integer = function(n) seq_len(n),
  factor = function(n) factor(paste0("value", seq_len(n))),
  factor_ordered = function(n) factor(paste0("value", seq_len(n)), ordered = TRUE),
  logical = function(n) sample(c(TRUE, FALSE), n, replace = TRUE),
  numeric = function(n) runif(n),
  character_with_nas = function(n, pct = .5) {
    x <- paste0("value", seq_len(n))
    x[sample.int(n, pct * n)] <- NA_character_
    x
  },
  integer_with_nas = function(n, pct = .5) {
    x <- seq_len(n)
    x[sample.int(n, pct * n)] <- NA_integer_
    x
  },
  factor_with_nas = function(n, pct = .5) {
    x <- factor(paste0("value", seq_len(n)))
    x[sample.int(n, pct * n)] <- NA_character_
    x
  },
  factor_ordered_with_nas = function(n, pct = .5) {
    x <- factor(paste0("value", seq_len(n)), ordered = TRUE)
    x[sample.int(n, pct * n)] <- NA_character_
    x
  },
  logical_with_nas = function(n, pct = .5) {
    x <- sample(c(TRUE, FALSE), n, replace = TRUE)
    x[sample.int(n, pct * n)] <- NA
    x
  },
  numeric_with_nas = function(n, pct = .5) {
    x <- runif(n)
    x[sample.int(n, pct * n)] <- NA_real_
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
