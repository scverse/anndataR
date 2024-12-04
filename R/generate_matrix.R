# nolint start
generate_numeric_matrix <- function(n_obs, n_vars, NAs = FALSE) {
  # byrow = TRUE to mimic the way a matrix gets filled in Python
  m <- matrix(seq(0.5, n_obs * n_vars), nrow = n_obs, ncol = n_vars, byrow = TRUE) 
  if (NAs) {
    m[1, 1] <- NA_real_
  }
  m
}

generate_integer_matrix <- function(n_obs, n_vars, NAs = FALSE) {
  # byrow = TRUE to mimic the way a matrix gets filled in Python
  m <- matrix(seq(0L, n_obs * n_vars), nrow = n_obs, ncol = n_vars, byrow = TRUE)
  if (NAs) {
    m[1, 1] <- NA_integer_
  }
  m
}

# nolint start
matrix_generators <- list(
  numeric_matrix = function(n_obs, n_vars) {
    generate_numeric_matrix(n_obs, n_vars)
  },
  numeric_dense = function(n_obs, n_vars) {
    m <- generate_numeric_matrix(n_obs, n_vars)
    as(m, "denseMatrix")
  },
  numeric_csparse = function(n_obs, n_vars) {
    m <- generate_numeric_matrix(n_obs, n_vars)
    as(m, "CsparseMatrix")
  },
  numeric_rsparse = function(n_obs, n_vars) {
    m <- generate_numeric_matrix(n_obs, n_vars)
    as(m, "RsparseMatrix")
  },
  numeric_matrix_with_nas = function(n_obs, n_vars) {
    generate_numeric_matrix(n_obs, n_vars, NAs = TRUE)
  },
  numeric_dense_with_nas = function(n_obs, n_vars) {
    m <- generate_numeric_matrix(n_obs, n_vars, NAs = TRUE)
    as(m, "denseMatrix")
  },
  numeric_csparse_with_nas = function(n_obs, n_vars) {
    m <- generate_numeric_matrix(n_obs, n_vars, NAs = TRUE)
    as(m, "CsparseMatrix")
  },
  numeric_rsparse_with_nas = function(n_obs, n_vars) {
    m <- generate_numeric_matrix(n_obs, n_vars, NAs = TRUE)
    as(m, "RsparseMatrix")
  },
  integer_matrix = function(n_obs, n_vars) {
    generate_integer_matrix(n_obs, n_vars)
  },
  integer_dense = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars)
    as(m, "denseMatrix")
  },
  integer_csparse = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars)
    as(m, "CsparseMatrix")
  },
  integer_rsparse = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars)
    as(m, "RsparseMatrix")
  },
  integer_matrix_with_nas = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars, NAs = TRUE)
    m
  },
  integer_dense_with_nas = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars, NAs = TRUE)
    as(m, "denseMatrix")
  },
  integer_csparse_with_nas = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars, NAs = TRUE)
    as(m, "CsparseMatrix")
  },
  integer_rsparse_with_nas = function(n_obs, n_vars) {
    m <- generate_integer_matrix(n_obs, n_vars, NAs = TRUE)
    as(m, "RsparseMatrix")
  }
)
# nolint end

#' Generate a matrix
#'
#' Generate a matrix of a given type
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#'
#' @return A matrix of the given type
#'
#' @noRd
#'
#' @examples
#' generate_matrix(10L, 20L)
generate_matrix <- function(n_obs, n_vars, type = names(matrix_generators)) {
  type <- match.arg(type)
  matrix_generators[[type]](n_obs, n_vars)
}
