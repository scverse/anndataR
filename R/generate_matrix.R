# nolint start
matrix_generators <- list(
  numeric_matrix = function(n_obs, n_vars) {
    matrix(runif(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  },
  numeric_dense = function(n_obs, n_vars) {
    m <- matrix(runif(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
    as(m, "denseMatrix")
  },
  numeric_csparse = function(n_obs, n_vars) {
    m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
    as(m, "CsparseMatrix")
  },
  numeric_rsparse = function(n_obs, n_vars) {
    m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
    as(m, "RsparseMatrix")
  },
  # TODO: re-enable
  # numeric_matrix_with_nas = function(n_obs, n_vars) {
  #   m <- matrix(runif(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_real_
  #   m
  # },
  # numeric_dense_with_nas = function(n_obs, n_vars) {
  #   m <- matrix(runif(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_real_
  #   as(m, "denseMatrix")
  # },
  # numeric_csparse_with_nas = function(n_obs, n_vars) {
  #   m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_real_
  #   as(m, "CsparseMatrix")
  # },
  # numeric_rsparse_with_nas = function(n_obs, n_vars) {
  #   m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_real_
  #   as(m, "RsparseMatrix")
  # },
  integer_matrix = function(n_obs, n_vars) {
    matrix(sample.int(100L, n_obs * n_vars, replace = TRUE), nrow = n_obs, ncol = n_vars)
  },
  integer_dense = function(n_obs, n_vars) {
    m <- matrix(sample.int(100L, n_obs * n_vars, replace = TRUE), nrow = n_obs, ncol = n_vars)
    as(m, "denseMatrix")
  },
  integer_csparse = function(n_obs, n_vars) {
    m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
    as(m, "CsparseMatrix")
  },
  integer_rsparse = function(n_obs, n_vars) {
    m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
    as(m, "RsparseMatrix")
  } # ,
  # TODO: re-enable
  # integer_matrix_with_nas = function(n_obs, n_vars) {
  #   m <- matrix(sample.int(100L, n_obs * n_vars, replace = TRUE), nrow = n_obs, ncol = n_vars)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_integer_
  #   m
  # },
  # integer_dense_with_nas = function(n_obs, n_vars) {
  #   m <- matrix(sample.int(100L, n_obs * n_vars, replace = TRUE), nrow = n_obs, ncol = n_vars)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_integer_
  #   as(m, "denseMatrix")
  # },
  # integer_csparse_with_nas = function(n_obs, n_vars) {
  #   m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_integer_
  #   as(m, "CsparseMatrix")
  # },
  # integer_rsparse_with_nas = function(n_obs, n_vars) {
  #   m <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)
  #   m[seq(1, n_obs * n_vars, by = 2)] <- NA_integer_
  #   as(m, "RsparseMatrix")
  # }
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
