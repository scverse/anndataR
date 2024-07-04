wrap_message <- function(...) {
  txt <- paste0(..., collapse = "")
  paste(strwrap(txt, exdent = 2L), collapse = "\n")
}

has_row_names <- function(x) {
  if (is.data.frame(x)) {
    .row_names_info(x) > 0
  } else {
    !is.null(dimnames(x)[[1]])
  }
}

get_shape <- function(obs, var, X, shape) {
  n_obs <-
    if (!is.null(obs)) {
      nrow(obs)
    } else if (!is.null(X)) {
      nrow(X)
    } else if (!is.null(shape)) {
      shape[[1]]
    } else {
      0L
    }
  n_vars <-
    if (!is.null(var)) {
      nrow(var)
    } else if (!is.null(X)) {
      ncol(X)
    } else if (!is.null(shape)) {
      shape[[2]]
    } else {
      0L
    }
  c(n_obs, n_vars)
}

get_initial_obs <- function(obs, X, shape) {
  if (is.null(obs)) {
    obs <- data.frame(matrix(NA, nrow = shape[[1]], ncol = 0))
    if (!is.null(X)) {
      rownames(obs) <- rownames(X)
    }
  }
  obs
}

get_initial_var <- function(var, X, shape) {
  if (is.null(var)) {
    var <- data.frame(matrix(NA, nrow = shape[[2]], ncol = 0))
    if (!is.null(X)) {
      colnames(var) <- colnames(X)
    }
  }
  var
}