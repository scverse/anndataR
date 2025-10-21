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
      rownames(var) <- colnames(X)
    }
  }
  var
}

to_py_matrix <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "RsparseMatrix")
  } else if (!inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }
  Matrix::t(mat)
}

# nolint start: object_name_linter
to_R_matrix <- function(mat) {
  # nolint end: object_name_linter
  if (inherits(mat, "dgRMatrix")) {
    mat <- as(mat, "CsparseMatrix")
  } else if (!inherits(mat, "dgCMatrix")) {
    mat <- as.matrix(mat)
  }
  Matrix::t(mat)
}

self_name <- function(x) {
  if (is.null(names(x))) {
    x <- setNames(x, x)
  } else if (any(names(x) == "")) {
    is_missing <- names(x) == ""
    names(x)[is_missing] <- x[is_missing]
  }

  x
}

#' Get mapping
#'
#' Get a mapping argument for a conversion function
#'
#' @param mapping The user-supplied mapping argument. Can be a named vector,
#'   `TRUE` or `FALSE`.
#' @param guesser A function that guesses the default mapping from `obj` if
#'   `mapping` is `TRUE`
#' @param obj The object that is being converted and is passed to `guesser`
#'   if needed
#' @param name The name of the mapping argument, used for error messages
#' @param ... Additional arguments passed to `guesser`
#'
#' @description
#' If `mapping` is `NULL` or empty it is set to `FALSE` with a warning. `FALSE`
#' values return an empty mapping.
#'
#' @returns A named mapping vector
#' @noRd
get_mapping <- function(mapping, guesser, obj, name, ...) {
  if (rlang::is_empty(mapping)) {
    cli_warn(c(
      "The {.arg {name}} argument is empty, setting it to {.val {FALSE}}"
    ))

    mapping <- FALSE
  }

  # If FALSE, return an empty mapping
  if (isFALSE(mapping)) {
    return(list())
  }

  # If TRUE, use the guesser function to get the default mapping
  if (isTRUE(mapping)) {
    return(guesser(obj, ...))
  }

  if (!is.vector(mapping)) {
    cli_abort(paste(
      "{.arg {name}} must be a vector, {.val {TRUE}} or {.val {FALSE}}, not",
      "{.cls {class(mapping)}}"
    ))
  }

  # Make sure provided mapping has names
  self_name(mapping)
}

#' Check dimensions and skip
#'
#' Check the dimensions of a matrix-like object and return `NULL` if they do not
#' match the expected dimensions, with a warning. For use in conversion
#' functions to skip items that do not match the required dimensions.
#'
#' @param x The object to check
#' @param field The field the object comes from, used in the warning message
#' @param name The name of the object, used in the warning message
#' @param expected_dims Expected dimensions
#' @param expected_rows Expected number of rows
#' @param expected_cols Expected number of columns
#'
#' @returns The object `x` if it matches the expected dimensions, otherwise
#'   `NULL`
#' @noRd
check_dims_and_skip <- function(
  x,
  field,
  name,
  expected_dims = NULL,
  expected_rows = NULL,
  expected_cols = NULL
) {
  msg <- NULL

  if (!is.null(expected_dims) && !identical(dim(x), expected_dims)) {
    expected_dims <- as.integer(expected_dims)
    msg <- c(
      "i" = paste0(
        "Expected [{style_vec(expected_dims)}], ",
        "got [{style_vec(as.integer(dim(x)))}]"
      )
    )
  } else if (!is.null(expected_rows) && nrow(x) != expected_rows) {
    msg <- c(
      "i" = paste0(
        "Expected {.val {expected_rows}} rows, got {.val {nrow(x)}}"
      )
    )
  } else if (!is.null(expected_cols) && ncol(x) != expected_cols) {
    msg <- c(
      "i" = paste0(
        "Expected {.val {expected_cols}} colums, got {.val {ncol(x)}}"
      )
    )
  }

  if (!is.null(msg)) {
    cli_warn(
      c(
        "Skipping {.field {field}} {.val {name}} with unexpected dimensions",
        msg
      ),
      call = NULL
    )

    return(NULL)
  } else {
    return(x)
  }
}
