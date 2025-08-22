#' Conversion utilities for AnnData objects between R and Python
#'
#' These functions provide utilities for converting between R and Python
#' representations of AnnData objects and their components.
#'
#' @name anndata_conversion_utils
NULL

#' Convert R objects to Python for use in AnnData
#'
#' This is a utility function for converting R objects to appropriate Python
#' representations for use in AnnData objects. For converting complete AnnData
#' objects, consider using `reticulate::r_to_py()` directly, which will automatically
#' use the appropriate S3 method when anndataR is loaded.
#'
#' @param x An R object to convert to Python
#' @param ... Additional arguments passed to reticulate conversion functions
#'
#' @return A Python object
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Requires Python and reticulate
#' mat <- matrix(1:12, 3, 4)
#' py_mat <- try_r_to_py(mat)
#'
#' # For AnnData objects, prefer the automatic method:
#' adata <- AnnData(X = mat)
#' py_adata <- r_to_py(adata)  # Automatically converts
#' }
try_r_to_py <- function(x, ...) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli_abort("The {.pkg reticulate} package is required")
  }

  if (is.null(x)) {
    return(reticulate::py_none())
  }

  # Handle different R object types appropriately for AnnData
  if (is.data.frame(x)) {
    # Convert data frame to pandas DataFrame
    reticulate::r_to_py(x)
  } else if (is.matrix(x) || inherits(x, "Matrix")) {
    # Handle matrices - may need special treatment for sparse matrices
    if (inherits(x, "dgCMatrix") || inherits(x, "dgRMatrix")) {
      # For sparse matrices, convert directly using reticulate
      reticulate::r_to_py(x)
    } else {
      # Regular dense matrix
      reticulate::r_to_py(x)
    }
  } else if (is.list(x) && !is.null(names(x))) {
    # Convert named list to Python dict
    py_dict <- reticulate::dict()
    for (name in names(x)) {
      py_dict[[name]] <- try_r_to_py(x[[name]], ...)
    }
    py_dict
  } else if (is.list(x)) {
    # Convert unnamed list to Python dict if empty, otherwise to Python list
    if (length(x) == 0) {
      # Empty list should be a dict for AnnData compatibility
      reticulate::dict()
    } else {
      # Non-empty unnamed list becomes Python list
      reticulate::r_to_py(x, ...)
    }
  } else {
    # For other types, use standard reticulate conversion
    reticulate::r_to_py(x, ...)
  }
}

#' Convert Python objects to R for use with AnnData
#'
#' This is a utility function for converting Python objects to appropriate R
#' representations for use with AnnData objects. For converting complete AnnData
#' objects from Python, consider using `reticulate::py_to_r()` directly, which will
#' automatically use the appropriate S3 method when anndataR is loaded.
#'
#' @param x A Python object to convert to R
#' @param ... Additional arguments passed to reticulate conversion functions
#'
#' @return An R object
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Requires Python and reticulate
#' # For individual objects:
#' py_mat <- r_to_py(matrix(1:12, 3, 4))
#' r_mat <- try_py_to_r(py_mat)
#'
#' # For AnnData objects, prefer the automatic method:
#' library(reticulate)
#' ad <- import("anndata", convert = FALSE)
#' py_adata <- ad$AnnData(X = py_mat)
#' r_adata <- py_to_r(py_adata)  # Automatically converts to ReticulateAnnData
#' }
try_py_to_r <- function(x, ...) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli_abort("The {.pkg reticulate} package is required")
  }

  if (is.null(x) || identical(x, reticulate::py_none())) {
    return(NULL)
  }

  # Check if it's actually a Python object first
  if (!inherits(x, "python.builtin.object")) {
    return(x) # Already an R object
  }

  # Handle different Python object types appropriately for AnnData
  if (inherits(x, "pandas.core.frame.DataFrame")) {
    # Convert pandas DataFrame to R data.frame
    reticulate::py_to_r(x)
  } else if (inherits(x, "numpy.ndarray")) {
    # Convert numpy array to R matrix/array
    reticulate::py_to_r(x)
  } else if (inherits(x, "scipy.sparse.base.spmatrix")) {
    # Convert scipy sparse matrix to R sparse matrix
    reticulate::py_to_r(x)
  } else {
    # Handle AnnData special mapping objects (layers, obsm, varm, etc.)
    if (grepl("anndata.*mapping", paste(class(x), collapse = " "))) {
      tryCatch(
        {
          # Convert to regular dict first, then to R list
          py_dict <- reticulate::py_eval("dict")(x)
          return(reticulate::py_to_r(py_dict))
        },
        error = function(e) {
          # Fallback: create empty list
          return(list())
        }
      )
    }

    # Handle other dictionary-like objects
    if (inherits(x, "python.builtin.dict")) {
      # Convert Python dict to R named list using bi$list() for keys
      bi <- reticulate::import_builtins()
      dict_keys <- reticulate::py_to_r(bi$list(x$keys()))
      result <- list()
      for (key in dict_keys) {
        result[[key]] <- try_py_to_r(x[[key]], ...)
      }
      result
    } else {
      # Try to convert other dictionary-like objects to named list
      tryCatch(
        {
          if (reticulate::py_has_attr(x, "keys")) {
            # Use import_builtins to convert keys that reticulate can't handle
            bi <- reticulate::import_builtins()
            keys <- tryCatch(
              {
                reticulate::py_to_r(bi$list(x$keys()))
              },
              error = function(e) {
                # Fallback to direct py_to_r
                reticulate::py_to_r(x$keys())
              }
            )

            result <- list()
            for (key in keys) {
              result[[key]] <- try_py_to_r(x[[key]], ...)
            }
            return(result)
          } else {
            # Fallback to direct conversion
            return(reticulate::py_to_r(x, ...))
          }
        },
        error = function(e) {
          # Final fallback to direct conversion
          return(reticulate::py_to_r(x, ...))
        }
      )
    }
  }
}
