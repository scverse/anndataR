#' Helper functions for concatenating AnnData objects
#'
#' This file contains internal functions used for concatenating AnnData objects,
#' including reindexing, matrix concatenation, and annotation merging.
#'
#' Related files:
#' - concat.R: Main concatenation interface (concat, rbind, cbind)
#' - concat_strategies.R: Merge strategies for handling non-aligned elements
#'
#' @keywords internal
#' @name concat_helpers

#' Reindexer class for handling dimension alignment during concatenation
#'
#' This class handles the complex logic of reindexing matrices and data frames
#' when concatenating AnnData objects with mismatched dimensions.
#'
#' @field old_idx Original index (character vector)
#' @field new_idx Target index (character vector)
#' @field old_pos Positions in original index that will be kept (integer vector)
#' @field new_pos Positions in new index where data will be placed (integer vector)
#' @field no_change Whether indices are identical (logical)
#'
#' @keywords internal
#' @noRd
Reindexer <- R6::R6Class(
  "Reindexer",
  public = list(
    old_idx = NULL,
    new_idx = NULL,
    old_pos = NULL,
    new_pos = NULL,
    no_change = NULL,

    #' @description
    #' Create a new Reindexer
    #' @param old_idx Original index
    #' @param new_idx Target index
    initialize = function(old_idx, new_idx) {
      self$old_idx <- as.character(old_idx)
      self$new_idx <- as.character(new_idx)
      self$no_change <- identical(self$old_idx, self$new_idx)

      if (!self$no_change) {
        # Find positions for reindexing (1-based indexing in R)
        new_pos <- match(self$old_idx, self$new_idx)
        old_pos <- seq_along(new_pos)

        # Keep only valid matches (remove NA positions)
        mask <- !is.na(new_pos)
        self$new_pos <- new_pos[mask]
        self$old_pos <- old_pos[mask]
      }
    },

    #' @description
    #' Apply reindexing to an element
    #' @param el Element to reindex (matrix, data.frame, etc.)
    #' @param axis Axis to reindex along (0 for rows, 1 for columns)
    #' @param fill_value Value to use for missing positions
    apply = function(el, axis = 1L, fill_value = NULL) {
      if (self$no_change) {
        return(el)
      }

      if (is.data.frame(el)) {
        return(self$apply_to_dataframe(el, axis, fill_value))
      } else if (inherits(el, "sparseMatrix") || methods::is(el, "Matrix")) {
        return(self$apply_to_sparse(el, axis, fill_value))
      } else if (is.matrix(el) || is.array(el)) {
        return(self$apply_to_dense(el, axis, fill_value))
      } else {
        cli::cli_abort("Cannot reindex object of class {class(el)[1]}")
      }
    },

    #' @description Apply reindexing to data.frame
    apply_to_dataframe = function(el, axis, fill_value = NULL) {
      if (is.null(fill_value)) {
        fill_value <- NA
      }

      if (axis == 0L) {
        # Reindex rows
        result <- el[rep(NA_integer_, length(self$new_idx)), , drop = FALSE]
        rownames(result) <- self$new_idx
        if (length(self$new_pos) > 0) {
          result[self$new_pos, ] <- el[self$old_pos, , drop = FALSE]
        }
      } else {
        # Reindex columns
        result <- data.frame(matrix(
          fill_value,
          nrow = nrow(el),
          ncol = length(self$new_idx)
        ))
        colnames(result) <- self$new_idx
        rownames(result) <- rownames(el)
        if (length(self$new_pos) > 0) {
          result[, self$new_pos] <- el[, self$old_pos, drop = FALSE]
        }
      }

      return(result)
    },

    #' @description Apply reindexing to sparse matrix
    apply_to_sparse = function(el, axis, fill_value = NULL) {
      if (is.null(fill_value)) {
        fill_value <- 0
      }

      if (axis == 0L) {
        # Reindex rows
        new_nrow <- length(self$new_idx)
        if (length(self$new_pos) == 0) {
          # No matching rows - return empty sparse matrix
          result <- Matrix::sparseMatrix(
            i = integer(0),
            j = integer(0),
            x = numeric(0),
            dims = c(new_nrow, ncol(el))
          )
        } else {
          # Create indexing matrix for efficient sparse reindexing
          idx_i <- rep(self$new_pos, each = ncol(el))
          idx_j <- rep(seq_len(ncol(el)), times = length(self$new_pos))

          # Extract values from original positions
          orig_vals <- as.matrix(el[self$old_pos, , drop = FALSE])
          vals <- as.numeric(orig_vals)

          result <- Matrix::sparseMatrix(
            i = idx_i,
            j = idx_j,
            x = vals,
            dims = c(new_nrow, ncol(el))
          )
        }
        rownames(result) <- self$new_idx
        colnames(result) <- colnames(el)
      } else {
        # Reindex columns
        new_ncol <- length(self$new_idx)
        if (length(self$new_pos) == 0) {
          # No matching columns - return empty sparse matrix
          result <- Matrix::sparseMatrix(
            i = integer(0),
            j = integer(0),
            x = numeric(0),
            dims = c(nrow(el), new_ncol)
          )
        } else {
          # Create indexing matrix for efficient sparse reindexing
          idx_i <- rep(seq_len(nrow(el)), times = length(self$new_pos))
          idx_j <- rep(self$new_pos, each = nrow(el))

          # Extract values from original positions
          orig_vals <- as.matrix(el[, self$old_pos, drop = FALSE])
          vals <- as.numeric(orig_vals)

          result <- Matrix::sparseMatrix(
            i = idx_i,
            j = idx_j,
            x = vals,
            dims = c(nrow(el), new_ncol)
          )
        }
        rownames(result) <- rownames(el)
        colnames(result) <- self$new_idx
      }

      return(result)
    },

    #' @description Apply reindexing to dense matrix/array
    apply_to_dense = function(el, axis, fill_value = NULL) {
      if (is.null(fill_value)) {
        fill_value <- NA
      }

      if (axis == 0L) {
        # Reindex rows
        result <- array(fill_value, dim = c(length(self$new_idx), ncol(el)))
        if (length(self$new_pos) > 0) {
          result[self$new_pos, ] <- el[self$old_pos, , drop = FALSE]
        }
        rownames(result) <- self$new_idx
        colnames(result) <- colnames(el)
      } else {
        # Reindex columns
        result <- array(fill_value, dim = c(nrow(el), length(self$new_idx)))
        if (length(self$new_pos) > 0) {
          result[, self$new_pos] <- el[, self$old_pos, drop = FALSE]
        }
        rownames(result) <- rownames(el)
        colnames(result) <- self$new_idx
      }

      return(result)
    }
  )
)

#' Create a reindexer for aligning dimensions
#' @param new_idx Target index
#' @param cur_idx Current index
#' @return A Reindexer object
#' @keywords internal
#' @noRd
gen_reindexer <- function(new_idx, cur_idx) {
  Reindexer$new(cur_idx, new_idx)
}

#' Merge indices using specified join strategy
#' @param indices_list List of index vectors
#' @param join Join strategy ("inner" or "outer")
#' @return Merged index vector
#' @keywords internal
#' @noRd
merge_indices <- function(indices_list, join = "outer") {
  if (join == "inner") {
    # Intersection of all indices
    result <- indices_list[[1]]
    for (i in seq_along(indices_list)[-1]) {
      result <- intersect(result, indices_list[[i]])
    }
    return(result)
  } else if (join == "outer") {
    # Union of all indices, preserving order
    result <- character(0)
    for (idx in indices_list) {
      result <- union(result, idx)
    }
    return(result)
  } else {
    cli::cli_abort("Join must be 'inner' or 'outer', not {join}")
  }
}

#' @noRd
gen_outer_reindexers <- function(indices_list, merged_index) {
  lapply(indices_list, function(idx) gen_reindexer(merged_index, idx))
}

#' @noRd
gen_inner_reindexers <- function(indices_list) {
  common_idx <- merge_indices(indices_list, join = "inner")
  lapply(indices_list, function(idx) gen_reindexer(common_idx, idx))
}

#' @noRd
default_fill_value <- function(elements) {
  # Check if any element is sparse - use 0 for sparse matrices
  if (
    any(sapply(elements, function(x) {
      if (is.null(x)) {
        return(FALSE)
      }
      inherits(x, "sparseMatrix") ||
        methods::is(x, "sparseMatrix") ||
        (methods::is(x, "Matrix") && attr(class(x), "package") == "Matrix")
    }))
  ) {
    return(0)
  } else {
    return(NA)
  }
}

#' @noRd
concat_X <- function(adatas, reindexers, axis, fill_value) {
  X_matrices <- lapply(adatas, function(ad) ad$X)

  # Check if all X are NULL
  all_null <- all(sapply(X_matrices, is.null))
  any_null <- any(sapply(X_matrices, is.null))

  if (all_null) {
    return(NULL)
  } else if (any_null) {
    cli::cli_abort(
      "Cannot concatenate: some (but not all) AnnData objects have X = NULL"
    )
  }

  # Apply reindexing to align dimensions
  reindexed <- Map(
    function(X, reindexer) {
      if (axis == 0L) {
        # Concatenating rows - reindex columns
        reindexer$apply(X, axis = 1L, fill_value = fill_value)
      } else {
        # Concatenating columns - reindex rows
        reindexer$apply(X, axis = 0L, fill_value = fill_value)
      }
    },
    X_matrices,
    reindexers
  )

  # Concatenate along specified axis
  if (axis == 0L) {
    # Stack rows
    do.call(rbind, reindexed)
  } else {
    # Stack columns
    do.call(cbind, reindexed)
  }
}

#' @noRd
concat_annotations <- function(annotations, new_index, join = "outer") {
  if (length(annotations) == 0) {
    return(data.frame())
  }

  # Remove empty annotations
  annotations <- annotations[sapply(annotations, function(x) nrow(x) > 0)]

  if (length(annotations) == 0) {
    return(data.frame(row.names = new_index))
  }

  # Concatenate data frames
  result <- do.call(rbind, annotations)

  # Set new index
  rownames(result) <- new_index

  return(result)
}

#' @noRd
merge_annotations <- function(annotations, new_index, strategy) {
  if (length(annotations) == 0) {
    return(data.frame())
  }

  # Apply merge strategy
  merged <- apply_merge_strategy(annotations, strategy)

  # Convert to data.frame with proper index
  if (length(merged) == 0) {
    result <- data.frame(row.names = new_index)
  } else {
    # Take first annotation as template for structure
    template <- annotations[[1]]
    result <- template[new_index, , drop = FALSE]

    # Replace columns with merged values
    for (col_name in names(merged)) {
      if (col_name %in% colnames(template)) {
        result[[col_name]] <- merged[[col_name]]
      }
    }
  }

  return(result)
}

#' @noRd
concat_layers <- function(adatas, reindexers, axis, join, fill_value, merge) {
  all_layers <- lapply(adatas, function(ad) ad$layers)

  # Get layer names based on join strategy
  layer_names <- merge_indices(lapply(all_layers, names), join = join)

  if (length(layer_names) == 0) {
    return(list())
  }

  result_layers <- list()

  for (layer_name in layer_names) {
    # Get matrices for this layer from all objects
    layer_matrices <- lapply(all_layers, function(layers) {
      if (layer_name %in% names(layers)) {
        layers[[layer_name]]
      } else {
        NULL
      }
    })

    # Handle case where not all objects have this layer
    if (any(sapply(layer_matrices, is.null))) {
      if (join == "inner") {
        next # Skip layers not present in all objects
      } else {
        # For outer join, create empty matrices for missing layers
        for (i in seq_along(layer_matrices)) {
          if (is.null(layer_matrices[[i]])) {
            # Create empty matrix with correct dimensions
            if (axis == 0L) {
              dims <- c(adatas[[i]]$n_obs(), length(reindexers[[i]]$new_idx))
            } else {
              dims <- c(length(reindexers[[i]]$new_idx), adatas[[i]]$n_vars())
            }
            layer_matrices[[i]] <- Matrix::sparseMatrix(
              i = integer(0),
              j = integer(0),
              x = numeric(0),
              dims = dims
            )
          }
        }
      }
    }

    # Apply reindexing and concatenate
    reindexed <- Map(
      function(mat, reindexer) {
        if (is.null(mat)) {
          return(NULL)
        }

        if (axis == 0L) {
          reindexer$apply(mat, axis = 1L, fill_value = fill_value)
        } else {
          reindexer$apply(mat, axis = 0L, fill_value = fill_value)
        }
      },
      layer_matrices,
      reindexers
    )

    # Remove NULL entries
    reindexed <- reindexed[!sapply(reindexed, is.null)]

    if (length(reindexed) > 0) {
      if (axis == 0L) {
        result_layers[[layer_name]] <- do.call(rbind, reindexed)
      } else {
        result_layers[[layer_name]] <- do.call(cbind, reindexed)
      }
    }
  }

  return(result_layers)
}

#' @noRd
concat_matrices_dict <- function(matrices_list, axis, join, fill_value) {
  # Get all matrix names
  all_names <- unique(unlist(lapply(matrices_list, names)))

  if (length(all_names) == 0) {
    return(list())
  }

  result <- list()

  for (mat_name in all_names) {
    # Get matrices for this name
    matrices <- lapply(matrices_list, function(mats) {
      if (mat_name %in% names(mats)) mats[[mat_name]] else NULL
    })

    # Skip if not all objects have this matrix for inner join
    if (join == "inner" && any(sapply(matrices, is.null))) {
      next
    }

    # Remove NULL entries and concatenate
    matrices <- matrices[!sapply(matrices, is.null)]
    if (length(matrices) > 0) {
      if (axis == 0L) {
        result[[mat_name]] <- do.call(rbind, matrices)
      } else {
        result[[mat_name]] <- do.call(cbind, matrices)
      }
    }
  }

  return(result)
}

#' @noRd
merge_matrices_dict <- function(matrices_list, reindexers, axis, strategy) {
  # Apply merge strategy to get surviving matrices
  merged <- apply_merge_strategy(matrices_list, strategy)

  if (length(merged) == 0) {
    return(list())
  }

  # Note: For merge operations, we typically don't need reindexing
  # as the merged matrices should already have compatible dimensions
  # But we include reindexers parameter for API consistency

  return(merged)
}

#' @noRd
concat_pairwise_dict <- function(
  pairwise_list,
  batch_sizes,
  axis = NULL,
  strategy
) {
  if (is.null(axis)) {
    # Not concatenating along this axis - merge using strategy
    return(apply_merge_strategy(pairwise_list, strategy))
  }

  # Get all matrix names
  all_names <- unique(unlist(lapply(pairwise_list, names)))

  if (length(all_names) == 0) {
    return(list())
  }

  result <- list()

  for (mat_name in all_names) {
    # Get matrices for this name
    matrices <- lapply(pairwise_list, function(mats) {
      if (mat_name %in% names(mats)) mats[[mat_name]] else NULL
    })

    # Remove NULL entries
    non_null_matrices <- matrices[!sapply(matrices, is.null)]
    non_null_sizes <- batch_sizes[!sapply(matrices, is.null)]

    if (length(non_null_matrices) > 0) {
      # Create block diagonal matrix
      result[[mat_name]] <- create_block_diagonal(
        non_null_matrices,
        non_null_sizes
      )
    }
  }

  return(result)
}

#' @noRd
create_block_diagonal <- function(matrices, sizes) {
  if (length(matrices) == 0) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(0, 0)
    ))
  }

  if (length(matrices) == 1) {
    return(matrices[[1]])
  }

  # Use Matrix::bdiag for efficient block diagonal construction
  Matrix::bdiag(matrices)
}
