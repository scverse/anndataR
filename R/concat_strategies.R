#' Concatenation strategies for anndataR
#'
#' Functions that implement different strategies for merging elements
#' that are not aligned to the concatenation axis when concatenating AnnData objects.
#'
#' @name concat_strategies
#' @keywords internal

#' @noRd
apply_merge_strategy <- function(mappings, strategy) {
  if (is.character(strategy)) {
    strategy <- get_merge_strategy(strategy)
  }

  # Get union of all keys
  all_keys <- unique(unlist(lapply(mappings, names)))

  result <- list()
  for (key in all_keys) {
    # Get values for this key from all mappings
    values <- lapply(mappings, function(m) {
      if (key %in% names(m)) m[[key]] else NULL
    })

    # Apply strategy to determine which value to keep
    merged_value <- strategy(values)
    if (!is.null(merged_value)) {
      result[[key]] <- merged_value
    }
  }

  return(result)
}

#' @noRd
get_merge_strategy <- function(strategy_name) {
  strategies <- list(
    "same" = merge_same,
    "unique" = merge_unique,
    "first" = merge_first,
    "only" = merge_only
  )

  if (!strategy_name %in% names(strategies)) {
    cli::cli_abort(
      "Unknown merge strategy: {strategy_name}. Must be one of: {names(strategies)}"
    )
  }

  return(strategies[[strategy_name]])
}

#' @noRd
merge_same <- function(values) {
  # Remove NULL values
  non_null_values <- values[!sapply(values, is.null)]

  if (length(non_null_values) == 0) {
    return(NULL)
  }

  # Check if all non-NULL values are identical
  first_value <- non_null_values[[1]]
  all_same <- all(sapply(non_null_values, function(x) {
    if (is.data.frame(x) && is.data.frame(first_value)) {
      return(identical(x, first_value))
    } else if (is.matrix(x) && is.matrix(first_value)) {
      return(identical(x, first_value))
    } else {
      return(identical(x, first_value))
    }
  }))

  if (all_same) {
    return(first_value)
  } else {
    return(NULL)
  }
}

#' @noRd
merge_unique <- function(values) {
  # Remove NULL values
  non_null_values <- values[!sapply(values, is.null)]

  if (length(non_null_values) == 0) {
    return(NULL)
  } else if (length(non_null_values) == 1) {
    return(non_null_values[[1]])
  } else {
    # Check if all values are identical
    return(merge_same(values))
  }
}

#' @noRd
merge_first <- function(values) {
  for (value in values) {
    if (!is.null(value)) {
      return(value)
    }
  }
  return(NULL)
}

#' @noRd
merge_only <- function(values) {
  # Remove NULL values
  non_null_values <- values[!sapply(values, is.null)]

  if (length(non_null_values) == 1) {
    return(non_null_values[[1]])
  } else {
    return(NULL)
  }
}

#' @noRd
merge_nested <- function(mappings, strategy) {
  if (is.character(strategy)) {
    strategy <- get_merge_strategy(strategy)
  }

  # Get all top-level keys
  all_keys <- unique(unlist(lapply(mappings, names)))

  result <- list()
  for (key in all_keys) {
    # Get values for this key from all mappings
    values <- lapply(mappings, function(m) {
      if (key %in% names(m)) m[[key]] else NULL
    })

    # Check if all non-NULL values are themselves mappings (lists)
    non_null_values <- values[!sapply(values, is.null)]
    if (
      length(non_null_values) > 0 &&
        all(sapply(non_null_values, function(x) {
          is.list(x) && !is.data.frame(x)
        }))
    ) {
      # Recursively merge nested mappings
      nested_result <- merge_nested(non_null_values, strategy)
      if (length(nested_result) > 0) {
        result[[key]] <- nested_result
      }
    } else {
      # Apply strategy to leaf values
      merged_value <- strategy(values)
      if (!is.null(merged_value)) {
        result[[key]] <- merged_value
      }
    }
  }

  return(result)
}

#' @noRd
merge_uns_with_batches <- function(uns_list, batch_keys, join_str = "_batch_") {
  if (length(uns_list) != length(batch_keys)) {
    cli::cli_abort("Length of uns_list and batch_keys must be equal")
  }

  # First apply normal merge strategy
  merged <- apply_merge_strategy(uns_list, "unique")

  # Then add batch-specific keys for elements that didn't merge
  all_keys <- unique(unlist(lapply(uns_list, names)))

  for (key in all_keys) {
    if (!key %in% names(merged)) {
      # This key didn't survive merging, add batch-specific versions
      for (i in seq_along(uns_list)) {
        if (key %in% names(uns_list[[i]])) {
          batch_key <- paste0(key, join_str, batch_keys[i])
          merged[[batch_key]] <- uns_list[[i]][[key]]
        }
      }
    }
  }

  return(merged)
}
