#' Concatenate AnnData objects
#'
#' Combine multiple AnnData objects along specified axis with flexible
#' strategies for handling mismatched dimensions and metadata.
#'
#' Related files:
#' - concat_helpers.R: Helper functions for matrix and annotation concatenation
#' - concat_strategies.R: Merge strategies for handling non-aligned elements
#'
#' @param adatas List of AnnData objects to concatenate, or named list for batch labels
#' @param axis Which axis to concatenate along. Either "obs"/0 for observations (rows)
#'   or "var"/1 for variables (columns)
#' @param join How to align the other axis. "outer" takes the union of indices,
#'   "inner" takes the intersection
#' @param merge Strategy for elements not aligned to concatenation axis. One of:
#'   \itemize{
#'     \item "same": Keep elements that are identical across all objects
#'     \item "unique": Keep elements that appear in only one object (or are identical)
#'     \item "first": Keep the first occurrence of each element
#'     \item "only": Keep elements that appear in exactly one object
#'   }
#' @param uns_merge Strategy for merging .uns metadata (same options as merge)
#' @param label Column name to add batch information to obs/var. NULL for no label
#' @param keys Batch labels. Defaults to names(adatas) if named list, otherwise sequential integers
#' @param index_unique Separator to make row/column names unique. NULL keeps original names
#' @param fill_value Value for missing data when join="outer". Default: 0 for sparse, NA for dense
#' @param backend Backend to use for result ("memory" or "hdf5")
#'
#' @return A new AnnData object containing the concatenated data
#'
#' @details
#' This function provides flexible concatenation of AnnData objects similar to
#' Python's anndata.concat(). It handles:
#'
#' - Mismatched observation/variable names through join strategies
#' - Complex metadata merging through merge strategies
#' - Batch tracking through labeling
#' - Memory-efficient operations for large datasets
#'
#' @examples
#' \dontrun{
#' # Create example datasets
#' ad1 <- generate_dataset(n_obs = 100, n_vars = 50, format = "AnnData")
#' ad2 <- generate_dataset(n_obs = 80, n_vars = 50, format = "AnnData")
#' ad3 <- generate_dataset(n_obs = 120, n_vars = 50, format = "AnnData")
#'
#' # Basic row concatenation
#' combined <- concat(list(ad1, ad2, ad3), axis = "obs")
#'
#' # With batch tracking
#' combined <- concat(list(ctrl = ad1, treat = ad2), axis = "obs", label = "condition")
#'
#' # Inner join (intersection of variables)
#' combined <- concat(list(ad1, ad2), axis = "obs", join = "inner")
#' }
#'
#' @export
concat <- function(
  adatas,
  axis = c("obs", "var", 0, 1),
  join = c("outer", "inner"),
  merge = c("unique", "same", "first", "only"),
  uns_merge = merge,
  label = NULL,
  keys = NULL,
  index_unique = NULL,
  fill_value = NULL,
  backend = c("memory", "hdf5")
) {
  # Input validation
  if (is.null(adatas) || length(adatas) == 0) {
    cli::cli_abort("adatas must be a non-empty list")
  }

  # Handle named lists
  if (is.null(keys) && !is.null(names(adatas))) {
    keys <- names(adatas)
    adatas <- unname(adatas)
  }

  # Convert to list if single object
  if (inherits(adatas, "AbstractAnnData")) {
    adatas <- list(adatas)
  }

  # Validate all objects are AnnData
  if (!all(sapply(adatas, inherits, "AbstractAnnData"))) {
    cli::cli_abort("All objects in adatas must inherit from AbstractAnnData")
  }

  if (length(adatas) == 1) {
    cli::cli_warn("Only one object provided, returning copy")
    return(as_InMemoryAnnData(adatas[[1]]))
  }

  # Resolve arguments
  axis <- resolve_axis(axis)
  join <- match.arg(join)
  merge <- match.arg(merge)
  uns_merge <- match.arg(uns_merge, c("unique", "same", "first", "only"))
  backend <- match.arg(backend)

  if (is.null(keys)) {
    keys <- as.character(seq_along(adatas))
  } else if (length(keys) != length(adatas)) {
    cli::cli_abort("Length of keys must match length of adatas")
  }

  # Check for empty objects
  empty_mask <- sapply(adatas, function(ad) any(dim(ad) == 0))
  if (any(empty_mask)) {
    cli::cli_warn("Removing {sum(empty_mask)} empty object{?s}")
    adatas <- adatas[!empty_mask]
    keys <- keys[!empty_mask]

    if (length(adatas) == 0) {
      cli::cli_abort("No non-empty objects remaining")
    }
  }

  # Set default fill value
  if (is.null(fill_value)) {
    all_X <- lapply(adatas, function(ad) ad$X)
    fill_value <- default_fill_value(all_X)
  }

  # Perform concatenation
  result <- concat_impl(
    adatas = adatas,
    axis = axis,
    join = join,
    merge = merge,
    uns_merge = uns_merge,
    label = label,
    keys = keys,
    index_unique = index_unique,
    fill_value = fill_value,
    backend = backend
  )

  return(result)
}

#' Internal implementation of concat
#' @inheritParams concat
#' @keywords internal
#' @noRd
concat_impl <- function(
  adatas,
  axis,
  join,
  merge,
  uns_merge,
  label,
  keys,
  index_unique,
  fill_value,
  backend
) {
  n_obs <- sapply(adatas, function(ad) ad$n_obs())
  n_vars <- sapply(adatas, function(ad) ad$n_vars())

  # Determine concatenation and alignment axes
  if (axis == 0L) {
    # Concatenating observations (rows)
    concat_axis <- 0L
    align_axis <- 1L
    concat_indices_list <- lapply(adatas, function(ad) ad$obs_names)
    align_indices_list <- lapply(adatas, function(ad) ad$var_names)
  } else {
    # Concatenating variables (columns)
    concat_axis <- 1L
    align_axis <- 0L
    concat_indices_list <- lapply(adatas, function(ad) ad$var_names)
    align_indices_list <- lapply(adatas, function(ad) ad$obs_names)
  }

  # Create batch labels
  if (axis == 0L) {
    batch_sizes <- n_obs
  } else {
    batch_sizes <- n_vars
  }

  batch_labels <- rep(keys, batch_sizes)
  batch_labels <- factor(batch_labels, levels = keys)

  # Merge and create indices for alignment axis
  merged_align_indices <- merge_indices(align_indices_list, join = join)
  reindexers <- lapply(align_indices_list, function(idx) {
    gen_reindexer(merged_align_indices, idx)
  })

  # Create concatenated indices
  concat_indices <- unlist(concat_indices_list)
  if (!is.null(index_unique)) {
    # Make indices unique by appending batch keys
    concat_indices <- paste(
      concat_indices,
      rep(keys, batch_sizes),
      sep = index_unique
    )
  } else {
    # Check if indices are already unique, if not make them unique
    if (any(duplicated(concat_indices))) {
      concat_indices <- paste(concat_indices, rep(keys, batch_sizes), sep = "_")
    }
  }

  # Concatenate main matrix (X)
  X <- concat_X(adatas, reindexers, concat_axis, fill_value)

  # Concatenate observation/variable annotations
  if (axis == 0L) {
    # Concatenating observations - combine obs, merge var
    obs <- concat_annotations(
      lapply(adatas, function(ad) ad$obs),
      concat_indices,
      join = "outer" # obs always outer join along concat axis
    )

    if (!is.null(label)) {
      obs[[label]] <- batch_labels
    }

    var <- merge_annotations(
      lapply(adatas, function(ad) ad$var),
      merged_align_indices,
      strategy = merge
    )
  } else {
    # Concatenating variables - merge obs, combine var
    obs <- merge_annotations(
      lapply(adatas, function(ad) ad$obs),
      merged_align_indices,
      strategy = merge
    )

    var <- concat_annotations(
      lapply(adatas, function(ad) ad$var),
      concat_indices,
      join = "outer" # var always outer join along concat axis
    )

    if (!is.null(label)) {
      var[[label]] <- batch_labels
    }
  }

  # Handle layers
  layers <- concat_layers(
    adatas,
    reindexers,
    concat_axis,
    join,
    fill_value,
    merge
  )

  # Handle obsm/varm (observation/variable matrices)
  if (axis == 0L) {
    # Concatenating obs - combine obsm, merge varm
    obsm <- concat_matrices_dict(
      lapply(adatas, function(ad) ad$obsm),
      axis = 0L,
      join = "outer",
      fill_value = fill_value
    )

    varm <- merge_matrices_dict(
      lapply(adatas, function(ad) ad$varm),
      reindexers,
      axis = 0L,
      strategy = merge
    )
  } else {
    # Concatenating var - merge obsm, combine varm
    obsm <- merge_matrices_dict(
      lapply(adatas, function(ad) ad$obsm),
      reindexers,
      axis = 0L,
      strategy = merge
    )

    varm <- concat_matrices_dict(
      lapply(adatas, function(ad) ad$varm),
      axis = 0L,
      join = "outer",
      fill_value = fill_value
    )
  }

  # Handle obsp/varp (pairwise matrices) - these need special block diagonal handling
  obsp <- concat_pairwise_dict(
    lapply(adatas, function(ad) ad$obsp),
    batch_sizes,
    axis = if (axis == 0L) 0L else NULL,
    strategy = merge
  )

  varp <- concat_pairwise_dict(
    lapply(adatas, function(ad) ad$varp),
    batch_sizes,
    axis = if (axis == 1L) 1L else NULL,
    strategy = merge
  )

  # Handle uns (unstructured metadata)
  uns <- merge_nested(
    lapply(adatas, function(ad) ad$uns),
    strategy = uns_merge
  )

  # Create result object
  if (backend == "memory") {
    result <- InMemoryAnnData$new(
      X = X,
      obs = obs,
      var = var,
      layers = layers,
      obsm = obsm,
      varm = varm,
      obsp = obsp,
      varp = varp,
      uns = uns
    )
  } else {
    cli::cli_abort("HDF5 backend not yet implemented for concat")
  }

  return(result)
}

#' Resolve axis argument to integer
#' @param axis Axis specification
#' @return Integer axis (0 or 1)
#' @keywords internal
#' @noRd
resolve_axis <- function(axis) {
  if (is.character(axis)) {
    axis <- match.arg(axis, c("obs", "var"))
    if (axis == "obs") {
      return(0L)
    }
    if (axis == "var") return(1L)
  } else if (is.numeric(axis)) {
    axis <- as.integer(axis)
    if (axis %in% c(0L, 1L)) return(axis)
  }

  cli::cli_abort("axis must be 'obs', 'var', 0, or 1")
}

#' @rdname concat
#' @method rbind AbstractAnnData
#' @export
rbind.AbstractAnnData <- function(
  ...,
  join = "outer",
  merge = "unique",
  label = "batch",
  fill_value = NULL
) {
  adatas <- list(...)

  # Create keys from argument names if available
  arg_names <- ...names()
  if (!is.null(arg_names) && !all(arg_names == "")) {
    keys <- arg_names
    keys[keys == ""] <- paste0("X", which(keys == ""))
  } else {
    keys <- paste0("X", seq_along(adatas))
  }

  concat(
    adatas,
    axis = "obs",
    join = join,
    merge = merge,
    label = label,
    keys = keys,
    fill_value = fill_value
  )
}

#' @rdname concat
#' @method cbind AbstractAnnData
#' @export
cbind.AbstractAnnData <- function(
  ...,
  join = "outer",
  merge = "unique",
  label = "batch",
  fill_value = NULL
) {
  adatas <- list(...)

  # Create keys from argument names if available
  arg_names <- ...names()
  if (!is.null(arg_names) && !all(arg_names == "")) {
    keys <- arg_names
    keys[keys == ""] <- paste0("X", which(keys == ""))
  } else {
    keys <- paste0("X", seq_along(adatas))
  }

  concat(
    adatas,
    axis = "var",
    join = join,
    merge = merge,
    label = label,
    keys = keys,
    fill_value = fill_value
  )
}
