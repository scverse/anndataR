#' Generate a mock dataset
#'
#' Generate a mock synthetic dataset with different types of columns and layers.
#' This is primarily designed for use in tests, examples, vignettes and other
#' documentation but is also provided to users for creating reproducible
#' examples.
#'
#' @param n_obs Number of observations to generate
#' @param n_vars Number of variables to generate
#' @param x_type Type of matrix to generate for `X`
#' @param layer_types Types of matrices to generate for `layers`
#' @param obs_types Types of vectors to generate for `obs`
#' @param var_types Types of vectors to generate for `var`
#' @param obsm_types Types of matrices to generate for `obsm`
#' @param varm_types Types of matrices to generate for `varm`
#' @param obsp_types Types of matrices to generate for `obsp`
#' @param varp_types Types of matrices to generate for `varp`
#' @param uns_types Types of objects to generate for `uns`
#' @param example If `TRUE`, the types will be overridden with a small subset of
#'   types. This is useful for documentation.
#' @param format Object type to output, one of "list", "AnnData",
#'   "SingleCellExperiment", or "Seurat".
#'
#' @return For `generate_dataset()`, an object as defined by `output` containing
#'   the generated dataset
#' @export
#'
#' @details
#' To generate no data for a given slot, set the matching type argument to
#' `NULL` or an empty vector, e.g. `obs_types = c()` will generate an empty
#' `obs` data frame.
#'
#' When generating `SingleCellExperiment` or `Seurat` objects, only some of the
#' generated slots will be included in the output object. To generate a more
#' complete object, use `format = "AnnData"` followed by
#' `adata$as_SingleCellExperiment()` or `adata$as_Seurat()`.
#'
#' @examples
#' # Generate all types as a list
#' dummy <- generate_dataset()
#'
#' # Generate the example types
#' dummy_example <- generate_dataset(example = TRUE)
#'
#' # Generate an AnnData
#' dummy_anndata <- generate_dataset(format = "AnnData", example = TRUE)
#'
#' # Generate a SingleCellExperiment
#' if (rlang::is_installed("SingleCellExperiment")) {
#'   dummy_sce <- generate_dataset(format = "SingleCellExperiment", example = TRUE)
#' }
#'
#' # Generate a Seurat object
#' if (rlang::is_installed("SeuratObject")) {
#'   dummy_seurat <- generate_dataset(format = "Seurat", example = TRUE)
#' }
#'
generate_dataset <- function(
  n_obs = 10L,
  n_vars = 20L,
  x_type = "numeric_matrix",
  layer_types = get_generator_types(slot = "layers"),
  obs_types = get_generator_types(slot = "obs"),
  var_types = get_generator_types(slot = "var"),
  obsm_types = get_generator_types(slot = "obsm"),
  varm_types = get_generator_types(slot = "varm"),
  obsp_types = get_generator_types(slot = "obsp"),
  varp_types = get_generator_types(slot = "varp"),
  uns_types = get_generator_types(slot = "uns"),
  example = FALSE,
  format = c("list", "AnnData", "SingleCellExperiment", "Seurat")
) {
  format <- match.arg(format)

  if (example) {
    example_generator_types <- get_generator_types(example = TRUE)

    x_type <- example_generator_types$X
    layer_types <- example_generator_types$layers
    obs_types <- example_generator_types$obs
    var_types <- example_generator_types$var
    obsm_types <- example_generator_types$obsm
    varm_types <- example_generator_types$varm
    obsp_types <- example_generator_types$obsp
    varp_types <- example_generator_types$varp
    uns_types <- example_generator_types$uns
  }

  if (!rlang::is_empty(x_type) && length(x_type) != 1) {
    cli_abort("{.arg x_type} must be a single type")
  }

  for (slot in .anndata_slots) {
    types_arg <- if (slot == "X") {
      "x_type"
    } else if (slot == "layers") {
      "layer_types"
    } else {
      paste0(slot, "_types")
    }
    types <- get(types_arg)
    if (!all(types %in% .generator_types[[slot]])) {
      invalid_types <- types[!types %in% .generator_types[[slot]]] # nolint: object_use_linter
      cli_abort(c(
        "Some {.arg {types_arg}} types are not valid: {.val {invalid_types}}.",
        "i" = "Valid types are: {.val {.generator_types[[slot]]}}"
      ))
    }
  }

  dataset_list <- .generate_dataset_as_list(
    n_obs = n_obs,
    n_vars = n_vars,
    x_type = x_type,
    layer_types = layer_types,
    obs_types = obs_types,
    var_types = var_types,
    obsm_types = obsm_types,
    varm_types = varm_types,
    obsp_types = obsp_types,
    varp_types = varp_types,
    uns_types = uns_types
  )

  conversion_fun <- switch(
    format,
    "list" = identity,
    "SingleCellExperiment" = .generate_dataset_as_sce,
    "Seurat" = .generate_dataset_as_seurat,
    "AnnData" = .generate_dataset_as_anndata
  )

  conversion_fun(dataset_list)
}

#' Generate a dummy dataset as a list
#'
#' @inheritParams generate_dataset
#'
#' @return A list with the generated dataset
#'
#' @noRd
.generate_dataset_as_list <- function(
  n_obs = 10L,
  n_vars = 20L,
  x_type = get_generator_types(slot = "X")[1],
  layer_types = get_generator_types(slot = "layers"),
  obs_types = get_generator_types(slot = "obs"),
  var_types = get_generator_types(slot = "var"),
  obsm_types = get_generator_types(slot = "obsm"),
  varm_types = get_generator_types(slot = "varm"),
  obsp_types = get_generator_types(slot = "obsp"),
  varp_types = get_generator_types(slot = "varp"),
  uns_types = get_generator_types(slot = "uns")
) {
  # generate X
  X <- generate_matrix(n_obs, n_vars, x_type)

  # generate layers
  layers <- lapply(layer_types, generate_matrix, n_obs = n_obs, n_vars = n_vars)
  names(layers) <- layer_types

  # generate obs
  obs <- generate_dataframe(n_obs, obs_types)

  # generate var
  var <- generate_dataframe(n_vars, var_types)

  # generate obs_names
  obs_names <- paste0("cell", seq_len(n_obs))
  rownames(obs) <- obs_names

  # generate var_names
  var_names <- paste0("gene", seq_len(n_vars))
  rownames(var) <- var_names

  # generate obsm
  obsm <- lapply(obsm_types, function(obsm_type) {
    if (obsm_type %in% names(vector_generators)) {
      generate_dataframe(n_obs, obsm_type)
    } else {
      # generate n_obs vars to stay aligned with dummy-anndata, see scverse/anndataR#286
      generate_matrix(n_obs, n_vars = n_obs, obsm_type)
    }
  })
  names(obsm) <- obsm_types

  # generate varm
  varm <- lapply(varm_types, function(varm_type) {
    if (varm_type %in% names(vector_generators)) {
      generate_dataframe(n_vars, varm_type)
    } else {
      # generate n_vars vars to stay aligned with dummy-anndata, see scverse/anndataR#286
      generate_matrix(n_vars, n_vars = n_vars, varm_type)
    }
  })
  names(varm) <- varm_types

  # generate obsp
  obsp <- lapply(obsp_types, generate_matrix, n_obs = n_obs, n_vars = n_obs)
  names(obsp) <- obsp_types

  # generate varp
  varp <- lapply(varp_types, generate_matrix, n_obs = n_vars, n_vars = n_vars)
  names(varp) <- varp_types

  # generate uns by combining other classes
  uns <- lapply(uns_types, function(uns_type) {
    if (uns_type == "list") {
      # TODO: add multiple data types here
      list(
        1L,
        1,
        "a",
        factor("a"),
        TRUE,
        list(nested = list(list = c(1, 2, 3)))
      )
    } else if (uns_type == "scalar_character_with_nas") {
      NA_character_
    } else if (uns_type == "scalar_integer_with_nas") {
      NA_integer_
    } else if (uns_type == "scalar_factor_with_nas") {
      factor(NA_character_)
    } else if (uns_type == "scalar_factor_ordered_with_nas") {
      factor(NA_character_, ordered = TRUE)
    } else if (uns_type == "scalar_logical_with_nas") {
      NA_real_
    } else if (grepl("scalar_", uns_type)) {
      generate_vector(1L, gsub("scalar_", "", uns_type))
    } else if (grepl("vec_", uns_type)) {
      generate_vector(10L, gsub("vec_", "", uns_type))
    } else if (grepl("df_", uns_type)) {
      generate_dataframe(10L, gsub("df_", "", uns_type))
    } else if (grepl("mat_", uns_type)) {
      generate_matrix(10L, 10L, gsub("mat_", "", uns_type))
    } else {
      cli_abort("Unknown {.field uns} type: {.val {uns_type}}")
    }
  })
  names(uns) <- uns_types

  # return list
  list(
    X = X,
    obs = obs,
    obsm = obsm,
    obsp = obsp,
    var = var,
    varm = varm,
    varp = varp,
    layers = layers,
    uns = uns
  )
}

#' Convert a dummy dataset to a SingleCellExperiment object
#'
#' @param dataset_list Output of `.generate_dataset_as_list()`
#'
#' @return SingleCellExperiment containing the generated data
#'
#' @noRd
.generate_dataset_as_sce <- function(dataset_list) {
  check_requires(
    "Creating a SingleCellExperiment",
    "SingleCellExperiment",
    "Bioc"
  )

  assays_list <- c(
    list(X = dataset_list$X),
    dataset_list$layers
  )
  assays_list <- lapply(assays_list, Matrix::t)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    rowData = dataset_list$var,
    colData = dataset_list$obs
  )
  colnames(sce) <- dataset_list$obs_names
  rownames(sce) <- dataset_list$var_names

  # TODO: add obsm, varm, obsp, varp, uns?

  sce
}

#' Convert a dummy dataset to a Seurat object
#'
#' @param dataset_list Output of `.generate_dataset_as_list()`
#'
#' @return Seurat containing the generated data
#'
#' @noRd
.generate_dataset_as_seurat <- function(dataset_list) {
  check_requires("Creating a SeuratObject", "SeuratObject")

  if (!is.null(dataset_list$X)) {
    X <- t(dataset_list$X)
  } else if (!is.null(dataset_list$layers)) {
    X <- t(dataset_list$layers[[1]])
  } else {
    cli_abort(
      "Creating a {.cls Seurat} requires {.arg x_type} or {.arg layer_types} to be set"
    )
  }

  X <- t(dataset_list$X)
  colnames(X) <- dataset_list$obs_names
  rownames(X) <- dataset_list$var_names

  seurat <- SeuratObject::CreateSeuratObject(X)

  for (layer in names(dataset_list$layers)) {
    layer_data <- Matrix::t(dataset_list$layers[[layer]])
    colnames(layer_data) <- dataset_list$obs_names
    rownames(layer_data) <- dataset_list$var_names
    seurat <- SeuratObject::SetAssayData(seurat, layer, layer_data)
  }

  seurat <- SeuratObject::AddMetaData(seurat, dataset_list$obs)

  # TODO: add obsm, varm, obsp, varp, uns?

  seurat
}

#' Convert a dummy dataset to an AnnData object
#'
#' @param list Output of `.generate_dataset_as_list()`
#'
#' @return SingleCellExperiment containing the generated data
#'
#' @noRd
# nolint start: object_name_linter object_length_linter
.generate_dataset_as_anndata <- function(dataset_list) {
  # nolint end: object_name_linter object_length_linter
  AnnData(
    X = dataset_list$X,
    obs = dataset_list$obs,
    obsm = dataset_list$obsm,
    obsp = dataset_list$obsp,
    var = dataset_list$var,
    varm = dataset_list$varm,
    varp = dataset_list$varp,
    layers = dataset_list$layers,
    uns = dataset_list$uns
  )
}
