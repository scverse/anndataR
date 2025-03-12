#' Convert a SingleCellExperiment object to an AnnData object
#'
#' `from_SingleCellExperiment()` converts a
#'   SingleCellExperiment to an AnnData object.
#'
#' @param sce An object inheriting from SingleCellExperiment.
#'
#' @param output_class Name of the AnnData class. Must be one of `"HDF5AnnData"`
#' or `"InMemoryAnnData"`.
#'
#' @param x_mapping Name of the assay in `sce` to use as the `X` matrix in the AnnData object.
#' @param layers_mapping A named list mapping `assay` names in `sce` to `layers` in the created AnnData object.
#' The names of the list should be the names of the `layers` in the resulting AnnData object, and the values should be
#' the names of the `assays` in the `sce` object.
#' @param obs_mapping A named list mapping `colData` in `sce` to `obs` in the created AnnData object.
#' The names of the list should be the names of the `obs` columns in the resulting AnnData object. The values of the
#' list should be the names of the `colData` columns in `sce`.
#' @param var_mapping A named list mapping `rowData` in `sce` to `var` in the created AnnData object.
#' The names of the list should be the names of the `var` columns in the resulting AnnData object. The values of the
#' list should be the names of the `rowData` columns in `sce`.
#' @param obsm_mapping A named list mapping `reducedDim` in `sce` to `obsm` in the created AnnData object.
#' The names of the list should be the names of the `obsm` in the resulting AnnData object. The values of the list
#' should be a named list with as key the name of the `obsm` slot in the resulting AnnData object, and as value a list
#' with the following elements
#' - `reducedDim`
#' - the name of the `reducedDim` in `sce`
#' @param varm_mapping A named list mapping `reducedDim` in `sce` to `varm` in the created AnnData object.
#' The names of the list should be the names of the `varm` in the resulting AnnData object. The values of the list
#' should be a named list with as key the name of the `varm` slot in the resulting AnnData object, and as value a
#' list with the following elements
#' - `reducedDim`
#' - the name of the `reducedDim` in `sce`, that is `LinearEmbeddingMatrix` of which you want the featureLoadings to
#' end up in the `varm` slot
#' @param obsp_mapping A named list mapping `colPairs` in `sce` to `obsp` in the created AnnData object.
#' The names of the list should be the names of the `obsp` in the resulting AnnData object. The values of the list
#' should be the names of the `colPairs` in `sce`.
#' @param varp_mapping A named list mapping `rowPairs` in `sce` to `varp` in the created AnnData object.
#' The names of the list should be the names of the `varp` in the resulting AnnData object. The values of the list
#' should be the names of the `rowPairs` in `sce`.
#' @param uns_mapping A named list mapping `metadata` in `sce` to `uns` in the created AnnData object.
#' The names of the list should be the names of the `uns` in the resulting AnnData object. The values of the list
#' should be the names of the `metadata` in `sce`.
#' @param ... Additional arguments to pass to the generator function.
#'
#' @return `from_SingleCellExperiment()` returns an AnnData object
#'   (e.g., InMemoryAnnData) representing the content of `sce`.
#'
#' @examples
#' ## construct an AnnData object from a SingleCellExperiment
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(
#'   assays = list(counts = matrix(1:5, 5L, 3L)),
#'   colData = DataFrame(cell = 1:3),
#'   rowData = DataFrame(gene = 1:5)
#' )
#' from_SingleCellExperiment(sce, "InMemory")
#'
#' @export
# nolint start: object_name_linter
from_SingleCellExperiment <- function(
  # nolint end: object_name_linter
  sce,
  output_class = c("InMemory", "HDF5AnnData"),
  x_mapping = NULL,
  layers_mapping = NULL,
  obs_mapping = NULL,
  var_mapping = NULL,
  obsm_mapping = NULL,
  varm_mapping = NULL,
  obsp_mapping = NULL,
  varp_mapping = NULL,
  uns_mapping = NULL,
  ...
) {
  check_requires(
    "Converting SingleCellExperiment to AnnData",
    "SingleCellExperiment",
    "Bioc"
  )

  output_class <- match.arg(output_class)

  if (!(inherits(sce, "SingleCellExperiment"))) {
    cli_abort(
      "{.arg sce} must be a {.cls SingleCellExperiment} but has class {.cls {sce}}"
    )
  }

  # For any mappings that are not set, using the guessing function
  layers_mapping <- layers_mapping %||% .from_SCE_guess_layers(sce, x_mapping)
  obs_mapping <- obs_mapping %||%
    .from_SCE_guess_all(
      sce,
      SingleCellExperiment::colData
    )
  var_mapping <- var_mapping %||%
    .from_SCE_guess_all(
      sce,
      SingleCellExperiment::rowData
    )
  obsm_mapping <- obsm_mapping %||% .from_SCE_guess_obsm(sce)
  varm_mapping <- varm_mapping %||% .from_SCE_guess_varm(sce)
  obsp_mapping <- obsp_mapping %||%
    .from_SCE_guess_obspvarp(
      sce,
      SingleCellExperiment::colPairs
    )
  varp_mapping <- varp_mapping %||%
    .from_SCE_guess_obspvarp(
      sce,
      SingleCellExperiment::rowPairs
    )
  uns_mapping <- uns_mapping %||% .from_SCE_guess_all(sce, S4Vectors::metadata)

  # Create list to store converted objects
  # Also contains additional arguments passed to the generator function
  adata_list <- list(...)
  # Fill in slots in the list
  # trackstatus: class=SingleCellExperiment, feature=set_obs, status=wip
  adata_list$obs <- .from_SCE_process_obsvar(
    sce,
    obs_mapping,
    SingleCellExperiment::colData
  )
  # trackstatus: class=SingleCellExperiment, feature=set_var, status=wip
  adata_list$var <- .from_SCE_process_obsvar(
    sce,
    var_mapping,
    SingleCellExperiment::rowData
  )
  # trackstatus: class=SingleCellExperiment, feature=set_X, status=wip
  if (!is.null(x_mapping)) {
    adata_list$X <- .from_SCE_convert(
      SummarizedExperiment::assay(sce, x_mapping, withDimnames = FALSE)
    )
  }
  # trackstatus: class=SingleCellExperiment, feature=set_layers, status=wip
  adata_list$layers <- .from_SCE_process_simple_mapping(
    sce,
    layers_mapping,
    SummarizedExperiment::assays
  )
  # trackstatus: class=SingleCellExperiment, feature=set_obsm, status=wip
  adata_list$obsm <- .from_SCE_process_obsm(sce, obsm_mapping)
  # trackstatus: class=SingleCellExperiment, feature=set_varm, status=wip
  adata_list$varm <- .from_SCE_process_varm(sce, varm_mapping)
  # trackstatus: class=SingleCellExperiment, feature=set_obsp, status=wip
  adata_list$obsp <- .from_SCE_process_pairs(
    sce,
    obsp_mapping,
    SingleCellExperiment::colPairs,
    asSparse = TRUE
  )
  # trackstatus: class=SingleCellExperiment, feature=set_varp, status=wip
  adata_list$varp <- .from_SCE_process_pairs(
    sce,
    varp_mapping,
    SingleCellExperiment::rowPairs,
    asSparse = TRUE
  )
  # trackstatus: class=SingleCellExperiment, feature=set_uns, status=wip
  adata_list$uns <- .from_SCE_process_simple_mapping(
    sce,
    uns_mapping,
    S4Vectors::metadata,
    convert = FALSE
  )

  # Get a generator to create a new AnnData object
  generator <- get_anndata_constructor(output_class)

  tryCatch(
    {
      do.call(generator$new, adata_list)
    },
    error = function(e) {
      if (output_class == "HDF5AnnData") {
        on.exit(cleanup_HDF5AnnData(adata_list$file))
      }
      cli_abort(
        c(
          "Failed to create a {.cls {output_class}} with error: {conditionMessage(e)}",
          "i" = "Error call: {.code {capture.output(print(conditionCall(e)))}}"
        ),
        call = rlang::caller_env(4)
      )
    }
  )
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_all <- function(sce, slot) {
  # nolint end: object_length_linter object_name_linter
  mapping <- names(slot(sce))
  names(mapping) <- names(slot(sce))

  mapping
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_layers <- function(sce, x_mapping) {
  # nolint end: object_length_linter object_name_linter
  layers_mapping <- list()

  for (assay_name in names(SummarizedExperiment::assays(sce))) {
    if (is.null(x_mapping) || assay_name != x_mapping) {
      layers_mapping[[assay_name]] <- assay_name
    }
  }

  layers_mapping
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_obsm <- function(sce) {
  # nolint end: object_length_linter object_name_linter
  if (!inherits(sce, "SingleCellExperiment")) {
    return(list())
  }
  obsm_mapping <- list()

  for (reduction_name in names(SingleCellExperiment::reducedDims(sce))) {
    dest_name <- paste0("X_", reduction_name)
    obsm_mapping[[dest_name]] <- c("reducedDim", reduction_name)
  }

  obsm_mapping
}

# nolint start: object_length_linter object_name_linter
.from_SCE_guess_varm <- function(sce) {
  # nolint end: object_length_linter object_name_linter
  if (!inherits(sce, "SingleCellExperiment")) {
    return(list())
  }
  varm_mapping <- list()

  for (reduction_name in names(SingleCellExperiment::reducedDims(sce))) {
    reduction <- SingleCellExperiment::reducedDim(sce, reduction_name)
    if (inherits(reduction, "LinearEmbeddingMatrix")) {
      dest_name <- paste0(reduction_name, "s") # this is for PCA, should this be generalized?
      varm_mapping[[dest_name]] <- c("reducedDim", reduction_name)
    }
  }

  varm_mapping
}

# If sce is a SummarizedExperiment, return an empty mapping
# nolint start: object_length_linter object_name_linter
.from_SCE_guess_obspvarp <- function(sce, slot) {
  # nolint end: object_length_linter object_name_linter
  if (!inherits(sce, "SingleCellExperiment")) {
    return(list())
  }
  .from_SCE_guess_all(sce, slot)
}

# Convert BioConductor-specific objects to base R objects
# Convert matrices
# nolint start: object_length_linter object_name_linter
.from_SCE_convert <- function(object) {
  # nolint end: object_length_linter object_name_linter
  if (inherits(object, "DataFrame")) {
    as.data.frame(object)
  } else if (inherits(object, "SimpleList")) {
    as.list(object)
  } else if (inherits(object, "matrix") || inherits(object, "Matrix")) {
    m <- t(object)
    if (inherits(m, "denseMatrix")) {
      m <- as.matrix(m)
    }
    m
  } else {
    object
  }
}

# nolint start: object_length_linter object_name_linter
.from_SCE_process_obsvar <- function(sce, mapping, slot, convert = TRUE) {
  # nolint end: object_length_linter object_name_linter
  mapped <- .from_SCE_process_simple_mapping(sce, mapping, slot, convert)

  # Ensure that the obs & var retain rownames, even if there are no columns
  # in the DataFrame.
  if (inherits(mapped, "list") && length(mapped) == 0) {
    if (!is.null(rownames(slot(sce)))) {
      mapped <- as.data.frame(rownames(slot(sce)))
      rownames(mapped) <- rownames(slot(sce))
      mapped[, 1] <- NULL
    } else {
      # We will infer rownames
      inferred_rownames <- paste0("", seq_len(nrow(slot(sce))))
      mapped <- as.data.frame(inferred_rownames)
      rownames(mapped) <- inferred_rownames
      mapped[, 1] <- NULL
    }
  } else {
    mapped <- as.data.frame(mapped)
    rownames(mapped) <- rownames(slot(sce))
  }

  mapped
}
# nolint start: object_length_linter object_name_linter
.from_SCE_process_simple_mapping <- function(
  # nolint end: object_length_linter object_name_linter
  sce,
  mapping,
  slot,
  convert = TRUE
) {
  mapped <- NULL
  mapped <- lapply(seq_along(mapping), function(i) {
    element <- slot(sce)[[mapping[[i]]]]
    if (convert) {
      element <- .from_SCE_convert(element)
    }
    element
  })
  names(mapped) <- names(mapping)

  .from_SCE_convert(mapped)
}

# nolint start: object_length_linter object_name_linter
.from_SCE_process_pairs <- function(sce, mapping, slot, asSparse = TRUE) {
  # nolint end: object_length_linter object_name_linter
  if (!inherits(sce, "SingleCellExperiment")) {
    return(list())
  }
  pairs <- NULL
  # check if mapping contains all columns of slot
  if (length(setdiff(names(slot(sce)), names(mapping))) == 0) {
    pairs <- slot(sce, asSparse = asSparse)
  } else {
    pairs <- lapply(seq_along(mapping), function(i) {
      .from_SCE_convert(slot(sce, asSparse = asSparse)[[mapping[[i]]]])
    })
    names(pairs) <- names(mapping)
  }

  .from_SCE_convert(pairs)
}

# nolint start: object_name_linter
.from_SCE_process_obsm <- function(sce, obsm_mapping) {
  # nolint end: object_name_linter
  purrr::imap(obsm_mapping, function(.mapping, .idx) {
    if (!is.character(.mapping) || length(.mapping) != 2) {
      cli_abort(c(
        paste(
          "Each item in {.arg obsm_mapping} must be a {.cls character}",
          "vector of length 2"
        ),
        "i" = paste(
          "{.code obsm_mapping[[{.val { .idx }}]]} is",
          "{.obj_type_friendly { .mapping }}"
        )
      ))
    }

    slot <- .mapping[1]
    key <- .mapping[2]

    if (slot != "reducedDim") {
      cli_abort(c(
        "The first item in each {.arg obsm_mapping} must be {.val reducedDim}",
        "i" = "{.code obsm_mapping[[{.val { .idx }}]][1]}: {.val {slot}}"
      ))
    }

    reduction <- SingleCellExperiment::reducedDim(sce, key)
    if (inherits(reduction, "LinearEmbeddingMatrix")) {
      SingleCellExperiment::sampleFactors(reduction)
    } else {
      reduction
    }
  })
}

# nolint start: object_name_linter
.from_SCE_process_varm <- function(sce, varm_mapping) {
  # nolint end: object_name_linter
  purrr::imap(varm_mapping, function(.mapping, .idx) {
    if (!is.character(.mapping) || length(.mapping) != 2) {
      cli_abort(c(
        paste(
          "Each item in {.arg varm_mapping} must be a {.cls character}",
          "vector of length 2"
        ),
        "i" = paste(
          "{.code varm_mapping[[{.val { .idx }}]]} is",
          "{.obj_type_friendly { .mapping }}"
        )
      ))
    }

    slot <- .mapping[1]
    key <- .mapping[2]

    if (slot != "reducedDim") {
      cli_abort(c(
        "The first item in each {.arg varm_mapping} must be {.val reducedDim}",
        "i" = "{.code varm_mapping[[{.val { .idx }}]][1]}: {.val {slot}}"
      ))
    }

    reduction <- SingleCellExperiment::reducedDim(sce, key)
    if (!inherits(reduction, "LinearEmbeddingMatrix")) {
      cli_abort(paste(
        "{.code reducedDim(sce, {.val key}} must be a {.cls LinearEmbeddingMatrix}",
        "but has class {.cls {reduction}}"
      ))
    }

    SingleCellExperiment::featureLoadings(reduction)
  })
}
