# Notes to self:
# https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#2_Creating_SingleCellExperiment_instances
# https://bioconductor.org/packages/3.20/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#column-sample-data
# https://github.com/ivirshup/sc-interchange/issues
# https://github.com/ivirshup/sc-interchange/issues/2
# https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#3_Adding_low-dimensional_representations

#' Convert an AnnData object to a SingleCellExperiment object
#'
#' `to_SingleCellExperiment()` converts an AnnData object
#'   to a SingleCellExperiment object.
#'
#' @param adata an AnnData object, e.g., InMemoryAnnData
#'
#' @param assays_mapping A named list mapping `layers` in `adata` to
#'   `assay` names in the created SingleCellExperiment object.
#' @param col_data_mapping a named list mapping `obs` in `adata` to
#'  `colData` in the created SingleCellExperiment object
#' @param row_data_mapping a named list mapping `var` names in `adata` to
#'   `rowData` in the created SingleCellExperiment object
#' @param reduction_mapping a named list mapping reduction names in `adata` to
#' reduction names in the created SingleCellExperiment object
#' @param colPairs_mapping a named list mapping obsp names in `adata` to
#' colPairs names in the created SingleCellExperiment object
#' @param rowPairs_mapping a named list mapping varp names in `adata` to
#' rowPairs names in the created SingleCellExperiment object
#' @param metadata_mapping a named list mapping uns names in `adata` to
#' metadata names in the created SingleCellExperiment object
#'
#' @return `to_SingleCellExperiment()` returns a SingleCellExperiment
#'   representing the content of `object`.
#'
#' @noRd
#' @examples
#' if (interactive()) {
#'   ## useful when interacting with the SingleCellExperiment !
#'   library(SingleCellExperiment)
#' }
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#'
#' ## construct a SingleCellExperiment from an AnnData object
#' sce <- to_SingleCellExperiment(ad)
#' sce
to_SingleCellExperiment <- function(
    adata,
    assays_mapping = NULL,
    col_data_mapping = NULL,
    row_data_mapping = NULL,
    reduction_mapping = NULL,
    colPairs_mapping = NULL,
    rowPairs_mapping = NULL,
    metadata_mapping = NULL) {

  check_requires("Converting AnnData to SingleCellExperiment", "SingleCellExperiment", "Bioc")

  stopifnot(inherits(adata, "AbstractAnnData"))

  # guess mappings if not provided
  if (is.null(assays_mapping)) {
    assays_mapping <- to_SCE_guess_assays(adata)
  }
  if (is.null(col_data_mapping)) {
    col_data_mapping <- to_SCE_guess_all(adata, "obs")
  }
  if(is.null(row_data_mapping)) {
    row_data_mapping <- to_SCE_guess_all(adata, "var")
  }
  if (is.null(reduction_mapping)) {
    reduction_mapping <- to_SCE_guess_reduction(adata)
  }
  if (is.null(colPairs_mapping)) {
    colPairs_mapping <- to_SCE_guess_all(adata, "obsp")
  }
  if (is.null(rowPairs_mapping)) {
    rowPairs_mapping <- to_SCE_guess_all(adata, "varp")
  }
  if (is.null(metadata_mapping)) {
    metadata_mapping <- to_SCE_guess_all(adata, "uns")
  }

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  sce_assays <- vector("list", length(assays_mapping))
  names(sce_assays) <- names(assays_mapping)
  for (i in seq_along(assays_mapping)){
    from <- assays_mapping[[i]]
    to <- names(assays_mapping)[[i]]
    if (from != "X") {
      sce_assays[[to]] <- t(adata$layers[[from]])
    } else {
      sce_assays[[to]] <- t(adata$X)
    }
  }

  # construct colData
  # FIXME: probably better way to make a dataframe from a list of vectors
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- .to_SCE_process_simple_mapping(adata, col_data_mapping, "obs")

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- .to_SCE_process_simple_mapping(adata, row_data_mapping, "var")

  # construct reducedDims
  # trackstatus: class=SingleCellExperiment, feature=get_reductions, status=wip
  reduceddims <- vector("list", length(reduction_mapping))
  names(reduceddims) <- names(reduction_mapping)
  for (i in seq_along(reduction_mapping)) {

    name <- names(reduction_mapping)[[i]]
    reduction <- reduction_mapping[[i]]

    obsm_key <- reduction$obsm
    varm_key <- reduction$varm
    uns_key <- reduction$uns

    reduceddims[[name]] <- .to_SingleCellExperiment_process_reduction(adata, name, obsm_key, varm_key, uns_key)
  }

  # construct colPairs
  # trackstatus: class=SingleCellExperiment, feature=get_obsp, status=done
  col_pairs <- .to_SCE_process_simple_mapping(adata, colPairs_mapping, "obsp")

  # construct rowPairs
  # trackstatus: class=SingleCellExperiment, feature=get_varp, status=done
  row_pairs <- .to_SCE_process_simple_mapping(adata, rowPairs_mapping, "varp")

  # construct metadata
  # trackstatus: class=SingleCellExperiment, feature=get_uns, status=done
  metadata <- .to_SCE_process_simple_mapping(adata, metadata_mapping, "uns")

  arguments <- list(
    assays = sce_assays,
    reducedDims = reduceddims,
    colPairs = col_pairs,
    rowPairs = row_pairs,
    metadata = metadata,
    checkDimnames = TRUE
  )
  # add col_data if not empty list
  if (length(col_data) > 0) {
    arguments$colData <- as(col_data, "DataFrame")
  }
  # add row_data if not empty list
  if (length(row_data) > 0) {
    arguments$rowData <- as(row_data, "DataFrame")
  }

  # construct output object
  sce <- do.call(SingleCellExperiment::SingleCellExperiment, arguments)

  sce
}

to_SCE_guess_assays <- function(adata) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  layers <- list()

  if (!is.null(adata$X)) {
    layer_name_for_x <- 
      if (!"counts" %in% names(adata$layers)) { # could expand checks, to check if integers
        "counts"
      } else {
        "data"
      }
    layers[[layer_name_for_x]] <- "X"
  }

  for (layer_name in names(adata$layers)) {
    layers[[layer_name]] <- layer_name
  }

  layers
}

to_SCE_guess_all <- function(adata, slot) {
  if (!inherits(adata, "AbstractAnnData")) {
    stop("adata must be an object inheriting from AbstractAnnData")
  }

  mapping <- names(adata[[slot]])
  names(mapping) <- names(adata[[slot]])

  mapping
}

to_SCE_guess_reduction <- function(adata) {
  to_Seurat_guess_reductions(adata)
}

.to_SCE_process_simple_mapping <- function(adata, mapping, slot) {
  # check if mapping contains all columns of slot
  if (length(setdiff(names(adata[[slot]]), names(mapping))) == 0) {
    adata[[slot]]
  } else {
    mapped <- lapply(seq_along(mapping), function(i) {
      adata[[slot]][[mapping[[i]]]]
    })
    names(mapped) <- names(mapping)
    mapped
  }
}

.to_SingleCellExperiment_process_reduction <- function(adata, key, obsm_key, varm_key, uns_key) {
  # check arguments

  embedding <- adata$obsm[[obsm_key]]

  if (is.null(embedding)) {
    stop(paste0("No embedding found for key '", obsm_key, "' in adata$obsm"))
  }

  rownames(embedding) <- adata$obs_names

  if (! is.null(varm_key) && varm_key %in% names(adata$varm)) {
    loadings <- adata$varm[[varm_key]]
    rownames(loadings) <- colnames(embedding)
    
    metadata <- list()
    if (! is.null(uns_key) && uns_key %in% names(adata$uns)) {
      metadata <- adata$uns[[uns_key]]
    }

    lem <- LinearEmbeddingMatrix(
      sampleFactors = embedding,
      featureLoadings = loadings,
      metadata = metadata
    )
    rownames(lem) <- rownames(embedding)
    colnames(lem) <- colnames(loadings)
    lem
  } else {
    embedding
  }
}

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
#' @param ... Additional arguments passed to the generator function.
#' See the "Details" section for more information on which parameters
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
    ...) {

  check_requires("Converting SingleCellExperiment to AnnData", "SingleCellExperiment", "Bioc")

  output_class <- match.arg(output_class)

  stopifnot(inherits(sce, "SingleCellExperiment"))

  if (is.null(layers_mapping)) {
    layers_mapping <- .from_SCE_guess_layers(sce, x_mapping)
  }
  if (is.null(obs_mapping)) {
    obs_mapping <- .from_SCE_guess_all(sce, colData)
  }
  if (is.null(var_mapping)) {
    var_mapping <- .from_SCE_guess_all(sce, rowData)
  }
  if (is.null(obsm_mapping)) {
    obsm_mapping <- .from_SCE_guess_obsm(sce)
  }
  if (is.null(varm_mapping)) {
    varm_mapping <- list()
  }
  if (is.null(obsp_mapping)) {
    obsp_mapping <- .from_SCE_guess_all(sce, colPairs)
  }
  if (is.null(varp_mapping)) {
    varp_mapping <- .from_SCE_guess_all(sce, rowPairs)
  }
  if (is.null(uns_mapping)) {
    uns_mapping <- .from_SCE_guess_all(sce, metadata)
  }

  # fetch generator
  generator <- get_anndata_constructor(output_class)

  tryCatch(
    {
      # get obs
      # trackstatus: class=SingleCellExperiment, feature=set_obs, status=wip
      obs <- .from_SCE_process_simple_mapping(sce, obs_mapping, colData)

      # get var
      # trackstatus: class=SingleCellExperiment, feature=set_var, status=wip
      var <- .from_SCE_process_simple_mapping(sce, var_mapping, rowData)

      adata <- generator$new(
        obs = obs,
        var = var
      )

      # fetch X
      # trackstatus: class=SingleCellExperiment, feature=set_X, status=wip
      X <- NULL
      if (!is.null(x_mapping)) {
        X <- .from_SCE_convert(SummarizedExperiment::assay(sce, x_mapping, withDimnames = FALSE))
      }
      adata$X <- X

      # fetch layers
      # trackstatus: class=SingleCellExperiment, feature=set_layers, status=wip
      # for (i in seq_along(layers_mapping)) {
      #   layer <- layers_mapping[[i]]
      #   layer_name <- names(layers_mapping)[[i]]

      #   adata$layers[[layer_name]] <- .to_anndata_matrix(assay(sce, layer))
      # }

      adata$layers <- .from_SCE_process_simple_mapping(sce, layers_mapping, assays)

      # trackstatus: class=SingleCellExperiment, feature=set_obsm, status=wip
      for (i in seq_along(obsm_mapping)) {
        obsm <- obsm_mapping[[i]]
        obsm_name <- names(obsm_mapping)[[i]]

        if (!is.character(obsm) || length(obsm) != 2) {
          stop("each obsm_mapping must be a character vector of length 2")
        }

        obsm_slot <- obsm[[1]]
        obsm_key <- obsm[[2]]

        adata$obsm[[obsm_name]] <- .from_SCE_process_obsm_reduction(sce, obsm_slot, obsm_key)
      }

      # fetch varm
      # trackstatus: class=SingleCellExperiment, feature=set_varm, status=wip
      for (i in seq_along(varm_mapping)) {
        varm <- varm_mapping[[i]]
        varm_name <- names(varm_mapping)[[i]]

        if (!is.character(varm) || length(varm) != 2) {
          stop("each varm_mapping must be a character vector of length 2")
        }

        varm_slot <- varm[[1]]
        varm_key <- varm[[2]]

        if (varm_slot != "reducedDim") {
          stop("varm_slot must be 'reducedDims'")
        }

        if (! inherits(reducedDims(sce)[[varm_key]], "LinearEmbeddingMatrix")) {
          stop("reducedDim must be a LinearEmbeddingMatrix")
        }

        adata$varm[[varm_name]] <- featureLoadings(reducedDim(sce, varm_key))

      }

      # fetch obsp
      # trackstatus: class=SingleCellExperiment, feature=set_obsp, status=wip
      obsp <- .from_SCE_process_pairs(sce, obsp_mapping, colPairs, asSparse = TRUE)
      
      adata$obsp <- .from_SCE_process_pairs(sce, obsp_mapping, colPairs, asSparse = TRUE)

      # fetch varp
      # trackstatus: class=SingleCellExperiment, feature=set_varp, status=wip
      adata$varp <- .from_SCE_process_pairs(sce, varp_mapping, rowPairs, asSparse = TRUE)

      # fetch uns
      # trackstatus: class=SingleCellExperiment, feature=set_uns, status=wip
      adata$uns <- .from_SCE_process_simple_mapping(sce, uns_mapping, metadata)

      adata
    },
    error = function(e) {
      if (output_class == "HDF5AnnData") {
        on.exit(cleanup_HDF5AnnData(adata))
      }
      stop(e)
    }
  )
}

.from_SCE_guess_all <- function(sce, slot) {
  mapping <- names(slot(sce))
  names(mapping) <- names(slot(sce))

  mapping
}

.from_SCE_guess_layers <- function(sce, x_mapping) {
  layers_mapping <- list()

  for (assay_name in names(assays(sce))) {
    if (!is.null(x_mapping) && assay_name != x_mapping) {
      layers_mapping[[assay_name]] <- assay_name
    }
  }

  layers_mapping
}

.from_SCE_guess_obsm <- function(sce) {
  obsm_mapping <- list()

  for (reduction_name in names(SingleCellExperiment::reducedDims(sce))) {
    dest_name <- paste0("X_", reduction_name)
    obsm_mapping[[dest_name]] <- c("reducedDim", reduction_name)
  }

  obsm_mapping
}

.from_SCE_guess_varm <- function(sce) {
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

.from_SCE_convert <- function(object) {
  if (inherits(object, "DataFrame")) {
    as.data.frame(object)
  } else if (inherits(object, "SimpleList")) {
    as.list(object)
  } else if (inherits(object, "matrix") || inherits(object, "sparseMatrix")) {
    m <- t(object)
    if (inherits(m, "denseMatrix")) {
      m <- as.matrix(m)
    }
    m
  } else {
    object
  }
}

.to_anndata_matrix <- function(mat) {
  m <- t(mat)
  if (inherits(m, "denseMatrix")) {
    m <- as.matrix(m)
  }
  m
}

.from_SCE_process_simple_mapping <- function(sce, mapping, slot, ...) {
  mapped <- NULL
  # check if mapping contains all columns of slot
  if (length(setdiff(names(slot(sce)), names(mapping))) == 0) {
    mapped <- slot(sce)
  } else {
    mapped <- lapply(seq_along(mapping), function(i) {

      .from_SCE_convert(slot(sce, ...)[[mapping[[i]]]])
    })
    names(mapped) <- names(mapping)
  }
  
  .from_SCE_convert(mapped)
}

.from_SCE_process_pairs <- function(sce, mapping, slot, asSparse = TRUE) {
  # check if slot is colPairs or rowPairs

  pairs <- NULL
  # check if mapping contains all columns of slot
  if (length(setdiff(names(slot(sce)), names(mapping))) == 0) {
    pairs <- slot(sce, asSparse=asSparse)
  } else {
    pairs <- lapply(seq_along(mapping), function(i) {
      .from_SCE_convert(slot(sce, asSparse=asSparse)[[mapping[[i]]]])
    })
    names(pairs) <- names(mapping)
  }
  
  .from_SCE_convert(pairs)
}

# .from_SCE_process_layers_mapping <- function(sce, mapping, slot, toDF = TRUE) {
#   layers <- list()
#   for (i in seq_along(mapping)) {
#     layer <- mapping[[i]]
#     layer_name <- names(mapping)[[i]]

#     layers[[layer_name]] <- .to_anndata_matrix(assay(sce, layer))
#   }
#   if (toDF == TRUE) {
#     layers <- as.data.frame(layers)
#   }
#   layers
# }

.from_SCE_process_obsm_reduction <- function(sce, slot, key) {
  if (slot != "reducedDim") {
    stop("slot must be 'reducedDim'")
  }
  
  reduction <- reducedDim(sce, key)
  if (inherits(reduction, "LinearEmbeddingMatrix")) {
    sampleFactors(reduction)
  } else {
    reduction
  }
}
