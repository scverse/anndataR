#' Convert to/from SingleCellExperiment objects
#'
#' Conversion between AnnData and `SingleCellExperiment` objects.
#'
#' @seealso
#' [`write_h5ad()`] and [`read_h5ad()`] for directly interacting with H5AD files.
#'
#' @name SingleCellExperiment-Conversion
#' @rdname SingleCellExperiment-Conversion
NULL

#' @rdname SingleCellExperiment-Conversion
#'
#' @param obj an AnnData object, e.g., InMemoryAnnData
#'
#' @return `to_SingleCellExperiment()` returns a `SingleCellExperiment`
#'   representing the content of `object`
to_SingleCellExperiment <- function(obj) { # nolint
  stopifnot(
    inherits(obj, "AbstractAnnData")
  )

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  ## FIXME: name of 'X' from metadata[["X_name"]]
  x_name <- "X"
  assay_names <- as.character(c(
    if (!is.null(obj$X)) x_name,
    obj$layers_keys()
  ))

  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  assays <- vector("list", length(assay_names))
  names(assays) <- assay_names
  for (assay in assay_names) {
    value <- if (identical(assay, x_name)) {
      obj$X
    } else {
      obj$layers[[assay]]
    }
    ## FIXME: is transposition robust & efficient here?
    assays[[assay]] <- t(value)
  }

  # construct colData
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- S4Vectors::DataFrame(
    obj$obs,
    row.names = obj$obs_names
  )

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- S4Vectors::DataFrame(
    obj$var,
    row.names = obj$var_names
  )

  # construct output object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = col_data,
    rowData = row_data,
    metadata = list(),
    ## FIXME: assign obj$uns to metadata
    checkDimnames = TRUE
  )

  ## reducedDims

  ## rowPairs

  sce
}

#' @rdname SingleCellExperiment-Conversion
#'
#' @param sce An object inheriting from SingleCellExperiment
#' @param output_class Name of the AnnData class. Must be one of `"HDF5AnnData"`
#'   or `"InMemoryAnnData"`.
#' @param ... Additional arguments passed to the object constructor
#'
#' @return `from_SingleCellExperiment()` returns an AnnData object
#'   (e.g., InMemoryAnnData) representing the content of `sce`
from_SingleCellExperiment <- function(sce, # nolint
                                      output_class = c("InMemory", "HDF5AnnData"),
                                      ...) {

  stopifnot(
    inherits(sce, "SingleCellExperiment")
  )

  # fetch generator
  generator <- get_anndata_constructor(output_class)

  # trackstatus: class=SingleCellExperiment, feature=set_obs, status=done
  obs <- as.data.frame(
    SummarizedExperiment::colData(sce)
  )
  rownames(obs) <- NULL

  # trackstatus: class=SingleCellExperiment, feature=set_var, status=done
  var <- as.data.frame(
    SummarizedExperiment::rowData(sce)
  )
  rownames(var) <- NULL

  # trackstatus: class=SingleCellExperiment, feature=set_obs_names, status=done
  obs_names <- colnames(sce)
  if (is.null(obs_names)) {
    warning(wrap_message("colnames(sce) should not be NULL"))
    obs_names <- as.character(seq_len(nrow(obs)))
  }

  # trackstatus: class=SingleCellExperiment, feature=set_var_names, status=done
  var_names <- rownames(sce)
  if (is.null(var_names)) {
    warning(wrap_message("rownames(sce) should not be NULL"))
    var_names <- as.character(seq_len(nrow(var)))
  }

  # trackstatus: class=SingleCellExperiment, feature=set_X, status=done
  # trackstatus: class=SingleCellExperiment, feature=set_layers, status=done
  x_and_layers <- lapply(
    SummarizedExperiment::assays(sce, withDimnames = FALSE),
    function(mat) {
      m <- t(mat)
      # nolint start
      # WORKAROUND: convert denseMatrix to matrix, because otherwise:
      # - Could not write element '/layers/integer_dense' of type 'dgeMatrix':
      #   no applicable method for 'h5writeDataset' applied to an object of class "c('dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'compMa
      # - Could not write element '/layers/integer_dense_with_nas' of type 'dgeMatrix':
      #   no applicable method for 'h5writeDataset' applied to an object of class "c('dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'compMatrix', 'Matrix', 'replValueSp')"
      # nolint end
      if (inherits(m, "denseMatrix")) {
        m <- as.matrix(m)
      }
      m
    }
  )
  if (length(x_and_layers) == 0L) {
    x <- NULL
    layers <- list()
    names(layers) <- character()
  } else {
    x <- x_and_layers[[1]]
    layers <- x_and_layers[-1]
  }

  generator$new(
    X = x,
    obs = obs,
    var = var,
    obs_names = obs_names,
    var_names = var_names,
    layers = layers,
    ...
  )
}
