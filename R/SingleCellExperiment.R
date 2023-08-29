#' @rdname SingleCellExperiment
#'
#' @title Convert Between AnnData and SingleCellExperiment
#'
#' @description `to_SingleCellExperiment()` converts an AnnData object
#'   to a SingleCellExperiment.
#'
#' @param object an AnnData object, e.g., InMemoryAnnData
#'
#' @return `to_SingleCellExperiment()` returns a SingleCellExperiment
#'   representing the content of `object`.
#'
#' @examples
#' if (interactive()) {
#'   ## useful when interacting with the SingleCellExperiment !
#'   library(SingleCellExperiment)
#' }
#' ad <- InMemoryAnnData$new(
#'   X = matrix(1:15, 3L, 5L),
#'   layers = list(
#'     A = matrix(15:1, 3L, 5L),
#'     B = matrix(letters[1:15], 3L, 5L)
#'   ),
#'   obs = data.frame(cell = 1:3),
#'   var = data.frame(gene = 1:5),
#'   obs_names = LETTERS[1:3],
#'   var_names = letters[1:5]
#' )
#'
#' ## construct a SingleCellExperiment from an AnnData object
#' sce <- to_SingleCellExperiment(ad)
#' sce
#'
#' @export
to_SingleCellExperiment <- function(object) { # nolint
  stopifnot(
    inherits(object, "AbstractAnnData")
  )

  # trackstatus: class=SingleCellExperiment, feature=get_X, status=done
  ## FIXME: name of 'X' from metadata[["X_name"]]
  x_name <- "X"
  assay_names <- as.character(c(
    if (!is.null(object$X)) x_name,
    object$layers_keys()
  ))

  # trackstatus: class=SingleCellExperiment, feature=get_layers, status=done
  assays <- vector("list", length(assay_names))
  names(assays) <- assay_names
  for (assay in assay_names) {
    value <- if (identical(assay, x_name)) {
      object$X
    } else {
      object$layers[[assay]]
    }
    ## FIXME: is transposition robust & efficient here?
    assays[[assay]] <- t(value)
  }

  # construct colData
  # trackstatus: class=SingleCellExperiment, feature=get_obs, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_obs_names, status=done
  col_data <- S4Vectors::DataFrame(
    object$obs,
    row.names = object$obs_names
  )

  # construct rowData
  # trackstatus: class=SingleCellExperiment, feature=get_var, status=done
  # trackstatus: class=SingleCellExperiment, feature=get_var_names, status=done
  row_data <- S4Vectors::DataFrame(
    object$var,
    row.names = object$var_names
  )

  # construct output object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = col_data,
    rowData = row_data,
    metadata = list(),
    ## FIXME: assign object$uns to metadata
    checkDimnames = TRUE
  )

  ## reducedDims

  ## rowPairs

  sce
}

#' @rdname SingleCellExperiment
#'
#' @description `from_SingleCellExperiment()` converts a
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
#' from_SingleCellExperiment(sce, "InMemory")
#'
#' @export
from_SingleCellExperiment <- function(sce, output_class = c("InMemory", "HDF5AnnData"), ...) { # nolint
  stopifnot(
    inherits(sce, "SingleCellExperiment")
  )

  # fetch generator
  generator <- get_generator(output_class)

  # trackstatus: class=SingleCellExperiment, feature=set_obs, status=done
  obs <- SummarizedExperiment::as.data.frame(
    SummarizedExperiment::colData(sce)
  )
  rownames(obs) <- NULL

  # trackstatus: class=SingleCellExperiment, feature=set_var, status=done
  var <- SummarizedExperiment::as.data.frame(
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
    t
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
