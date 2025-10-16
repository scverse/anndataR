#' S3 Methods for AbstractAnnData Objects
#'
#' These S3 methods provide standard R interfaces for AbstractAnnData objects,
#' making them behave like native R objects with familiar syntax.
#'
#' @details
#' **Subsetting behaviour**: The `[` method supports logical, integer, and character
#' subsetting for both observations (rows) and variables (columns). However, unlike
#' standard R behaviour:
#'
#' - Logical vectors are **not recycled** and must have the exact same length as
#'   the dimension being subset
#' - **Negative indices are not supported** (R's "exclude these" syntax)
#'
#' These design choices ensure clear and predictable subsetting behaviour for
#' biological data matrices, avoiding potential confusion from accidental recycling
#' or exclusion patterns.
#'
#' @name AbstractAnnData-s3methods
#' @param x An AbstractAnnData object
#' @param value For `dimnames<-`: A list of two character vectors (obs_names, var_names)
#' @param i Row indices (observations). Can be numeric, logical, or character.
#' @param j Column indices (variables). Can be numeric, logical, or character.
#' @param drop Ignored (for compatibility with generic)
#' @param ... Additional arguments passed to methods
#'
#' @details
#' The following S3 methods are available:
#'
#'   * `dim(x)`: Get dimensions (n_obs, n_vars), equivalent to `x$shape()`
#'   * `nrow(x)`: Get number of observations, equivalent to `x$n_obs()`
#'   * `ncol(x)`: Get number of variables, equivalent to `x$n_vars()`
#'   * `dimnames(x)`: Get dimension names, (obs_names, var_names)
#'   * `rownames(x)`: Get observation names, equivalent to `x$obs_names()`
#'   * `colnames(x)`: Get variable names, equivalent to `x$var_names()`
#'   * `dimnames(x) <- value`: Set dimension names
#'   * `rownames(x) <- value`: Set observation names, equivalent to `x$obs_names() <- ...`
#'   * `colnames(x) <- value`: Set variable names, equivalent to `x$var_names() <- ...`
#'   * `x[i, j]`: Subset observations and/or variables
#'
#' @return
#'
#'   * `dim`: Numeric vector of length 2 (n_obs, n_vars)
#'   * `nrow`, `ncol`: Integer count
#'   * `dimnames`: List with obs_names and var_names
#'   * `rownames`, `colnames`: Character vector
#'   * `dimnames<-`, `rownames<-`, `colnames<-`: The modified object (invisibly)
#'   * `[`: A AnnDataView object with the specified subset
#'
#' @examples
#' # Create example data
#' ad <- generate_dataset(n_obs = 100, n_vars = 50, format = "AnnData")
#'
#' # Standard R methods work
#' dim(ad)
#' nrow(ad)
#' ncol(ad)
#' dimnames(ad)
#' rownames(ad)
#' colnames(ad)
#'
#' # Set names using dimnames
#' dimnames(ad) <- list(
#'   paste0("cell_", 1:nrow(ad)),
#'   paste0("gene_", 1:ncol(ad))
#' )
#'
#' # Or set names individually (uses dimnames<- internally)
#' rownames(ad) <- paste0("cell_", 1:nrow(ad))
#' colnames(ad) <- paste0("gene_", 1:ncol(ad))
#'
#' # Subsetting creates AnnDataView
#' subset_ad <- ad[1:10, 1:5]
#' subset_ad <- ad[rep(c(TRUE, FALSE), length.out = nrow(ad)), ]  # logical subsetting (no recycling)
#' subset_ad <- ad[c("cell_1", "cell_2"), c("gene_1", "gene_2")]  # name subsetting
NULL

#' @rdname AbstractAnnData-s3methods
#' @method dim AbstractAnnData
#' @export
dim.AbstractAnnData <- function(x) {
  x$shape()
}

#' @rdname AbstractAnnData-s3methods
#' @method nrow AbstractAnnData
#' @export
nrow.AbstractAnnData <- function(x) {
  x$n_obs()
}

#' @rdname AbstractAnnData-s3methods
#' @method ncol AbstractAnnData
#' @export
ncol.AbstractAnnData <- function(x) {
  x$n_vars()
}

#' @rdname AbstractAnnData-s3methods
#' @method dimnames AbstractAnnData
#' @export
dimnames.AbstractAnnData <- function(x) {
  list(x$obs_names, x$var_names)
}

#' @rdname AbstractAnnData-s3methods
#' @method dimnames<- AbstractAnnData
#' @export
`dimnames<-.AbstractAnnData` <- function(x, value) {
  if (is.null(value)) {
    # When dimnames(x) <- NULL, generate default sequential names
    # This mimics Python AnnData behaviour but uses 1-based indexing (R convention)
    # instead of 0-based indexing (Python convention)
    x$obs_names <- as.character(seq_len(x$n_obs()))
    x$var_names <- as.character(seq_len(x$n_vars()))
  } else {
    if (is.null(value[[1]])) {
      x$obs_names <- as.character(seq_len(x$n_obs()))
    } else {
      x$obs_names <- value[[1]]
    }

    if (is.null(value[[2]])) {
      x$var_names <- as.character(seq_len(x$n_vars()))
    } else {
      x$var_names <- value[[2]]
    }
  }
  x
}

#' @rdname AbstractAnnData-s3methods
#' @method [ AbstractAnnData
#' @export
`[.AbstractAnnData` <- function(x, i, j, drop = TRUE, ...) {
  if (inherits(x, "AnnDataView")) {
    # If x is already a view, we need to update the view with new indices
    return(x$subset(i, j))
  }

  # Create AnnDataView with both subsets at once
  AnnDataView$new(x, i, j) # nolint: object_usage_linter
}
