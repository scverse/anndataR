#' S3 Methods for AbstractAnnData Objects
#'
#' These S3 methods provide standard R interfaces for AbstractAnnData objects,
#' making them behave like native R objects with familiar syntax.
#'
#' @name AbstractAnnData-s3methods
#' @param x An AbstractAnnData object
#' @param i Row indices (observations). Can be numeric, logical, or character.
#' @param j Column indices (variables). Can be numeric, logical, or character.
#' @param drop Ignored (for compatibility with generic)
#' @param ... Additional arguments passed to methods
#'
#' @details
#' The following S3 methods are available:
#'
#'   * `dim(x)`: Get dimensions (n_obs, n_vars)
#'   * `nrow(x)`: Get number of observations
#'   * `ncol(x)`: Get number of variables
#'   * `rownames(x)`: Get observation names
#'   * `colnames(x)`: Get variable names
#'   * `x[i, j]`: Subset observations and/or variables
#'
#' @return
#'
#'   * `dim`: Numeric vector of length 2 (n_obs, n_vars)
#'   * `nrow`, `ncol`: Integer count
#'   * `rownames`, `colnames`: Character vector
#'   * `[`: A AnnDataView object with the specified subset
#'
#' @examples
#' \dontrun{
#' # Create example data
#' ad <- generate_dataset(n_obs = 100, n_vars = 50)
#'
#' # Standard R methods work
#' dim(ad)
#' nrow(ad)
#' ncol(ad)
#' rownames(ad)
#' colnames(ad)
#'
#' # Subsetting creates AnnDataView
#' subset_ad <- ad[1:10, 1:5]
#' subset_ad <- ad[c(TRUE, FALSE), ]  # logical subsetting
#' subset_ad <- ad[c("obs_1", "obs_2"), c("var_1", "var_2")]  # name subsetting
#' }
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
#' @method rownames AbstractAnnData
#' @export
rownames.AbstractAnnData <- function(x) {
  x$obs_names
}

#' @rdname AbstractAnnData-s3methods
#' @method colnames AbstractAnnData
#' @export
colnames.AbstractAnnData <- function(x) {
  x$var_names
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
  AnnDataView$new(x, i, j)
}
